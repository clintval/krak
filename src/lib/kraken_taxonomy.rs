//! Parse Kraken taxonomy from database files.
//!
//! Loads the Kraken taxonomy tree from either a Kraken v1 or Kraken v2
//! database, returning the same [`KrakenReportEntry`] structure in both
//! cases so downstream consumers (e.g. embedding the taxonomy in an
//! annotated BAM) work uniformly across versions and avoid the need for a
//! per-run `kraken2 --report` file.
//!
//! - **Kraken v1** ships the NCBI taxonomy as flat files under the database
//!   directory (`taxonomy/nodes.dmp` and `taxonomy/names.dmp`). See
//!   [`read_taxonomy_dmp`].
//! - **Kraken v2** stores the pruned NCBI taxonomy in a single binary
//!   `taxo.k2d` file inside the database directory. See [`read_taxo_k2d`].
//!
//! # Kraken v2 `taxo.k2d` binary layout
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────┐
//! │ Header (32 bytes)                                               │
//! │   [0..8]   magic: "K2TAXDAT" (8 ASCII bytes, no null term.)     │
//! │   [8..16]  node_count:     u64 LE                               │
//! │   [16..24] name_data_len:  u64 LE                               │
//! │   [24..32] rank_data_len:  u64 LE                               │
//! ├─────────────────────────────────────────────────────────────────┤
//! │ Node array (node_count x 56 bytes)                              │
//! │   Node 0: invalid/placeholder sentinel (all zeros).             │
//! │   Node 1: root of the taxonomy tree.                            │
//! │   Each node: 7 x u64 LE; see TaxoNode.                         │
//! ├─────────────────────────────────────────────────────────────────┤
//! │ Name data (name_data_len bytes)                                 │
//! │   Null-terminated C-strings, concatenated.                      │
//! ├─────────────────────────────────────────────────────────────────┤
//! │ Rank data (rank_data_len bytes)                                 │
//! │   Null-terminated C-strings, concatenated (deduplicated).       │
//! └─────────────────────────────────────────────────────────────────┘
//! ```

use std::io::Read;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::kraken_report::{KrakenReportEntry, SPACES_PER_DEPTH};

const MAGIC: &[u8; 8] = b"K2TAXDAT";
const NODE_BYTE_LEN: usize = 56; // 7 of u64

/// A single node as stored in the `taxo.k2d` node array.
#[derive(Debug, Clone, Copy)]
struct TaxoNode {
    _parent_id: u64,
    first_child: u64,
    child_count: u64,
    name_offset: u64,
    rank_offset: u64,
    external_id: u64,
    // godparent_id  // is reserved/unused; so let's ignore it too
}

impl TaxoNode {
    fn from_bytes(buf: &[u8; NODE_BYTE_LEN]) -> Self {
        let mut off = 0;
        let mut next = || {
            let v = u64::from_le_bytes(buf[off..off + 8].try_into().unwrap());
            off += 8;
            v
        };
        let p = next();
        let fc = next();
        let cc = next();
        let no = next();
        let ro = next();
        let ei = next();
        TaxoNode {
            _parent_id: p,
            first_child: fc,
            child_count: cc,
            name_offset: no,
            rank_offset: ro,
            external_id: ei,
        }
    }
}

/// Read a null-terminated C-string from `data` at `offset`.
///
/// Errors if `offset` is out of bounds or string bytes are not valid UTF-8.
fn cstr_at(data: &[u8], offset: usize) -> Result<&str> {
    let rest = data.get(offset..).ok_or_else(|| {
        anyhow::anyhow!(
            "string offset {offset} is out of bounds (data length {})",
            data.len()
        )
    })?;
    let len = rest.iter().position(|&b| b == 0).unwrap_or(rest.len());
    std::str::from_utf8(&rest[..len])
        .with_context(|| format!("string at offset {offset} contains invalid UTF-8"))
}

/// Map a Kraken NCBI rank name to its abbreviated rank code.
fn rank_code(rank: &str) -> &'static str {
    match rank {
        "root" | "Root" => "R",
        "domain" | "superkingdom" => "D",
        "kingdom" => "K",
        "subkingdom" => "K1",
        "phylum" => "P",
        "subphylum" => "P1",
        "class" => "C",
        "subclass" => "C1",
        "infraclass" => "C2",
        "order" => "O",
        "superorder" => "O1",
        "suborder" => "O2",
        "infraorder" => "O3",
        "parvorder" => "O4",
        "family" => "F",
        "superfamily" => "F1",
        "subfamily" => "F2",
        "tribe" => "F3",
        "subtribe" => "F4",
        "genus" => "G",
        "subgenus" => "G1",
        "species group" => "G2",
        "species subgroup" => "G3",
        "species" => "S",
        "subspecies" => "S1",
        "varietas" | "variety" => "S2",
        "forma" | "form" => "S3",
        other => {
            log::debug!("unmapped Kraken rank: {other:?}");
            "--"
        }
    }
}

/// Read an exact u64 from the buffer
fn read_u64_le(f: &mut impl Read) -> Result<u64> {
    let mut buf = [0u8; 8];
    f.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

/// Parse a `taxo.k2d` file and return its taxonomy as a list of
/// [`KrakenReportEntry`] values in pre-order DFS order, with indent values
/// matching the `krak:report:` BAM header embedding convention (2 spaces per
/// depth level).
///
/// Nodes with `external_id == 0` (the placeholder sentinel plus any internal
/// Kraken-only nodes) are skipped; their children are promoted to the current
/// depth so the visible tree structure is preserved.
pub fn read_taxo_k2d(path: &Path) -> Result<Vec<KrakenReportEntry>> {
    let mut handle =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;

    let mut magic = [0u8; 8];
    handle
        .read_exact(&mut magic)
        .context("failed to read taxo.k2d magic")?;
    if &magic != MAGIC {
        bail!(
            "{}: not a valid taxo.k2d file (bad magic bytes)",
            path.display()
        );
    }

    let node_count = read_u64_le(&mut handle).context("failed to read node_count")? as usize;
    let name_data_len = read_u64_le(&mut handle).context("failed to read name_data_len")? as usize;
    let rank_data_len = read_u64_le(&mut handle).context("failed to read rank_data_len")? as usize;

    // Reject declared sizes that cannot fit in the file before allocating, so a
    // corrupt/hostile header can't trigger a giant allocation that aborts the
    // process instead of producing a clean error. The `read_exact` calls below
    // re-validate the exact byte counts.
    let file_len = handle
        .metadata()
        .with_context(|| format!("failed to stat {}", path.display()))?
        .len();
    const HEADER_LEN: u64 = 8 + 8 * 3; // magic + node_count + name/rank lengths
    let body_len = file_len.saturating_sub(HEADER_LEN);
    let max_nodes = body_len / NODE_BYTE_LEN as u64;
    if node_count as u64 > max_nodes {
        bail!(
            "{}: taxo.k2d declares node_count={node_count} but the file holds room for at \
             most {max_nodes} nodes; file is truncated or corrupt",
            path.display()
        );
    }
    if name_data_len as u64 > body_len || rank_data_len as u64 > body_len {
        bail!(
            "{}: taxo.k2d declares name_data_len={name_data_len}/rank_data_len={rank_data_len} \
             exceeding the {body_len}-byte body; file is truncated or corrupt",
            path.display()
        );
    }

    let mut nodes = Vec::with_capacity(node_count);
    for i in 0..node_count {
        let mut buf = [0u8; NODE_BYTE_LEN];
        handle
            .read_exact(&mut buf)
            .with_context(|| format!("failed to read node {i}"))?;
        nodes.push(TaxoNode::from_bytes(&buf));
    }

    let mut name_data = vec![0u8; name_data_len];
    handle
        .read_exact(&mut name_data)
        .context("failed to read name_data")?;
    let mut rank_data = vec![0u8; rank_data_len];
    handle
        .read_exact(&mut rank_data)
        .context("failed to read rank_data")?;

    if node_count < 2 {
        bail!(
            "{}: taxo.k2d has node_count={node_count}; expected at least a placeholder and a root",
            path.display()
        );
    }

    // Validate that every node's child range is in-bounds before traversal.
    for (i, node) in nodes.iter().enumerate() {
        let first = node.first_child as usize;
        let count = node.child_count as usize;
        if count > 0 {
            let end = first.checked_add(count).ok_or_else(|| {
                anyhow::anyhow!("node {i}: first_child + child_count overflows usize")
            })?;
            if end > node_count {
                bail!("node {i}: child range [{first}, {end}) exceeds node_count {node_count}");
            }
        }
    }

    // Pre-order the DFS from node 1 (root), stack entries: (node_index, depth).
    // Nodes with external_id == 0 are skipped; their children are emitted at
    // the same depth as the skipped node (promotion).
    let mut entries: Vec<KrakenReportEntry> = Vec::new();
    let mut stack: Vec<(usize, usize)> = vec![(1, 0)];
    let mut visited: Vec<bool> = vec![false; node_count];

    while let Some((idx, depth)) = stack.pop() {
        if idx == 0 || idx >= node_count {
            continue;
        }
        if visited[idx] {
            bail!("cycle detected in taxo.k2d at node index {idx}");
        }
        visited[idx] = true;
        let node = nodes[idx];
        let first = node.first_child as usize;
        let count = node.child_count as usize;

        if node.external_id == 0 {
            // Invisible node: promote children to the current depth.
            for i in (0..count).rev() {
                stack.push((first + i, depth));
            }
            continue;
        }

        entries.push(KrakenReportEntry {
            pct_fragments: 0.0,
            num_fragments_clade: 0,
            num_fragments_direct: 0,
            rank_code: rank_code(cstr_at(&rank_data, node.rank_offset as usize)?).to_owned(),
            taxon_id: u32::try_from(node.external_id).map_err(|_| {
                anyhow::anyhow!(
                    "taxo.k2d node {idx}: external_id {} exceeds u32 range",
                    node.external_id
                )
            })?,
            name: cstr_at(&name_data, node.name_offset as usize)?.to_owned(),
            indent: depth * SPACES_PER_DEPTH,
            minimizer_count: None,
            distinct_minimizer_count: None,
        });

        // Push children in reverse order so the first child is processed first.
        for i in (0..count).rev() {
            stack.push((first + i, depth + 1));
        }
    }

    Ok(entries)
}

/// Parse NCBI taxonomy flat files (`taxonomy/nodes.dmp` and `taxonomy/names.dmp`)
/// from a Kraken v1 database directory and return a [`Vec<KrakenReportEntry>`] in
/// pre-order DFS order; the same structure produced by [`read_taxo_k2d`].
///
/// Each line in `nodes.dmp` is `|`-delimited with surrounding whitespace; only
/// the first three fields (`taxon_id`, `parent_id`, `rank`) are used. Each line
/// in `names.dmp` uses the same delimiter; only entries whose fourth field is
/// `"scientific name"` are kept.
pub fn read_taxonomy_dmp(db_path: &Path) -> Result<Vec<KrakenReportEntry>> {
    use ahash::{AHashMap, AHashSet};
    use std::io::BufRead as _;

    let nodes_path = db_path.join("taxonomy").join("nodes.dmp");
    let names_path = db_path.join("taxonomy").join("names.dmp");

    let mut names: AHashMap<u32, String> = AHashMap::new();
    let names_file = std::fs::File::open(&names_path)
        .with_context(|| format!("failed to open {}", names_path.display()))?;
    for line in std::io::BufReader::new(names_file).lines() {
        let line = line.context("failed to read names.dmp")?;
        let parts: Vec<&str> = line.split('|').map(str::trim).collect();
        if parts.len() >= 4 && parts[3] == "scientific name" {
            if let Ok(id) = parts[0].parse::<u32>() {
                names.insert(id, parts[1].to_owned());
            }
        }
    }

    let mut rank_of: AHashMap<u32, String> = AHashMap::new();
    let mut children_of: AHashMap<u32, Vec<u32>> = AHashMap::new();

    let nodes_file = std::fs::File::open(&nodes_path)
        .with_context(|| format!("failed to open {}", nodes_path.display()))?;

    for line in std::io::BufReader::new(nodes_file).lines() {
        let line = line.context("failed to read nodes.dmp")?;
        let parts: Vec<&str> = line.split('|').map(str::trim).collect();
        if parts.len() < 3 {
            continue;
        }
        let (Ok(taxon_id), Ok(parent_id)) = (parts[0].parse::<u32>(), parts[1].parse::<u32>())
        else {
            continue;
        };
        rank_of.insert(taxon_id, parts[2].to_owned());
        // The root (taxon 1) has parent 1; skip the self-referential edge.
        if taxon_id != parent_id {
            children_of.entry(parent_id).or_default().push(taxon_id);
        }
    }

    // Sort children for deterministic pre-order output.
    for children in children_of.values_mut() {
        children.sort_unstable();
    }

    if !rank_of.contains_key(&1) {
        bail!(
            "{}: nodes.dmp has no entry for root (taxon 1)",
            nodes_path.display()
        );
    }

    // Pre-order the DFS from root (taxon 1)
    let mut entries: Vec<KrakenReportEntry> = Vec::new();
    let mut stack: Vec<(u32, usize)> = vec![(1, 0)];
    let mut visited: AHashSet<u32> = AHashSet::new();

    while let Some((id, depth)) = stack.pop() {
        if !visited.insert(id) {
            bail!("cycle detected in nodes.dmp while visiting taxon {id}");
        }
        let name = names.get(&id).cloned().unwrap_or_default();
        let rank = rank_of.get(&id).map(String::as_str).unwrap_or("no rank");
        entries.push(KrakenReportEntry {
            pct_fragments: 0.0,
            num_fragments_clade: 0,
            num_fragments_direct: 0,
            rank_code: rank_code(rank).to_owned(),
            taxon_id: id,
            name,
            indent: depth * SPACES_PER_DEPTH,
            minimizer_count: None,
            distinct_minimizer_count: None,
        });
        if let Some(children) = children_of.get(&id) {
            // Reverse so the smallest taxon ID is processed first.
            for &child in children.iter().rev() {
                stack.push((child, depth + 1));
            }
        }
    }

    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn minimal_taxo_k2d_bytes() -> Vec<u8> {
        // Name data:
        //   offset  0: "root\0"              (5 bytes)
        //   offset  5: "Rodentia\0"          (9 bytes)
        //   offset 14: "Rattus norvegicus\0" (18 bytes)
        let name_data: &[u8] = b"root\0Rodentia\0Rattus norvegicus\0";
        // Rank data:
        //   offset  0: "no rank\0"  (8 bytes)
        //   offset  8: "order\0"    (6 bytes)
        //   offset 14: "species\0"  (8 bytes)
        let rank_data: &[u8] = b"no rank\0order\0species\0";

        // (parent_id, first_child, child_count, name_offset, rank_offset, external_id, godparent_id)
        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0),       // Node 0: placeholder
            (0, 2, 1, 0, 0, 1, 0),       // Node 1: root (ext_id=1)
            (1, 3, 1, 5, 8, 9989, 0),    // Node 2: Rodentia (ext_id=9989)
            (2, 0, 0, 14, 14, 10116, 0), // Node 3: Rattus norvegicus (ext_id=10116)
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);
        buf
    }

    fn write_tmp(bytes: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(bytes).unwrap();
        f
    }

    #[test]
    fn test_produces_entries_in_dfs_order() {
        let tmp = write_tmp(&minimal_taxo_k2d_bytes());
        let entries = read_taxo_k2d(tmp.path()).unwrap();

        assert_eq!(entries.len(), 3, "expected root + Rodentia + Rattus");

        assert_eq!(entries[0].taxon_id, 1);
        assert_eq!(entries[0].name, "root");
        assert_eq!(entries[0].indent, 0);
        assert_eq!(entries[0].rank_code, "--"); // "no rank" -> "--"

        assert_eq!(entries[1].taxon_id, 9989);
        assert_eq!(entries[1].name, "Rodentia");
        assert_eq!(entries[1].indent, 2);
        assert_eq!(entries[1].rank_code, "O"); // "order" -> "O"

        assert_eq!(entries[2].taxon_id, 10116);
        assert_eq!(entries[2].name, "Rattus norvegicus");
        assert_eq!(entries[2].indent, 4);
        assert_eq!(entries[2].rank_code, "S"); // "species" -> "S"
    }

    #[test]
    fn test_builds_valid_taxonomy_tree() {
        use crate::kraken_report::KrakenTaxonomyTree;

        let tmp = write_tmp(&minimal_taxo_k2d_bytes());
        let entries = read_taxo_k2d(tmp.path()).unwrap();
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();

        assert!(tree.contains(1));
        assert!(tree.contains(9989));
        assert!(tree.contains(10116));
        assert_eq!(tree.parent_of(9989), Some(1));
        assert_eq!(tree.parent_of(10116), Some(9989));
        assert!(tree.descendants_of(9989).contains(&10116));
    }

    #[test]
    fn test_bad_magic_errors() {
        let mut bytes = minimal_taxo_k2d_bytes();
        bytes[0] = b'X';
        let tmp = write_tmp(&bytes);
        assert!(read_taxo_k2d(tmp.path()).is_err());
    }

    #[test]
    fn test_taxo_k2d_too_small_errors() {
        // node_count must be >= 2 (placeholder + root); a single-node file
        // is structurally invalid and must produce a clear error.
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&1u64.to_le_bytes()); // node_count = 1
        buf.extend_from_slice(&0u64.to_le_bytes()); // name_data_len
        buf.extend_from_slice(&0u64.to_le_bytes()); // rank_data_len
                                                    // one TaxoNode of zeros
        buf.extend_from_slice(&[0u8; NODE_BYTE_LEN]);
        let tmp = write_tmp(&buf);
        let err = read_taxo_k2d(tmp.path()).unwrap_err();
        assert!(format!("{err:#}").contains("expected at least"));
    }

    #[test]
    fn test_taxo_k2d_implausible_node_count_errors_cleanly() {
        // A corrupt/hostile header declaring an enormous node_count must
        // produce a clean error rather than abort the process inside a giant
        // `Vec::with_capacity`.
        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&u64::MAX.to_le_bytes()); // node_count = u64::MAX
        buf.extend_from_slice(&0u64.to_le_bytes()); // name_data_len
        buf.extend_from_slice(&0u64.to_le_bytes()); // rank_data_len
        let tmp = write_tmp(&buf);
        let err = read_taxo_k2d(tmp.path()).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("node_count") || msg.contains("truncated") || msg.contains("corrupt"),
            "expected a size-sanity error, got: {msg}"
        );
    }

    #[test]
    fn test_taxo_k2d_external_id_exceeding_u32_errors() {
        // An external_id above u32::MAX must error, not silently truncate to a
        // wrong (or zero) taxon id.
        let name_data: &[u8] = b"root\0";
        let rank_data: &[u8] = b"no rank\0";
        let big: u64 = u32::MAX as u64 + 1; // 4_294_967_296
        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0),   // placeholder
            (0, 0, 0, 0, 0, big, 0), // root with out-of-range external_id
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);

        let tmp = write_tmp(&buf);
        let err = read_taxo_k2d(tmp.path()).unwrap_err();
        assert!(
            format!("{err:#}").contains("external_id"),
            "expected an external_id range error, got: {err:#}"
        );
    }

    #[test]
    fn test_taxo_k2d_child_range_overflow_errors() {
        // Node 1 declares first_child=0, child_count=999; exceeds node_count.
        // Validation must catch it BEFORE the DFS so we never index out of bounds.
        let name_data: &[u8] = b"root\0";
        let rank_data: &[u8] = b"no rank\0";
        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0),   // placeholder
            (0, 0, 999, 0, 0, 1, 0), // root with bogus child_count
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);

        let tmp = write_tmp(&buf);
        let err = read_taxo_k2d(tmp.path()).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("exceeds node_count"),
            "expected child-range error, got: {msg}"
        );
    }

    #[test]
    fn test_taxo_k2d_cycle_detected() {
        // Construct a node that points back to itself as its child range.
        // The visited-set guard must bail rather than loop forever.
        let name_data: &[u8] = b"root\0a\0";
        let rank_data: &[u8] = b"no rank\0";
        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0),    // placeholder
            (0, 2, 1, 0, 0, 1, 0),    // root -> child node 2
            (1, 1, 1, 5, 0, 9989, 0), // a -> root (cycle)
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);

        let tmp = write_tmp(&buf);
        let err = read_taxo_k2d(tmp.path()).unwrap_err();
        assert!(format!("{err:#}").contains("cycle"));
    }

    #[test]
    fn test_cstr_at_out_of_bounds_errors() {
        let data = b"hello\0";
        let result = cstr_at(data, 100);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("out of bounds"));
    }

    #[test]
    fn test_cstr_at_invalid_utf8_errors() {
        // 0xFF is not valid UTF-8.
        let data = &[0xFF, 0x00];
        let result = cstr_at(data, 0);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("invalid UTF-8"));
    }

    #[test]
    fn test_skips_nodes_with_external_id_zero() {
        // Tree: root(1) -> invisible(ext=0) -> Rodentia(9989) -> Rattus(10116)
        // The invisible node should be skipped; Rodentia promoted to depth 1.
        let name_data: &[u8] = b"root\0invisible\0Rodentia\0Rattus norvegicus\0";
        let rank_data: &[u8] = b"no rank\0order\0species\0";

        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0),       // Node 0: placeholder
            (0, 2, 1, 0, 0, 1, 0),       // Node 1: root
            (1, 3, 1, 5, 0, 0, 0),       // Node 2: invisible (ext_id=0)
            (2, 4, 1, 15, 8, 9989, 0),   // Node 3: Rodentia
            (3, 0, 0, 23, 14, 10116, 0), // Node 4: Rattus
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);

        let tmp = write_tmp(&buf);
        let entries = read_taxo_k2d(tmp.path()).unwrap();

        assert_eq!(entries.len(), 3); // root, Rodentia, Rattus (invisible skipped)
        assert_eq!(entries[1].taxon_id, 9989);
        assert_eq!(entries[1].indent, 2); // promoted from depth 2 -> depth 1
        assert_eq!(entries[2].taxon_id, 10116);
        assert_eq!(entries[2].indent, 4);
    }

    #[test]
    fn test_taxonomy_dmp_produces_entries_in_dfs_order() {
        let tmp_dir = tempfile::TempDir::new().unwrap();
        crate::write_minimal_taxonomy_dmp(tmp_dir.path());
        let entries = read_taxonomy_dmp(tmp_dir.path()).unwrap();

        assert_eq!(entries.len(), 3, "expected root + Rodentia + Rattus");

        assert_eq!(entries[0].taxon_id, 1);
        assert_eq!(entries[0].name, "root");
        assert_eq!(entries[0].indent, 0);
        assert_eq!(entries[0].rank_code, "--"); // "no rank" -> "--"

        assert_eq!(entries[1].taxon_id, 9989);
        assert_eq!(entries[1].name, "Rodentia");
        assert_eq!(entries[1].indent, 2);
        assert_eq!(entries[1].rank_code, "O"); // "order" -> "O"

        assert_eq!(entries[2].taxon_id, 10116);
        assert_eq!(entries[2].name, "Rattus norvegicus");
        assert_eq!(entries[2].indent, 4);
        assert_eq!(entries[2].rank_code, "S"); // "species" -> "S"
    }

    #[test]
    fn test_taxonomy_dmp_builds_valid_taxonomy_tree() {
        use crate::kraken_report::KrakenTaxonomyTree;

        let tmp_dir = tempfile::TempDir::new().unwrap();
        crate::write_minimal_taxonomy_dmp(tmp_dir.path());
        let entries = read_taxonomy_dmp(tmp_dir.path()).unwrap();
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();

        assert!(tree.contains(1));
        assert!(tree.contains(9989));
        assert!(tree.contains(10116));
        assert_eq!(tree.parent_of(9989), Some(1));
        assert_eq!(tree.parent_of(10116), Some(9989));
        assert!(tree.descendants_of(9989).contains(&10116));
    }

    #[test]
    fn test_taxonomy_dmp_ignores_non_scientific_names() {
        let tmp_dir = tempfile::TempDir::new().unwrap();
        let tax_dir = tmp_dir.path().join("taxonomy");
        std::fs::create_dir_all(&tax_dir).unwrap();

        std::fs::write(
            tax_dir.join("nodes.dmp"),
            "1\t|\t1\t|\tno rank\t|\n9989\t|\t1\t|\torder\t|\n",
        )
        .unwrap();
        // Provide a synonym before the scientific name; only scientific name kept.
        std::fs::write(
            tax_dir.join("names.dmp"),
            concat!(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n",
                "9989\t|\tRodents\t|\t\t|\tcommon name\t|\n",
                "9989\t|\tRodentia\t|\t\t|\tscientific name\t|\n",
            ),
        )
        .unwrap();

        let entries = read_taxonomy_dmp(tmp_dir.path()).unwrap();
        let rodentia = entries.iter().find(|e| e.taxon_id == 9989).unwrap();
        assert_eq!(rodentia.name, "Rodentia");
    }

    #[test]
    fn test_taxonomy_dmp_missing_nodes_file_errors() {
        let tmp_dir = tempfile::TempDir::new().unwrap();
        let tax_dir = tmp_dir.path().join("taxonomy");
        std::fs::create_dir_all(&tax_dir).unwrap();
        std::fs::write(
            tax_dir.join("names.dmp"),
            "1\t|\troot\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        assert!(read_taxonomy_dmp(tmp_dir.path()).is_err());
    }

    #[test]
    fn test_taxonomy_dmp_skips_short_or_unparseable_lines() {
        // Lines with fewer than 3 fields, or non-integer taxon IDs in the
        // first two fields, are silently skipped. Verify the rest of the
        // file still parses correctly.
        let tmp_dir = tempfile::TempDir::new().unwrap();
        let tax_dir = tmp_dir.path().join("taxonomy");
        std::fs::create_dir_all(&tax_dir).unwrap();
        std::fs::write(
            tax_dir.join("nodes.dmp"),
            // First line has only 2 fields (skipped).
            // Second line has "abc" as taxon_id (skipped).
            // Third line is valid root.
            "garbage\t|\nabc\t|\t1\t|\torder\t|\n1\t|\t1\t|\tno rank\t|\n",
        )
        .unwrap();
        std::fs::write(
            tax_dir.join("names.dmp"),
            "1\t|\troot\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        let entries = read_taxonomy_dmp(tmp_dir.path()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].taxon_id, 1);
    }

    #[test]
    fn test_taxonomy_dmp_missing_root_errors() {
        // nodes.dmp with no entry for taxon 1 must error out clearly rather
        // than silently producing an empty result.
        let tmp_dir = tempfile::TempDir::new().unwrap();
        let tax_dir = tmp_dir.path().join("taxonomy");
        std::fs::create_dir_all(&tax_dir).unwrap();
        std::fs::write(
            tax_dir.join("nodes.dmp"),
            // Only taxon 2 with parent 2 (self-edge, no taxon 1 anywhere).
            "2\t|\t2\t|\torder\t|\n",
        )
        .unwrap();
        std::fs::write(tax_dir.join("names.dmp"), "").unwrap();
        let err = read_taxonomy_dmp(tmp_dir.path()).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("no entry for root"), "got: {msg}");
    }

    #[test]
    fn test_taxo_k2d_skips_idx_zero_and_out_of_bounds() {
        // Hand-craft a taxo.k2d where node 1 declares its first_child
        // pointing at an idx >= node_count via a separate stack push path
        // is not possible (validation catches it), but we can exercise the
        // `idx == 0` skip arm by having node 1 declare itself as a child.
        //
        // Layout: 2 nodes (placeholder + root), root self-references.
        //   node 0: all zeros (placeholder)
        //   node 1: root, first_child=0, child_count=1 → child idx 0 (placeholder)
        // The DFS pushes (0, 1), then loops; idx==0 branch continues, draining
        // the stack. Only the root entry is emitted.
        let name_data: &[u8] = b"root\0";
        let rank_data: &[u8] = b"no rank\0";
        let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
            (0, 0, 0, 0, 0, 0, 0), // placeholder (idx 0)
            (0, 0, 1, 0, 0, 1, 0), // root, first_child=0 (placeholder), count=1
        ];

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
        buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
        for &(p, fc, cc, no, ro, ei, gp) in node_specs {
            for v in [p, fc, cc, no, ro, ei, gp] {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf.extend_from_slice(name_data);
        buf.extend_from_slice(rank_data);

        let tmp = write_tmp(&buf);
        let entries = read_taxo_k2d(tmp.path()).unwrap();
        assert_eq!(entries.len(), 1, "only root should emit; idx 0 is skipped");
        assert_eq!(entries[0].taxon_id, 1);
    }
}
