//! Kraken report parsing and in-memory taxonomy tree.

use std::io::BufRead;

use ahash::{AHashMap, AHashSet};
use std::path::Path;

use anyhow::{bail, Context, Result};

/// Number of leading spaces per depth level in the indent-prefixed name column
/// of a Kraken report (or `taxo.k2d` rendering).
pub const SPACES_PER_DEPTH: usize = 2;

/// A single entry from a `kraken-report` output file.
///
/// This represents one row of `kraken2 --report` output: six tab-delimited
/// columns (percent, clade fragment count, direct fragment count, rank code,
/// taxon ID, and indent-prefixed scientific name). The leading whitespace on
/// the name column encodes tree depth; exactly [`SPACES_PER_DEPTH`] spaces
/// per level; and is preserved here in [`indent`](Self::indent) for
/// [`KrakenTaxonomyTree`] reconstruction.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct KrakenReportEntry {
    /// Percentage of fragments in the clade rooted at this taxon.
    pub pct_fragments: f64,
    /// Number of fragments in the clade rooted at this taxon.
    pub num_fragments_clade: u64,
    /// Number of fragments directly assigned to this taxon.
    pub num_fragments_direct: u64,
    /// Rank code: U/R/D/K/P/C/O/F/G/S or an intermediate code (e.g. `G2`, `--`).
    pub rank_code: String,
    /// NCBI taxonomic ID.
    pub taxon_id: u32,
    /// Scientific name (trimmed of leading whitespace).
    pub name: String,
    /// Leading space count in the name column; encodes tree depth
    /// ([`SPACES_PER_DEPTH`] spaces per level).
    pub indent: usize,
}

impl KrakenReportEntry {
    /// Parse a single tab-delimited line from a Kraken report file.
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            bail!(
                "expected at least 6 tab-delimited fields, found {}: {:?}",
                fields.len(),
                line
            );
        }
        if fields.len() > 6 {
            log::debug!(
                "ignoring {} trailing tab-delimited field(s) in Kraken report line",
                fields.len() - 6
            );
        }

        let pct_fragments = fields[0]
            .trim()
            .parse::<f64>()
            .with_context(|| format!("failed to parse pct_fragments: {:?}", fields[0]))?;
        let num_fragments_clade = fields[1]
            .trim()
            .parse::<u64>()
            .with_context(|| format!("failed to parse num_fragments_clade: {:?}", fields[1]))?;
        let num_fragments_direct = fields[2]
            .trim()
            .parse::<u64>()
            .with_context(|| format!("failed to parse num_fragments_direct: {:?}", fields[2]))?;
        let rank_code = fields[3].trim().to_owned();
        let taxon_id = fields[4]
            .trim()
            .parse::<u32>()
            .with_context(|| format!("failed to parse taxon_id: {:?}", fields[4]))?;

        // Indent is the number of leading ASCII spaces (Kraken uses spaces only).
        let name_field = fields[5];
        let indent = name_field.bytes().take_while(|&b| b == b' ').count();
        if indent % SPACES_PER_DEPTH != 0 {
            bail!(
                "indent of {} space(s) is not a multiple of {} (SPACES_PER_DEPTH); \
                 line is likely tab-indented or malformed: {:?}",
                indent,
                SPACES_PER_DEPTH,
                line
            );
        }
        let name = name_field.trim().to_owned();

        Ok(KrakenReportEntry {
            pct_fragments,
            num_fragments_clade,
            num_fragments_direct,
            rank_code,
            taxon_id,
            name,
            indent,
        })
    }

    /// Read all entries from a Kraken report file.
    pub fn read_file(path: &Path) -> Result<Vec<KrakenReportEntry>> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open Kraken report file: {}", path.display()))?;
        let mut reader = std::io::BufReader::new(file);
        let mut entries = Vec::new();
        let mut buf = String::new();
        let mut line_no: usize = 0;
        loop {
            buf.clear();
            line_no += 1;
            let n = reader
                .read_line(&mut buf)
                .with_context(|| format!("failed to read line {line_no}"))?;
            if n == 0 {
                break;
            }
            let line = buf.trim_end_matches(['\n', '\r']);
            if line.trim().is_empty() {
                continue;
            }
            let entry = KrakenReportEntry::from_line(line).with_context(|| {
                format!("failed to parse Kraken report entry at line {line_no}")
            })?;
            entries.push(entry);
        }
        Ok(entries)
    }
}

/// An in-memory NCBI taxonomy tree reconstructed from a Kraken report.
///
/// Tree structure is derived from the indentation levels in the report file.
/// The tree is not thread-safe; clone for parallel use.
#[derive(Debug)]
pub struct KrakenTaxonomyTree {
    /// Maps taxon ID to its parent taxon ID (`None` for root).
    parents: AHashMap<u32, Option<u32>>,
    /// Maps taxon ID to its set of child taxon IDs.
    children: AHashMap<u32, AHashSet<u32>>,
}

impl KrakenTaxonomyTree {
    /// Build a taxonomy tree from a slice of report entries in report-file order.
    ///
    /// Entries with `taxon_id == 0` (Kraken's "unclassified" marker) are
    /// silently filtered out before tree construction.
    pub fn from_entries(entries: &[KrakenReportEntry]) -> Result<Self> {
        let mut parents: AHashMap<u32, Option<u32>> = AHashMap::new();
        let mut children: AHashMap<u32, AHashSet<u32>> = AHashMap::new();
        let mut seen: AHashSet<u32> = AHashSet::new();

        // Skip unclassified (taxon_id == 0).
        let mut iter = entries.iter().filter(|e| e.taxon_id != 0).peekable();

        let root = match iter.next() {
            Some(e) => e,
            None => return Ok(KrakenTaxonomyTree { parents, children }),
        };

        seen.insert(root.taxon_id);
        parents.insert(root.taxon_id, None);
        children.insert(root.taxon_id, AHashSet::new());

        // Map from indent -> parent taxon_id, updated as we walk the entries.
        let mut parent_per_indent: AHashMap<usize, u32> = AHashMap::new();
        let mut last = root;

        for entry in iter {
            if !seen.insert(entry.taxon_id) {
                bail!(
                    "duplicate taxon ID {} (name={:?}) encountered while building taxonomy tree",
                    entry.taxon_id,
                    entry.name,
                );
            }

            if entry.indent > last.indent {
                parent_per_indent.insert(entry.indent, last.taxon_id);
            }

            let &parent_id = parent_per_indent.get(&entry.indent).ok_or_else(|| {
                anyhow::anyhow!(
                    "could not find parent for taxon {} (name={:?}) at indent level {}; \
                     expected an ancestor at a lower indent level to have been seen first. \
                     Indent levels must increase by exactly {} per depth level ({} spaces per \
                     level). Check that this report was produced by a standard Kraken2 run.",
                    entry.taxon_id,
                    entry.name,
                    entry.indent,
                    SPACES_PER_DEPTH,
                    SPACES_PER_DEPTH,
                )
            })?;

            parents.insert(entry.taxon_id, Some(parent_id));
            children
                .entry(parent_id)
                .or_default()
                .insert(entry.taxon_id);

            last = entry;
        }

        let max_depth = parents.len();
        for &start in parents.keys() {
            let mut cur = parents.get(&start).copied().flatten();
            let mut steps = 0usize;
            while let Some(id) = cur {
                steps += 1;
                if steps > max_depth {
                    bail!(
                        "cycle detected in taxonomy tree while walking ancestors of taxon {start}"
                    );
                }
                cur = parents.get(&id).copied().flatten();
            }
        }

        Ok(KrakenTaxonomyTree { parents, children })
    }

    /// Build a taxonomy tree from a Kraken report file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let entries = KrakenReportEntry::read_file(path)?;
        Self::from_entries(&entries)
    }

    /// Returns `true` if the taxon ID is present in the tree.
    pub fn contains(&self, taxon_id: u32) -> bool {
        self.parents.contains_key(&taxon_id)
    }

    /// Returns the parent taxon ID, or `None` if this is the root or unknown.
    pub fn parent_of(&self, taxon_id: u32) -> Option<u32> {
        self.parents.get(&taxon_id).copied().flatten()
    }

    /// Returns the children of the given taxon ID.
    pub fn children_of(&self, taxon_id: u32) -> &AHashSet<u32> {
        static EMPTY: std::sync::OnceLock<AHashSet<u32>> = std::sync::OnceLock::new();
        self.children
            .get(&taxon_id)
            .unwrap_or_else(|| EMPTY.get_or_init(AHashSet::new))
    }

    /// Returns `true` if `taxon_id` is an ancestor of any taxon in `targets`.
    ///
    /// Walks up the ancestry chain of each target without allocating.
    ///
    /// Time complexity: O(|targets| × tree_depth).
    pub fn is_ancestor_of_any(&self, taxon_id: u32, targets: &AHashSet<u32>) -> bool {
        for &target in targets {
            let mut cur = self.parent_of(target);
            while let Some(id) = cur {
                if id == taxon_id {
                    return true;
                }
                cur = self.parent_of(id);
            }
        }
        false
    }

    /// Returns the set of all descendant taxon IDs (not including `taxon_id` itself).
    ///
    /// Uses a depth-first traversal of the children map.
    /// Returns an empty set if `taxon_id` is not in the tree.
    pub fn descendants_of(&self, taxon_id: u32) -> AHashSet<u32> {
        let mut result = AHashSet::new();
        let mut stack = vec![taxon_id];
        while let Some(id) = stack.pop() {
            for &child in self.children_of(id) {
                if result.insert(child) {
                    stack.push(child);
                }
            }
        }
        result
    }

    /// Returns the number of nodes in the tree.
    pub fn len(&self) -> usize {
        self.parents.len()
    }

    /// Returns `true` if the tree contains no nodes.
    pub fn is_empty(&self) -> bool {
        self.parents.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(taxon_id: u32, indent: usize, name: &str) -> KrakenReportEntry {
        KrakenReportEntry {
            pct_fragments: 0.0,
            num_fragments_clade: 0,
            num_fragments_direct: 0,
            rank_code: "S".to_owned(),
            taxon_id,
            name: name.to_owned(),
            indent,
        }
    }

    /// Build a small tree: root(1) -> A(2) -> B(3), root(1) -> C(4)
    fn small_tree() -> KrakenTaxonomyTree {
        let entries = vec![
            make_entry(1, 0, "root"),
            make_entry(2, 2, "A"),
            make_entry(3, 4, "B"),
            make_entry(4, 2, "C"),
        ];
        KrakenTaxonomyTree::from_entries(&entries).unwrap()
    }

    #[test]
    fn test_tree_contains() {
        let tree = small_tree();
        assert!(tree.contains(1));
        assert!(tree.contains(2));
        assert!(tree.contains(3));
        assert!(!tree.contains(99));
    }

    #[test]
    fn test_tree_parent_of() {
        let tree = small_tree();
        assert_eq!(tree.parent_of(1), None);
        assert_eq!(tree.parent_of(2), Some(1));
        assert_eq!(tree.parent_of(3), Some(2));
        assert_eq!(tree.parent_of(4), Some(1));
    }

    #[test]
    fn test_tree_is_ancestor_of_any() {
        let tree = small_tree();
        // root(1) -> A(2) -> B(3), root(1) -> C(4)
        let targets: AHashSet<u32> = [3].into_iter().collect(); // B
        assert!(tree.is_ancestor_of_any(2, &targets)); // A is ancestor of B
        assert!(tree.is_ancestor_of_any(1, &targets)); // root is ancestor of B
        assert!(!tree.is_ancestor_of_any(4, &targets)); // C is NOT ancestor of B
        assert!(!tree.is_ancestor_of_any(3, &targets)); // B is not its own ancestor

        let two_targets: AHashSet<u32> = [3, 4].into_iter().collect(); // B and C
        assert!(tree.is_ancestor_of_any(1, &two_targets)); // root is ancestor of both
        assert!(tree.is_ancestor_of_any(2, &two_targets)); // A is ancestor of B
        assert!(!tree.is_ancestor_of_any(4, &two_targets)); // C is not ancestor of B or C
    }

    #[test]
    fn test_tree_children_of() {
        let tree = small_tree();
        assert!(tree.children_of(1).contains(&2));
        assert!(tree.children_of(1).contains(&4));
        assert!(tree.children_of(2).contains(&3));
        assert!(tree.children_of(3).is_empty());
    }

    #[test]
    fn test_read_file_skips_blank_lines() {
        // Blank lines and pure-whitespace lines must be skipped; only real
        // entries appear in the result.
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("report.k2");
        std::fs::write(
            &p,
            // blank, whitespace-only, valid root, blank, valid leaf, blank
            "\n   \n100.00\t1\t1\tR\t1\troot\n\n  0.00\t1\t1\tS\t9606\t  Homo sapiens\n\n",
        )
        .unwrap();
        let entries = KrakenReportEntry::read_file(&p).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].taxon_id, 1);
        assert_eq!(entries[1].taxon_id, 9606);
    }

    #[test]
    fn test_from_line_accepts_trailing_extra_fields() {
        // A line with more than 6 fields (e.g. --report-minimizer-data has 8)
        // is parsed by taking the first 6 fields and emitting a debug log.
        let line = "100.00\t1\t1\tR\t1\troot\textra1\textra2";
        let e = KrakenReportEntry::from_line(line).unwrap();
        assert_eq!(e.taxon_id, 1);
        assert_eq!(e.name, "root");
    }

    #[test]
    fn test_tree_skips_unclassified() {
        let mut entries = vec![make_entry(0, 0, "unclassified")];
        entries.extend(vec![make_entry(1, 0, "root"), make_entry(2, 2, "A")]);
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        assert!(!tree.contains(0));
        assert!(tree.contains(1));
    }

    #[test]
    fn test_tree_len_and_is_empty() {
        let tree = small_tree();
        assert_eq!(tree.len(), 4);
        assert!(!tree.is_empty());

        let empty = KrakenTaxonomyTree::from_entries(&[]).unwrap();
        assert!(empty.is_empty());
    }

    #[test]
    fn test_report_entry_from_line() {
        let line = " 50.00\t1000\t100\tS\t9606\t  Homo sapiens";
        let entry = KrakenReportEntry::from_line(line).unwrap();
        assert_eq!(entry.taxon_id, 9606);
        assert_eq!(entry.name, "Homo sapiens");
        assert_eq!(entry.indent, 2);
    }

    #[test]
    fn test_report_entry_from_line_bad_fields() {
        let line = "50.00\t1000\t100\tS\t9606";
        assert!(KrakenReportEntry::from_line(line).is_err());
    }

    #[test]
    fn test_descendants_of_root() {
        // small_tree(): root(1) -> A(2) -> B(3), root(1) -> C(4)
        let tree = small_tree();
        let desc = tree.descendants_of(1);
        assert!(desc.contains(&2), "A should be a descendant of root");
        assert!(desc.contains(&3), "B should be a descendant of root");
        assert!(desc.contains(&4), "C should be a descendant of root");
        assert!(!desc.contains(&1), "root should not be its own descendant");
        assert_eq!(desc.len(), 3);
    }

    #[test]
    fn test_descendants_of_leaf_is_empty() {
        let tree = small_tree();
        assert!(tree.descendants_of(3).is_empty());
    }

    #[test]
    fn test_descendants_of_unknown_taxon_is_empty() {
        let tree = small_tree();
        assert!(tree.descendants_of(99).is_empty());
    }

    #[test]
    fn test_descendants_of_interior_node() {
        // small_tree(): root(1) -> A(2) -> B(3), root(1) -> C(4)
        // descendants_of(2) should return {3} only; not bleed into sibling subtree {4}.
        let tree = small_tree();
        let desc = tree.descendants_of(2);
        assert!(desc.contains(&3), "B should be a descendant of A");
        assert!(
            !desc.contains(&4),
            "C is in a sibling subtree, not a descendant of A"
        );
        assert_eq!(desc.len(), 1);
    }

    #[test]
    fn test_read_file() {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "100.00\t2000\t0\tR\t1\t  root").unwrap();
        writeln!(tmp, " 50.00\t1000\t100\tS\t9606\t    Homo sapiens").unwrap();
        let entries = KrakenReportEntry::read_file(tmp.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].taxon_id, 1);
        assert_eq!(entries[1].taxon_id, 9606);
    }

    #[test]
    fn test_from_line_tab_in_indent_errors() {
        // A leading tab inside the name field is not an ASCII space, so the
        // indent count is 0; but the field also embeds another tab, which
        // produces 7 fields (still parseable). The real B15 case: a line
        // whose name field starts with an odd number of leading spaces is
        // not a valid Kraken indent and must error.
        let line = "100.00\t2000\t0\tR\t1\t Bacteria";
        let err = KrakenReportEntry::from_line(line).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("not a multiple"),
            "expected odd-indent error, got: {msg}"
        );
    }

    #[test]
    fn test_from_line_accepts_extra_columns() {
        // Eight tab-delimited fields: the original 6 plus two extras. C9
        // says we should accept these and ignore trailing fields.
        let line = "100.00\t2000\t0\tR\t1\t  root\textra1\textra2";
        let entry = KrakenReportEntry::from_line(line).unwrap();
        assert_eq!(entry.taxon_id, 1);
        assert_eq!(entry.name, "root");
        assert_eq!(entry.indent, 2);
    }

    #[test]
    fn test_tree_duplicate_taxon_id_errors() {
        let entries = vec![
            make_entry(1, 0, "root"),
            make_entry(2, 2, "A"),
            make_entry(2, 2, "A-again"),
        ];
        let err = KrakenTaxonomyTree::from_entries(&entries).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("duplicate taxon ID"), "got: {msg}");
    }

    #[test]
    fn test_tree_orphan_indent_errors() {
        // Third entry sits at indent 2, but no entry at indent 2 has been
        // seen as a parent yet; root is at 0, second went to 4. The 2-level
        // is an orphan.
        let entries = vec![
            make_entry(1, 0, "root"),
            make_entry(2, 4, "deep"),
            make_entry(3, 2, "orphan"),
        ];
        let err = KrakenTaxonomyTree::from_entries(&entries).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("could not find parent"), "got: {msg}");
    }

    #[test]
    fn test_entry_json_round_trip() {
        let entry = KrakenReportEntry {
            pct_fragments: 50.0,
            num_fragments_clade: 1000,
            num_fragments_direct: 100,
            rank_code: "S".to_owned(),
            taxon_id: 9606,
            name: "Homo sapiens".to_owned(),
            indent: 4,
        };
        let json = serde_json::to_string(&entry).unwrap();
        let recovered: KrakenReportEntry = serde_json::from_str(&json).unwrap();
        assert_eq!(recovered.taxon_id, 9606);
        assert_eq!(recovered.name, "Homo sapiens");
        assert_eq!(recovered.indent, 4);
    }
}
