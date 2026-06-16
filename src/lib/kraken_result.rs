//! Kraken classification output parsing (v1 and v2).

use std::io::BufRead;

use ahash::AHashMap;
use log::warn;
use std::path::Path;

use anyhow::{bail, Context, Result};

/// A single result from a Kraken classification output file.
///
/// Compatible with both Kraken v1 and Kraken v2 output formats. The parser
/// strips trailing `/1`/`/2` suffixes from query names; a Kraken v1 quirk
/// for paired-end reads. This stripping is harmless on Kraken v2 output,
/// where these suffixes are already absent from the query column.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KrakenResult {
    /// `'C'` for classified, `'U'` for unclassified.
    pub classification_code: char,
    /// Read/query name, stripped of Kraken v1 `/1` or `/2` suffixes.
    pub query_name: String,
    /// NCBI taxonomic ID (`0` if unclassified).
    pub taxon_id: u32,
}

impl KrakenResult {
    /// Returns `true` if the read was classified (taxon ID is non-zero).
    pub fn is_classified(&self) -> bool {
        self.taxon_id != 0
    }

    /// Parse a single tab-delimited line from a Kraken output file.
    ///
    /// Expects exactly 5 tab-delimited fields. Fields 4 and 5 (query length
    /// and k-mer counts) are validated for presence but not parsed.
    pub fn from_line(line: &str) -> Result<Self> {
        let mut it = line.splitn(5, '\t');
        let (Some(code_field), Some(name_field), Some(taxon_field), Some(_), Some(_)) =
            (it.next(), it.next(), it.next(), it.next(), it.next())
        else {
            bail!("expected 5 tab-delimited fields: {:?}", line);
        };

        let classification_code = match code_field {
            "C" => 'C',
            "U" => 'U',
            code => bail!("unknown classification code: {code:?}"),
        };

        let query_name = crate::strip_pair_suffix(name_field).to_owned();

        let taxon_id: u32 = taxon_field
            .parse()
            .with_context(|| format!("failed to parse taxon ID: {taxon_field:?}"))?;

        match (classification_code, taxon_id) {
            ('U', id) if id != 0 => {
                bail!("unclassified read has non-zero taxon ID: {id}")
            }
            ('C', 0) => bail!("classified read has taxon ID 0"),
            _ => {}
        }

        Ok(KrakenResult {
            classification_code,
            query_name,
            taxon_id,
        })
    }

    /// Read all results from a Kraken output file, returning them as an `AHashMap<query_name, taxon_id>`.
    ///
    /// This loads the full file into memory, which is not appropriate for typical sequencing run sizes.
    pub fn load_as_map(path: &Path) -> Result<AHashMap<String, u32>> {
        use std::collections::hash_map::Entry;

        let file = std::fs::File::open(path).with_context(|| {
            format!("failed to open Kraken assignments file: {}", path.display())
        })?;
        let mut reader = std::io::BufReader::new(file);
        let mut map = AHashMap::new();
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
            let result = KrakenResult::from_line(line)
                .with_context(|| format!("failed to parse Kraken result at line {line_no}"))?;
            // After stripping /1//2, the same base name may appear twice in paired Kraken v1
            // output. That is only OK when both entries assign the same taxon ID.
            match map.entry(result.query_name) {
                Entry::Vacant(slot) => {
                    slot.insert(result.taxon_id);
                }
                Entry::Occupied(slot) => {
                    let existing = *slot.get();
                    if existing != result.taxon_id {
                        anyhow::bail!(
                            "duplicate read name {:?} at line {} with conflicting taxon IDs {} and {}; \
                             ensure Kraken output has one classification per template",
                            slot.key(),
                            line_no,
                            existing,
                            result.taxon_id,
                        );
                    }
                }
            }
        }
        Ok(map)
    }
}

/// Streaming merge-join with a self-sizing lookahead buffer for mostly-ordered
/// Kraken assignment streams. Used by both `annotate` (driving from
/// SAM/BAM/CRAM records) and `filter` (driving from FASTX records).
///
/// Advances the Kraken reader once per new query name, reusing the cached
/// taxon for consecutive lookups sharing the same name (e.g. interleaved
/// /1 + /2 FASTQ pairs, or paired SAM segments). When the Kraken stream runs
/// ahead of the input; as happens with Kraken v1 using multiple threads,
/// which flushes per-thread work-unit buffers in completion order;
/// mismatched entries are placed in a `HashMap` lookahead buffer. On the
/// next lookup the buffer is checked first; a hit removes the entry, so the
/// buffer grows only as deep as the actual disorder and shrinks back as
/// matching input records arrive.
///
/// On EOF, missing names return `Ok(None)`. Callers decide whether absence
/// is fatal (annotate) or treated as taxon 0 (filter).
///
/// Disagreeing duplicate entries in the buffer (same query name with
/// different taxon IDs) are reported as an error, matching the behaviour of
/// `KrakenResult::load_as_map`. Same-base-name pairs after `/1 //2` stripping
/// (Kraken v1 quirk) are fine as long as their taxons agree.
pub(crate) struct StreamingLookup<R: BufRead> {
    reader: R,
    line_buf: String,
    line_no: usize,
    current_name: String,
    current_taxon_id: Option<u32>,
    /// `true` once `next_kraken` has reached EOF; subsequent lookups consult
    /// only the buffer.
    exhausted: bool,
    lookahead: AHashMap<String, u32>,
    warned: bool,
}

/// Emit a single warning when the lookahead buffer reaches this size.
/// Kraken v1 with N threads produces disorder bounded by ~(N-1) × work-unit
/// size; for N=32 with short reads this stays well below 100,000. A buffer
/// this large almost always means the input files are mismatched.
const WARN_LOOKAHEAD: usize = 1_000_000;

impl<R: BufRead> StreamingLookup<R> {
    pub(crate) fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: String::new(),
            line_no: 0,
            current_name: String::new(),
            current_taxon_id: None,
            exhausted: false,
            lookahead: AHashMap::new(),
            warned: false,
        }
    }

    /// Read the next non-blank Kraken result line. Returns `Ok(None)` at EOF
    /// and sets the `exhausted` flag.
    fn next_kraken(&mut self) -> Result<Option<KrakenResult>> {
        loop {
            self.line_buf.clear();
            let n = self
                .reader
                .read_line(&mut self.line_buf)
                .context("failed to read Kraken assignments line")?;
            if n == 0 {
                self.exhausted = true;
                return Ok(None);
            }
            self.line_no += 1;
            let trimmed = self.line_buf.trim_end_matches(['\n', '\r']);
            if trimmed.trim().is_empty() {
                continue;
            }
            return KrakenResult::from_line(trimmed)
                .with_context(|| format!("failed to parse Kraken result at line {}", self.line_no))
                .map(Some);
        }
    }

    /// Resolve the taxon id for `name`, advancing/buffering the Kraken stream
    /// as needed. Returns `Ok(Some(_))` on hit, `Ok(None)` if `name` is absent
    /// from the entire assignments stream.
    pub(crate) fn lookup(&mut self, name: &str) -> Result<Option<u32>> {
        if name == self.current_name {
            return Ok(self.current_taxon_id);
        }

        // Reset cache so a miss for this name doesn't leak forward.
        self.current_name.clear();
        self.current_name.push_str(name);
        self.current_taxon_id = None;

        // Buffer hit: entry was buffered while resolving an earlier record.
        if let Some(tid) = self.lookahead.remove(name) {
            self.current_taxon_id = Some(tid);
            return Ok(Some(tid));
        }

        // If the stream has already been fully consumed, the name is missing.
        if self.exhausted {
            return Ok(None);
        }

        // Advance, buffering out-of-order entries, until found or EOF.
        loop {
            let Some(kr) = self.next_kraken()? else {
                return Ok(None);
            };
            if kr.query_name == name {
                self.current_taxon_id = Some(kr.taxon_id);
                return Ok(Some(kr.taxon_id));
            }
            // Same-base-name duplicates with disagreeing taxons are an error
            // (matches the unordered-load behaviour of `load_as_map`); equal
            // taxons collapse silently; common with Kraken v1 paired output
            // where /1 and /2 share a base name and the same classification.
            if let Some(prev) = self.lookahead.get(&kr.query_name) {
                if *prev != kr.taxon_id {
                    bail!(
                        "duplicate query name {:?} with disagreeing taxon ({prev} vs {})",
                        kr.query_name,
                        kr.taxon_id,
                    );
                }
            } else {
                self.lookahead.insert(kr.query_name, kr.taxon_id);
            }

            if !self.warned && self.lookahead.len() >= WARN_LOOKAHEAD {
                warn!(
                    "Kraken assignments lookahead buffer reached {WARN_LOOKAHEAD} entries while \
                     resolving read {name:?}; the assignments file may be significantly \
                     disordered relative to the input, or the files may be mismatched."
                );
                self.warned = true;
            }
        }
    }

    /// Count assignments that were left in the buffer or never read off the
    /// stream after the input was fully consumed. Drains the remaining stream.
    pub(crate) fn count_unconsumed(&mut self) -> Result<usize> {
        let mut n = self.lookahead.len();
        self.lookahead.clear();
        while !self.exhausted {
            // `next_kraken` sets `exhausted = true` and returns `Ok(None)` on
            // EOF; a malformed line propagates rather than being swallowed
            // (which would also undercount the leftovers).
            if self.next_kraken()?.is_some() {
                n += 1;
            }
        }
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_classified_v1() {
        let line = "C\tread1\t9606\t150\t9606:10 1:5 9606:20";
        let result = KrakenResult::from_line(line).unwrap();
        assert_eq!(result.classification_code, 'C');
        assert_eq!(result.query_name, "read1");
        assert_eq!(result.taxon_id, 9606);
        assert!(result.is_classified());
    }

    #[test]
    fn test_parse_unclassified() {
        let line = "U\tread2\t0\t100\t0:70";
        let result = KrakenResult::from_line(line).unwrap();
        assert_eq!(result.classification_code, 'U');
        assert_eq!(result.taxon_id, 0);
        assert!(!result.is_classified());
    }

    #[test]
    fn test_parse_v1_read_pair_suffix_stripped() {
        let line = "C\tread3/1\t9606\t150\t9606:75";
        let result = KrakenResult::from_line(line).unwrap();
        assert_eq!(result.query_name, "read3");
    }

    #[test]
    fn test_parse_v2_paired_format_accepted() {
        // Kraken v2 paired format (length as `150|150`, kmer column with `|:|` separator).
        let line = "C\tread4\t9606\t150|150\t9606:50 |:| 9606:50";
        let result = KrakenResult::from_line(line).unwrap();
        assert_eq!(result.taxon_id, 9606);
    }

    #[test]
    fn test_invalid_classification_code() {
        let line = "X\tread6\t9606\t100\t9606:10";
        assert!(KrakenResult::from_line(line).is_err());
    }

    #[test]
    fn test_unclassified_nonzero_taxon_id_is_error() {
        let line = "U\tread7\t9606\t100\t9606:10";
        assert!(KrakenResult::from_line(line).is_err());
    }

    #[test]
    fn test_classified_zero_taxon_id_is_error() {
        let line = "C\tread8\t0\t100\t0:10";
        assert!(KrakenResult::from_line(line).is_err());
    }

    #[test]
    fn test_too_few_fields_is_error() {
        let line = "C\tread9\t9606\t100";
        assert!(KrakenResult::from_line(line).is_err());
    }

    #[test]
    fn test_load_as_map() {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "C\tread1\t9606\t150\t9606:75").unwrap();
        writeln!(tmp, "U\tread2\t0\t100\t0:50").unwrap();
        let map = KrakenResult::load_as_map(tmp.path()).unwrap();
        assert_eq!(map["read1"], 9606);
        assert_eq!(map["read2"], 0);
    }

    #[test]
    fn test_load_as_map_skips_blank_lines() {
        // Blank and whitespace-only lines must be skipped without errors.
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp).unwrap();
        writeln!(tmp, "   ").unwrap();
        writeln!(tmp, "C\tread1\t9606\t150\t9606:75").unwrap();
        writeln!(tmp).unwrap();
        let map = KrakenResult::load_as_map(tmp.path()).unwrap();
        assert_eq!(map.len(), 1);
        assert_eq!(map["read1"], 9606);
    }

    #[test]
    fn test_load_as_map_missing_file_errors() {
        let err =
            KrakenResult::load_as_map(std::path::Path::new("/nonexistent/file.tsv")).unwrap_err();
        assert!(format!("{err:#}").contains("failed to open"));
    }

    #[test]
    fn test_duplicate_name_same_taxon_is_ok() {
        // Kraken v1 paired output: read1/1 and read1/2 both map to 9606.
        // After suffix stripping they share the same base name; this is fine.
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "C\tread1/1\t9606\t150\t9606:75").unwrap();
        writeln!(tmp, "C\tread1/2\t9606\t150\t9606:75").unwrap();
        let map = KrakenResult::load_as_map(tmp.path()).unwrap();
        assert_eq!(map["read1"], 9606);
    }

    #[test]
    fn test_duplicate_name_different_taxon_is_error() {
        // Same base name (after suffix stripping) but different taxon IDs: must error.
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "C\tread1/1\t9606\t150\t9606:75").unwrap();
        writeln!(tmp, "C\tread1/2\t1234\t150\t1234:75").unwrap();
        let result = KrakenResult::load_as_map(tmp.path());
        assert!(result.is_err());
        let msg = format!("{:#}", result.unwrap_err());
        assert!(msg.contains("duplicate read name"));
        assert!(msg.contains("conflicting taxon IDs"));
    }

    fn streaming(input: &str) -> StreamingLookup<std::io::Cursor<&[u8]>> {
        StreamingLookup::new(std::io::Cursor::new(input.as_bytes()))
    }

    #[test]
    fn test_streaming_lookup_in_order() {
        let kraken = "C\tr1\t9606\t100\t9606:5\nC\tr2\t1234\t100\t1234:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(9606));
        assert_eq!(s.lookup("r2").unwrap(), Some(1234));
    }

    #[test]
    fn test_streaming_lookup_caches_repeated_name() {
        // /1 + /2 paired records lookup the same base name twice; the second
        // lookup must hit the cache and not advance the stream.
        let kraken = "C\tr1\t9606\t100\t9606:5\nC\tr2\t1234\t100\t1234:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(9606));
        assert_eq!(s.lookup("r1").unwrap(), Some(9606));
        assert_eq!(s.lookup("r2").unwrap(), Some(1234));
    }

    #[test]
    fn test_streaming_lookup_buffers_out_of_order() {
        // Input drives r1, r2, r3, but stream is r2, r1, r3 (disordered).
        let kraken = "C\tr2\t2\t100\t2:5\nC\tr1\t1\t100\t1:5\nC\tr3\t3\t100\t3:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(s.lookup("r2").unwrap(), Some(2));
        assert_eq!(s.lookup("r3").unwrap(), Some(3));
    }

    #[test]
    fn test_streaming_lookup_missing_returns_none() {
        let kraken = "C\tr1\t1\t100\t1:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(s.lookup("missing").unwrap(), None);
    }

    #[test]
    fn test_streaming_lookup_caches_misses_after_exhaustion() {
        let kraken = "C\tr1\t1\t100\t1:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(s.lookup("missing").unwrap(), None);
        // Repeated miss for same name still returns None without re-walking.
        assert_eq!(s.lookup("missing").unwrap(), None);
    }

    #[test]
    fn test_streaming_lookup_skips_blank_lines() {
        let kraken = "\n  \nC\tr1\t1\t100\t1:5\n\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
    }

    #[test]
    fn test_streaming_lookup_disagreeing_buffer_duplicate_errors() {
        // Two buffered entries with the same name but disagreeing taxons;
        // matches the load_as_map invariant.
        let kraken = "C\tr1\t1\t100\t1:5\nC\tr1\t2\t100\t2:5\nC\tr2\t99\t100\t99:5\n";
        let mut s = streaming(kraken);
        let err = s.lookup("r2").unwrap_err();
        assert!(format!("{err:#}").contains("disagreeing taxon"));
    }

    #[test]
    fn test_streaming_lookup_same_taxon_buffer_duplicate_collapses() {
        // Same name twice with the SAME taxon; silently collapses, no error.
        let kraken =
            "C\tdup\t9606\t100\t9606:5\nC\tdup\t9606\t100\t9606:5\nC\ttarget\t5\t100\t5:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("target").unwrap(), Some(5));
        assert_eq!(s.lookup("dup").unwrap(), Some(9606));
    }

    #[test]
    fn test_streaming_lookup_strips_v1_pair_suffix() {
        // KrakenResult::from_line strips trailing /1 or /2; streaming lookup
        // inherits that and matches the base name from the FASTQ side.
        let kraken = "C\tread1/1\t9606\t100\t9606:5\nC\tread2/2\t1234\t100\t1234:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("read1").unwrap(), Some(9606));
        assert_eq!(s.lookup("read2").unwrap(), Some(1234));
    }

    #[test]
    fn test_streaming_lookup_buffer_resolved_after_exhaustion() {
        // After exhaustion, a buffered entry from earlier walks is still
        // resolvable. Models the case where input drives "missing" first
        // (drains the stream, buffering r1), then asks for r1.
        let kraken = "C\tr1\t9606\t100\t9606:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("missing").unwrap(), None);
        assert_eq!(s.lookup("r1").unwrap(), Some(9606));
        assert_eq!(s.lookup("also_missing").unwrap(), None);
    }

    #[test]
    fn test_streaming_lookup_count_unconsumed_drains_tail() {
        let kraken = "C\tr1\t1\t100\t1:5\nC\tr2\t2\t100\t2:5\nC\tr3\t3\t100\t3:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(s.count_unconsumed().unwrap(), 2);
    }

    #[test]
    fn test_streaming_lookup_count_unconsumed_includes_buffered() {
        let kraken = "C\tr2\t2\t100\t2:5\nC\tr1\t1\t100\t1:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(
            s.count_unconsumed().unwrap(),
            1,
            "r2 should still be in buffer"
        );
    }

    #[test]
    fn test_streaming_lookup_count_unconsumed_propagates_parse_error() {
        // A malformed line in the unconsumed tail must surface as an error, not
        // be silently swallowed (which would also undercount the leftovers).
        let kraken = "C\tr1\t1\t100\t1:5\nTHIS_LINE_IS_NOT_VALID\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert!(
            s.count_unconsumed().is_err(),
            "a malformed tail line must propagate as an error"
        );
    }

    #[test]
    fn test_streaming_lookup_count_unconsumed_zero_when_drained() {
        let kraken = "C\tr1\t1\t100\t1:5\n";
        let mut s = streaming(kraken);
        assert_eq!(s.lookup("r1").unwrap(), Some(1));
        assert_eq!(s.lookup("missing").unwrap(), None);
        assert_eq!(s.count_unconsumed().unwrap(), 0);
    }
}
