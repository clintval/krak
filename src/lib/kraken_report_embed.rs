//! Embed and extract the Kraken taxonomy report in SAM/BAM/CRAM `@CO` header lines.
//!
//! # Encoding choice
//!
//! The Kraken report is serialized as JSON and then base64-encoded into one or
//! more `@CO` header comments. JSON-then-base64 is chosen for human
//! debuggability: a developer can read a header with `samtools view -H`, copy
//! the base64 payload, decode it once, and inspect the resulting JSON with any
//! standard tool. A binary format such as bincode would be roughly 5x smaller
//! but would render the embedded data opaque to inspection. The chunk encoding
//! described below keeps the size cost of the verbose representation bounded
//! when reports are large.
//!
//! # Chunk format
//!
//! Each comment uses the form `krak:report:<i>/<n>:<b64>`, where `<i>` is a
//! 1-based chunk index and `<n>` is the total number of chunks. A reader
//! collects every comment matching the prefix, groups by `<n>`, sorts by
//! `<i>`, concatenates the base64 payloads, and decodes once. The legacy
//! single-comment form `krak:report:<b64>` (no `<i>/<n>:` header) is also
//! accepted on read for backward compatibility and is treated as a single
//! chunk.

use anyhow::{Context, Result};
use base64::{engine::general_purpose::STANDARD, Engine as _};
use noodles::sam;

use crate::kraken_report::KrakenReportEntry;

/// The `@CO` comment prefix that identifies an embedded Kraken report.
pub(crate) const HEADER_COMMENT_PREFIX: &str = "krak:report:";

/// Default chunk size for base64 payload (characters per chunk).
#[cfg(not(test))]
const CHUNK_SIZE: usize = 60 * 1024;

/// Test-only override of the chunk threshold so the multi-chunk path can be
/// exercised without producing huge fixtures.
#[cfg(test)]
const CHUNK_SIZE: usize = 32;

/// Serialize `entries` to one or more strings suitable for SAM `@CO` header lines.
///
/// Pipeline: `Vec<KrakenReportEntry>` -> JSON -> base64 -> chunked into
/// `"krak:report:<i>/<n>:<b64>"` comments.
pub(crate) fn entries_to_header_comment(entries: &[KrakenReportEntry]) -> Result<Vec<String>> {
    let json = serde_json::to_string(entries)
        .context("failed to serialize Kraken report entries to JSON")?;
    let b64 = STANDARD.encode(json.as_bytes());

    let bytes = b64.as_bytes();
    let total = bytes.len().div_ceil(CHUNK_SIZE).max(1);
    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        let start = i * CHUNK_SIZE;
        let end = (start + CHUNK_SIZE).min(bytes.len());
        let slice = &b64[start..end];
        out.push(format!(
            "{HEADER_COMMENT_PREFIX}{}/{}:{}",
            i + 1,
            total,
            slice
        ));
    }
    Ok(out)
}

/// Parsed components of one matching `@CO` comment.
#[derive(Debug)]
struct ParsedChunk {
    index: usize,
    total: usize,
    b64: String,
}

/// Parse a `krak:report:` comment body (the part after the prefix). Accepts
/// both the chunked form `<i>/<n>:<b64>` and the legacy `<b64>` form (treated
/// as `1/1`).
fn parse_chunk_body(body: &str) -> Result<ParsedChunk> {
    if let Some((header, payload)) = body.split_once(':') {
        if let Some((i_str, n_str)) = header.split_once('/') {
            if let (Ok(i), Ok(n)) = (i_str.parse::<usize>(), n_str.parse::<usize>()) {
                if i == 0 || n == 0 || i > n {
                    anyhow::bail!(
                        "invalid embedded Kraken report chunk header: {header:?} \
                         (index/total out of range)"
                    );
                }
                return Ok(ParsedChunk {
                    index: i,
                    total: n,
                    b64: payload.to_owned(),
                });
            }
        }
    }
    Ok(ParsedChunk {
        index: 1,
        total: 1,
        b64: body.to_owned(),
    })
}

/// Scan `@CO` lines in `header` for an embedded Kraken report.
///
/// Returns `None` if no `krak:report:` comment is present.
/// Returns `Err` if matching comments are found but base64 or JSON decoding fails,
/// or if more than one embedded report is found (indicating the file was annotated
/// more than once).
pub(crate) fn entries_from_header(header: &sam::Header) -> Result<Option<Vec<KrakenReportEntry>>> {
    let mut chunks: Vec<ParsedChunk> = Vec::new();
    let mut parse_err: Option<anyhow::Error> = None;

    for comment in header.comments() {
        // BString derefs through BStr to [u8]; skip any non-UTF-8 comment lines.
        if let Ok(s) = std::str::from_utf8(comment) {
            if let Some(body) = s.strip_prefix(HEADER_COMMENT_PREFIX) {
                match parse_chunk_body(body) {
                    Ok(c) => chunks.push(c),
                    Err(e) if parse_err.is_none() => parse_err = Some(e),
                    Err(_) => {}
                }
            }
        }
    }

    if chunks.is_empty() && parse_err.is_none() {
        return Ok(None);
    }

    // Group by total to detect duplicate annotations regardless of decode validity.
    let mut by_total: std::collections::BTreeMap<usize, Vec<ParsedChunk>> =
        std::collections::BTreeMap::new();
    for c in chunks {
        by_total.entry(c.total).or_default().push(c);
    }

    if by_total.len() > 1 {
        anyhow::bail!(
            "found embedded Kraken report chunks with mismatched totals; \
             this file may have been annotated more than once. Re-annotate \
             from the original unannotated file."
        );
    }

    let Some((total, mut group)) = by_total.into_iter().next() else {
        // No valid chunks but we did see a malformed header line.
        return Err(parse_err.unwrap_or_else(|| {
            anyhow::anyhow!("failed to parse embedded Kraken report comments")
        }));
    };

    group.sort_by_key(|c| c.index);

    // Detect duplicate annotation: more chunks than the declared total means
    // at least one chunk index repeats.
    if group.len() > total {
        anyhow::bail!(
            "found {} embedded Kraken report chunks but total is {total}; \
             the SAM/BAM/CRAM file appears to have been annotated more than \
             once. Re-annotate from the original unannotated file.",
            group.len()
        );
    }

    if group.len() < total {
        anyhow::bail!(
            "embedded Kraken report is missing chunks: have {} of {total}",
            group.len()
        );
    }

    // Verify chunk indices are exactly 1..=total with no duplicates.
    for (expected, c) in (1..=total).zip(group.iter()) {
        if c.index != expected {
            // Duplicate index means the same chunk appeared twice -> annotated
            // more than once.
            if group.iter().filter(|x| x.index == c.index).count() > 1 {
                anyhow::bail!(
                    "duplicate chunk index {} in embedded Kraken report; \
                     the SAM/BAM/CRAM file appears to have been annotated \
                     more than once. Re-annotate from the original \
                     unannotated file.",
                    c.index
                );
            }
            anyhow::bail!(
                "embedded Kraken report has out-of-order or missing chunk: \
                 expected index {expected}, got {}",
                c.index
            );
        }
    }

    let mut joined = String::new();
    for c in &group {
        joined.push_str(&c.b64);
    }

    let bytes = STANDARD
        .decode(&joined)
        .context("failed to base64-decode embedded Kraken report")?;
    let entries: Vec<KrakenReportEntry> =
        serde_json::from_slice(&bytes).context("failed to JSON-decode embedded Kraken report")?;
    Ok(Some(entries))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam;

    fn make_entry(taxon_id: u32, name: &str, indent: usize) -> KrakenReportEntry {
        KrakenReportEntry {
            pct_fragments: 50.0,
            num_fragments_clade: 100,
            num_fragments_direct: 10,
            rank_code: "S".to_owned(),
            taxon_id,
            name: name.to_owned(),
            indent,
            minimizer_count: None,
            distinct_minimizer_count: None,
        }
    }

    #[test]
    fn test_round_trip_recovers_entries() {
        let entries = vec![
            make_entry(1, "root", 0),
            make_entry(9606, "Homo sapiens", 4),
        ];
        let comments = entries_to_header_comment(&entries).unwrap();
        assert!(!comments.is_empty());
        for c in &comments {
            assert!(c.starts_with(HEADER_COMMENT_PREFIX));
        }

        let mut header = sam::Header::default();
        for c in comments {
            header.add_comment(c);
        }

        let recovered = entries_from_header(&header).unwrap().unwrap();
        assert_eq!(recovered.len(), 2);
        assert_eq!(recovered[0].taxon_id, 1);
        assert_eq!(recovered[1].taxon_id, 9606);
        assert_eq!(recovered[1].name, "Homo sapiens");
        assert_eq!(recovered[1].indent, 4);
    }

    #[test]
    fn test_round_trip_forces_multi_chunk() {
        // Build enough entries that the JSON+base64 payload exceeds the
        // test-only CHUNK_SIZE and is split into multiple comments.
        let entries: Vec<_> = (0..50)
            .map(|i| make_entry(i, &format!("taxon_{i}"), 0))
            .collect();
        let comments = entries_to_header_comment(&entries).unwrap();
        assert!(
            comments.len() > 1,
            "expected multi-chunk encoding, got {}",
            comments.len()
        );
        // All chunks share the same total.
        let total: usize = comments.len();
        for (i, c) in comments.iter().enumerate() {
            let body = c.strip_prefix(HEADER_COMMENT_PREFIX).unwrap();
            let header = body.split(':').next().unwrap();
            assert_eq!(header, format!("{}/{total}", i + 1));
        }

        let mut header = sam::Header::default();
        for c in comments {
            header.add_comment(c);
        }
        let recovered = entries_from_header(&header).unwrap().unwrap();
        assert_eq!(recovered.len(), entries.len());
        assert_eq!(recovered[10].taxon_id, 10);
    }

    #[test]
    fn test_returns_none_when_no_comment() {
        let header = sam::Header::default();
        assert!(entries_from_header(&header).unwrap().is_none());
    }

    #[test]
    fn test_returns_none_when_other_comments_present() {
        let mut header = sam::Header::default();
        header.add_comment("created by some-other-tool".to_owned());
        assert!(entries_from_header(&header).unwrap().is_none());
    }

    #[test]
    fn test_returns_err_on_corrupt_base64() {
        let mut header = sam::Header::default();
        header.add_comment("krak:report:1/1:!!!not-valid-base64!!!".to_owned());
        assert!(entries_from_header(&header).is_err());
    }

    #[test]
    fn test_legacy_single_comment_format_is_accepted() {
        // Legacy form with no <i>/<n>: prefix.
        let entries = vec![make_entry(1, "root", 0)];
        let json = serde_json::to_string(&entries).unwrap();
        let b64 = STANDARD.encode(json.as_bytes());
        let mut header = sam::Header::default();
        header.add_comment(format!("{HEADER_COMMENT_PREFIX}{b64}"));
        let recovered = entries_from_header(&header).unwrap().unwrap();
        assert_eq!(recovered.len(), 1);
        assert_eq!(recovered[0].taxon_id, 1);
    }

    #[test]
    fn test_errors_on_multiple_reports_legacy_form() {
        let entries = vec![make_entry(1, "root", 0)];
        let json = serde_json::to_string(&entries).unwrap();
        let b64 = STANDARD.encode(json.as_bytes());
        let comment = format!("{HEADER_COMMENT_PREFIX}{b64}");

        let mut header = sam::Header::default();
        header.add_comment(comment.clone());
        header.add_comment(comment);

        let result = entries_from_header(&header);
        assert!(result.is_err());
        let msg = format!("{:#}", result.unwrap_err());
        assert!(msg.contains("annotated more than once") || msg.contains("more than once"));
    }

    #[test]
    fn test_errors_on_corrupt_plus_valid_picks_duplicate_error() {
        // Valid encoding (any number of chunks N) plus a corrupt comment that
        // declares total=N+1 introduces a second "total" group, which the
        // reader treats as more-than-one report regardless of whether the
        // corrupt payload would decode.
        let entries = vec![make_entry(1, "root", 0)];
        let comments = entries_to_header_comment(&entries).unwrap();
        let n = comments.len();
        let foreign_total = n + 1;

        let mut header = sam::Header::default();
        for c in comments {
            header.add_comment(c);
        }
        header.add_comment(format!("krak:report:1/{foreign_total}:!!!corrupt!!!"));

        let result = entries_from_header(&header);
        assert!(result.is_err());
        let msg = format!("{:#}", result.unwrap_err());
        assert!(
            msg.contains("annotated more than once")
                || msg.contains("more than once")
                || msg.contains("mismatched totals"),
            "expected duplicate-annotation error, got: {msg}"
        );
    }

    #[test]
    fn test_parse_chunk_body_rejects_zero_index() {
        // Index/total of 0 (or i > n) is structurally invalid.
        let err = parse_chunk_body("0/3:abc").unwrap_err();
        assert!(format!("{err:#}").contains("out of range"));
    }

    #[test]
    fn test_parse_chunk_body_rejects_index_greater_than_total() {
        let err = parse_chunk_body("4/3:abc").unwrap_err();
        assert!(format!("{err:#}").contains("out of range"));
    }

    #[test]
    fn test_parse_chunk_body_legacy_form_when_header_unparseable() {
        // A body whose pre-colon portion is not `<digits>/<digits>` falls
        // back to legacy single-chunk form rather than erroring.
        let p = parse_chunk_body("abc").unwrap();
        assert_eq!(p.index, 1);
        assert_eq!(p.total, 1);
        assert_eq!(p.b64, "abc");
    }

    #[test]
    fn test_only_malformed_comment_surfaces_parse_error() {
        // A header with a single malformed `krak:report:` comment (and no
        // valid chunks anywhere) must surface the parse error rather than
        // silently returning Ok(None).
        let mut header = sam::Header::default();
        header.add_comment("krak:report:0/0:bogus".to_owned());
        let err = entries_from_header(&header).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("out of range") || msg.contains("invalid"),
            "expected parse-error surface, got: {msg}"
        );
    }

    #[test]
    fn test_duplicate_chunk_index_detected() {
        // Two chunks both claim index 1/2; the third claim 2/2. The reader
        // sees 3 chunks but total=2 → duplicate-annotation error.
        let mut header = sam::Header::default();
        header.add_comment("krak:report:1/2:AAAA".to_owned());
        header.add_comment("krak:report:1/2:BBBB".to_owned());
        header.add_comment("krak:report:2/2:CCCC".to_owned());
        let err = entries_from_header(&header).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("annotated more than once") || msg.contains("more than once"),
            "expected duplicate-annotation error, got: {msg}"
        );
    }

    #[test]
    fn test_duplicate_chunk_index_when_len_equals_total() {
        // total=3, chunks at indices [1, 1, 2]: group.len()==total but the
        // index sequence has a duplicate at slot 2. Exercises the duplicate-
        // detection inside the index-walk loop (lines ~178).
        let mut header = sam::Header::default();
        header.add_comment("krak:report:1/3:AAAA".to_owned());
        header.add_comment("krak:report:1/3:BBBB".to_owned());
        header.add_comment("krak:report:2/3:CCCC".to_owned());
        let err = entries_from_header(&header).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("annotated more than once")
                || msg.contains("more than once")
                || msg.contains("duplicate chunk index"),
            "expected duplicate-annotation error, got: {msg}"
        );
    }

    #[test]
    fn test_two_malformed_comments_keeps_first_error() {
        // First malformed comment populates parse_err. The second malformed
        // comment exercises the `Err(_) => {}` arm (second-parse-error
        // swallow). With no valid chunks, the surfaced error is the first.
        let mut header = sam::Header::default();
        header.add_comment("krak:report:0/0:first-bad".to_owned());
        header.add_comment("krak:report:0/0:second-bad".to_owned());
        let err = entries_from_header(&header).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("out of range") || msg.contains("invalid"));
    }

    #[test]
    fn test_errors_on_mismatched_totals() {
        let mut header = sam::Header::default();
        header.add_comment("krak:report:1/2:AAAA".to_owned());
        header.add_comment("krak:report:1/3:AAAA".to_owned());
        let result = entries_from_header(&header);
        assert!(result.is_err());
        let msg = format!("{:#}", result.unwrap_err());
        assert!(msg.contains("mismatched totals") || msg.contains("more than once"));
    }

    #[test]
    fn test_errors_on_missing_index() {
        // Declare 2 chunks but supply only chunk 1.
        let mut header = sam::Header::default();
        header.add_comment("krak:report:1/2:AAAA".to_owned());
        let result = entries_from_header(&header);
        assert!(result.is_err());
        let msg = format!("{:#}", result.unwrap_err());
        assert!(msg.contains("missing chunks"));
    }
}
