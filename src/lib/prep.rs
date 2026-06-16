//! Convert FASTX/SAM/BAM/CRAM for Kraken classification.

use std::collections::VecDeque;
use std::io::Write;

use anyhow::{Context, Result};

/// Arguments for the `prep` command.
pub struct PrepArgs {
    /// Primary input file (FASTA, FASTQ, SAM, BAM, or CRAM).
    pub input: std::path::PathBuf,
    /// R2 input file. Only valid with FASTQ/FASTA primary input; mutually
    /// exclusive with `per_record`.
    pub input2: Option<std::path::PathBuf>,
    /// Disable auto pair-detection: emit each FASTX record (or SAM/BAM/CRAM
    /// primary record) as its own single-end template.
    pub per_record: bool,
    /// Output FASTA file.
    pub output: std::path::PathBuf,
    /// Optional reference FASTA for CRAM decompression (requires `.fai` index).
    pub cram_reference: Option<std::path::PathBuf>,
}

/// One DNA template: a name, a required R1 sequence, and an optional R2 sequence.
#[derive(Debug)]
pub(crate) struct Template {
    pub name: String,
    pub r1: Vec<u8>,
    pub r2: Option<Vec<u8>>,
}

/// Pair suffix observed on a raw read name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairSuffix {
    None,
    Slash1,
    Slash2,
}

impl PairSuffix {
    fn detect(name: &str) -> Self {
        if name.ends_with("/1") {
            Self::Slash1
        } else if name.ends_with("/2") {
            Self::Slash2
        } else {
            Self::None
        }
    }
}

/// Validate that two raw mate names carry compatible pair suffixes.
///
/// Accepts only `(None, None)` (Casava-style mates that already share a token)
/// or `(Slash1, Slash2)` (Kraken v1 / classic Illumina). Anything else; two
/// `/1`s, two `/2`s, mixed presence, or `/2` then `/1`; indicates the pair is
/// out of order or duplicated.
fn check_pair_suffixes(name1_raw: &str, name2_raw: &str) -> anyhow::Result<()> {
    let s1 = PairSuffix::detect(name1_raw);
    let s2 = PairSuffix::detect(name2_raw);
    match (s1, s2) {
        (PairSuffix::None, PairSuffix::None) | (PairSuffix::Slash1, PairSuffix::Slash2) => Ok(()),
        _ => Err(anyhow::anyhow!(
            "mismatched pair suffixes: R1 {:?} and R2 {:?}; \
             expected either both with no /1 /2 suffix or R1 with /1 and R2 with /2. \
             Use --per-record to disable interleaved/paired detection.",
            name1_raw,
            name2_raw
        )),
    }
}

/// Reverse-complement a DNA sequence in place.
///
/// Maps `A<->T`, `C<->G`, `N<->N`. Bytes outside that alphabet pass through
/// unchanged (so e.g. IUPAC ambiguity codes are preserved as-is rather than
/// silently corrupted).
fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        };
    }
}

/// Parse a FASTQ record's name as the first whitespace-delimited token, owned.
///
/// Errors on non-UTF-8 bytes. Does not strip pair suffixes; callers apply
/// `strip_pair_suffix` when comparing mate names.
fn parse_fastq_name(bytes: &[u8]) -> anyhow::Result<String> {
    let tok = bytes.split(|&b| b == b' ').next().unwrap_or(bytes);
    std::str::from_utf8(tok).map(|s| s.to_owned()).map_err(|_| {
        anyhow::anyhow!(
            "non-UTF-8 FASTQ read name: {}",
            String::from_utf8_lossy(tok)
        )
    })
}

/// Validate paired-record names produced by parsing R1 and R2 names. Returns
/// the common base name (with pair suffixes stripped) on success.
///
/// Errors if the suffixes disagree (`check_pair_suffixes`), if either name is
/// empty after stripping, or if the names don't match. `fmt_label` is `FASTQ`
/// or `FASTA` and appears in the suffix-conflict context and the empty-name
/// errors so anyhow stack traces disambiguate the format.
fn validate_pair_names(
    raw1: &str,
    raw2: &str,
    pair_num: usize,
    fmt_label: &str,
) -> anyhow::Result<String> {
    check_pair_suffixes(raw1, raw2)
        .map_err(|e| e.context(format!("paired {fmt_label} at pair {pair_num}")))?;
    let name1 = crate::strip_pair_suffix(raw1).to_owned();
    if name1.is_empty() {
        anyhow::bail!("R1 {fmt_label} record has empty name at pair {pair_num}");
    }
    let name2 = crate::strip_pair_suffix(raw2);
    if name2.is_empty() {
        anyhow::bail!("R2 {fmt_label} record has empty name at pair {pair_num}");
    }
    if name1 != name2 {
        anyhow::bail!("R1 name {name1:?} does not match R2 name {name2:?} (at pair {pair_num})");
    }
    Ok(name1)
}

fn iter_fastq_auto(
    path: &std::path::Path,
    per_record: bool,
) -> anyhow::Result<impl Iterator<Item = anyhow::Result<Template>>> {
    let reader = crate::open_fastx_reader(path)
        .with_context(|| format!("failed to open FASTQ: {}", path.display()))?;
    Ok(iter_fastq_auto_from_reader(reader, per_record))
}

fn iter_fasta_auto(
    path: &std::path::Path,
    per_record: bool,
) -> anyhow::Result<impl Iterator<Item = anyhow::Result<Template>>> {
    let reader = crate::open_fastx_reader(path)
        .with_context(|| format!("failed to open FASTA: {}", path.display()))?;
    Ok(iter_fasta_auto_from_reader(reader, per_record))
}

/// Stream FASTA records as templates, auto-detecting interleaved layout from
/// the first one or two records. Detection rule mirrors
/// `iter_fastq_auto_from_reader`.
fn iter_fasta_auto_from_reader<R: std::io::BufRead + 'static>(
    reader: R,
    per_record: bool,
) -> impl Iterator<Item = anyhow::Result<Template>> {
    let mut reader = noodles::fasta::io::Reader::new(reader);
    let mut detection_done = false;
    let mut interleaved = false;
    let mut buffered: VecDeque<(String, Vec<u8>)> = VecDeque::with_capacity(2);
    let mut pair_num = 0usize;

    fn read_one(
        reader: &mut noodles::fasta::io::Reader<impl std::io::BufRead>,
    ) -> anyhow::Result<Option<(String, Vec<u8>)>> {
        let mut def = String::new();
        let mut seq = Vec::new();
        match reader.read_definition(&mut def)? {
            0 => Ok(None),
            _ => {
                reader.read_sequence(&mut seq)?;
                let raw = def
                    .trim()
                    .trim_start_matches('>')
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_owned();
                Ok(Some((raw, seq)))
            }
        }
    }

    std::iter::from_fn(move || {
        if !detection_done {
            detection_done = true;
            let first = match read_one(&mut reader) {
                Ok(Some(rec)) => rec,
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            };
            if !per_record {
                let first_name = &first.0;
                if first_name.is_empty() {
                    return Some(Err(anyhow::anyhow!("FASTA record has empty name")));
                }
                if first_name.ends_with("/2") {
                    return Some(Err(anyhow::anyhow!(
                        "interleaved FASTA: first record '{}' ends with /2; mates are out of order. \
                         Use --per-record to disable interleaved detection.",
                        first_name
                    )));
                }
                if first_name.ends_with("/1") {
                    interleaved = true;
                    buffered.push_back(first);
                } else {
                    match read_one(&mut reader) {
                        Ok(Some(second)) => {
                            if !second.0.is_empty() && first.0 == second.0 {
                                interleaved = true;
                            }
                            buffered.push_back(first);
                            buffered.push_back(second);
                        }
                        Ok(None) => {
                            buffered.push_back(first);
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
            } else {
                buffered.push_back(first);
            }
        }

        if interleaved {
            pair_num += 1;
            let r1 = match buffered.pop_front() {
                Some(r) => r,
                None => match read_one(&mut reader) {
                    Ok(Some(r)) => r,
                    Ok(None) => return None,
                    Err(e) => return Some(Err(e)),
                },
            };
            let r2 = match buffered.pop_front() {
                Some(r) => r,
                None => match read_one(&mut reader) {
                    Ok(Some(r)) => r,
                    Ok(None) => {
                        return Some(Err(anyhow::anyhow!(
                            "interleaved FASTA has an odd number of records (pair {} has no R2)",
                            pair_num
                        )))
                    }
                    Err(e) => return Some(Err(e)),
                },
            };
            if let Err(e) = check_pair_suffixes(&r1.0, &r2.0) {
                return Some(Err(
                    e.context(format!("interleaved FASTA at pair {}", pair_num))
                ));
            }
            let name1 = crate::strip_pair_suffix(&r1.0).to_owned();
            if name1.is_empty() {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTA R1 has empty name at pair {}",
                    pair_num
                )));
            }
            let name2 = crate::strip_pair_suffix(&r2.0).to_owned();
            if name2.is_empty() {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTA R2 has empty name at pair {}",
                    pair_num
                )));
            }
            if name1 != name2 {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTA: consecutive records have different names: {:?} then {:?} (at pair {}). \
                     Use --per-record to disable interleaved detection.",
                    name1,
                    name2,
                    pair_num
                )));
            }
            return Some(Ok(Template {
                name: name1,
                r1: r1.1,
                r2: Some(r2.1),
            }));
        }

        let r = match buffered.pop_front() {
            Some(r) => r,
            None => match read_one(&mut reader) {
                Ok(Some(r)) => r,
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            },
        };
        if r.0.is_empty() {
            return Some(Err(anyhow::anyhow!("FASTA record has empty name")));
        }
        Some(Ok(Template {
            name: r.0,
            r1: r.1,
            r2: None,
        }))
    })
}

/// Stream FASTQ records as templates, auto-detecting interleaved layout
/// from the first one or two records.
///
/// Detection rule (skipped when `per_record` is `true`):
/// 1. First parsed name ends with `/2` -> error (out of order).
/// 2. First parsed name ends with `/1` -> interleaved.
/// 3. Otherwise peek the second record; if its parsed name equals the first -> interleaved.
/// 4. Otherwise -> single-end.
///
/// `per_record` short-circuits to single-end without reading a second record.
fn iter_fastq_auto_from_reader<R: std::io::BufRead + 'static>(
    reader: R,
    per_record: bool,
) -> impl Iterator<Item = anyhow::Result<Template>> {
    let mut reader = noodles::fastq::io::Reader::new(reader);
    let mut detection_done = false;
    let mut interleaved = false;
    let mut buffered: VecDeque<noodles::fastq::Record> = VecDeque::with_capacity(2);
    let mut pair_num = 0usize;

    std::iter::from_fn(move || {
        if !detection_done {
            detection_done = true;
            let mut first = noodles::fastq::Record::default();
            match reader.read_record(&mut first) {
                Ok(0) => return None,
                Err(e) => return Some(Err(anyhow::Error::from(e))),
                Ok(_) => {}
            }
            if !per_record {
                let first_name = match parse_fastq_name(first.name()) {
                    Ok(n) => n,
                    Err(e) => return Some(Err(e)),
                };
                if first_name.is_empty() {
                    return Some(Err(anyhow::anyhow!("FASTQ record has empty name")));
                }
                if first_name.ends_with("/2") {
                    return Some(Err(anyhow::anyhow!(
                        "interleaved FASTQ: first record '{}' ends with /2; mates are out of order. \
                         Use --per-record to disable interleaved detection.",
                        first_name
                    )));
                }
                if first_name.ends_with("/1") {
                    interleaved = true;
                    buffered.push_back(first);
                } else {
                    let mut second = noodles::fastq::Record::default();
                    match reader.read_record(&mut second) {
                        Ok(0) => {
                            buffered.push_back(first);
                        }
                        Err(e) => return Some(Err(anyhow::Error::from(e))),
                        Ok(_) => {
                            let second_name = match parse_fastq_name(second.name()) {
                                Ok(n) => n,
                                Err(e) => return Some(Err(e)),
                            };
                            if !second_name.is_empty() && first_name == second_name {
                                interleaved = true;
                            }
                            buffered.push_back(first);
                            buffered.push_back(second);
                        }
                    }
                }
            } else {
                buffered.push_back(first);
            }
        }

        if interleaved {
            pair_num += 1;
            let r1 = match buffered.pop_front() {
                Some(r) => r,
                None => {
                    let mut r = noodles::fastq::Record::default();
                    match reader.read_record(&mut r) {
                        Ok(0) => return None,
                        Err(e) => return Some(Err(anyhow::Error::from(e))),
                        Ok(_) => r,
                    }
                }
            };
            let r2 = match buffered.pop_front() {
                Some(r) => r,
                None => {
                    let mut r = noodles::fastq::Record::default();
                    match reader.read_record(&mut r) {
                        Ok(0) => {
                            return Some(Err(anyhow::anyhow!(
                            "interleaved FASTQ has an odd number of records (pair {} has no R2)",
                            pair_num
                        )))
                        }
                        Err(e) => return Some(Err(anyhow::Error::from(e))),
                        Ok(_) => r,
                    }
                }
            };
            let raw1 = match parse_fastq_name(r1.name()) {
                Ok(n) => n,
                Err(e) => return Some(Err(e)),
            };
            let raw2 = match parse_fastq_name(r2.name()) {
                Ok(n) => n,
                Err(e) => return Some(Err(e)),
            };
            if let Err(e) = check_pair_suffixes(&raw1, &raw2) {
                return Some(Err(
                    e.context(format!("interleaved FASTQ at pair {}", pair_num))
                ));
            }
            let name1 = crate::strip_pair_suffix(&raw1).to_owned();
            if name1.is_empty() {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTQ R1 has empty name at pair {}",
                    pair_num
                )));
            }
            let name2 = crate::strip_pair_suffix(&raw2).to_owned();
            if name2.is_empty() {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTQ R2 has empty name at pair {}",
                    pair_num
                )));
            }
            if name1 != name2 {
                return Some(Err(anyhow::anyhow!(
                    "interleaved FASTQ: consecutive records have different names: {:?} then {:?} (at pair {}). \
                     Use --per-record to disable interleaved detection.",
                    name1,
                    name2,
                    pair_num
                )));
            }
            return Some(Ok(Template {
                name: name1,
                r1: r1.sequence().to_vec(),
                r2: Some(r2.sequence().to_vec()),
            }));
        }

        let r = match buffered.pop_front() {
            Some(r) => r,
            None => {
                let mut r = noodles::fastq::Record::default();
                match reader.read_record(&mut r) {
                    Ok(0) => return None,
                    Err(e) => return Some(Err(anyhow::Error::from(e))),
                    Ok(_) => r,
                }
            }
        };
        let name = match parse_fastq_name(r.name()) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        if name.is_empty() {
            return Some(Err(anyhow::anyhow!("FASTQ record has empty name")));
        }
        Some(Ok(Template {
            name,
            r1: r.sequence().to_vec(),
            r2: None,
        }))
    })
}

/// Stream paired FASTQ records from two readers (R1 and R2).
///
/// Errors if names do not match at any position or if the files have different lengths.
fn iter_paired_fastq_from_readers<R1, R2>(
    r1_reader: R1,
    r2_reader: R2,
) -> impl Iterator<Item = anyhow::Result<Template>>
where
    R1: std::io::BufRead + 'static,
    R2: std::io::BufRead + 'static,
{
    let mut r1_reader = noodles::fastq::io::Reader::new(r1_reader);
    let mut r2_reader = noodles::fastq::io::Reader::new(r2_reader);
    let mut pair_num = 0usize;
    std::iter::from_fn(move || {
        pair_num += 1;
        let mut r1 = noodles::fastq::Record::default();
        let mut r2 = noodles::fastq::Record::default();
        let r1_n = match r1_reader.read_record(&mut r1) {
            Ok(n) => n,
            Err(e) => return Some(Err(e.into())),
        };
        let r2_n = match r2_reader.read_record(&mut r2) {
            Ok(n) => n,
            Err(e) => return Some(Err(e.into())),
        };
        match (r1_n, r2_n) {
            (0, 0) => return None,
            (0, _) | (_, 0) => {
                return Some(Err(anyhow::anyhow!(
                    "R1 and R2 FASTQ files have unequal record counts (at pair {})",
                    pair_num
                )))
            }
            _ => {}
        }
        let parse_raw = |rec: &noodles::fastq::Record| {
            let bytes = rec.name();
            let tok = bytes.split(|&b| b == b' ').next().unwrap_or(bytes);
            std::str::from_utf8(tok)
                .map(|s| s.to_owned())
                .map_err(|_| anyhow::anyhow!("non-UTF-8 FASTQ read name at pair {}", pair_num))
        };
        let raw1 = match parse_raw(&r1) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        let raw2 = match parse_raw(&r2) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        let name = match validate_pair_names(&raw1, &raw2, pair_num, "FASTQ") {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        Some(Ok(Template {
            name,
            r1: r1.sequence().to_vec(),
            r2: Some(r2.sequence().to_vec()),
        }))
    })
}

fn iter_paired_fastq(
    r1_path: &std::path::Path,
    r2_path: &std::path::Path,
) -> anyhow::Result<impl Iterator<Item = anyhow::Result<Template>>> {
    let r1 = crate::open_fastx_reader(r1_path)
        .with_context(|| format!("failed to open R1 FASTQ: {}", r1_path.display()))?;
    let r2 = crate::open_fastx_reader(r2_path)
        .with_context(|| format!("failed to open R2 FASTQ: {}", r2_path.display()))?;
    Ok(iter_paired_fastq_from_readers(r1, r2))
}

/// Stream paired FASTA records by zipping two readers.
fn iter_paired_fasta_from_readers<R1, R2>(
    reader1: R1,
    reader2: R2,
) -> impl Iterator<Item = anyhow::Result<Template>>
where
    R1: std::io::BufRead + 'static,
    R2: std::io::BufRead + 'static,
{
    let mut r1 = noodles::fasta::io::Reader::new(reader1);
    let mut r2 = noodles::fasta::io::Reader::new(reader2);
    let mut def1 = String::new();
    let mut seq1 = Vec::new();
    let mut def2 = String::new();
    let mut seq2 = Vec::new();
    let mut pair_num = 0usize;
    std::iter::from_fn(move || {
        pair_num += 1;
        def1.clear();
        seq1.clear();
        def2.clear();
        seq2.clear();
        let n1 = match r1.read_definition(&mut def1) {
            Ok(n) => n,
            Err(e) => return Some(Err(e.into())),
        };
        if n1 > 0 {
            if let Err(e) = r1.read_sequence(&mut seq1) {
                return Some(Err(e.into()));
            }
        }
        let n2 = match r2.read_definition(&mut def2) {
            Ok(n) => n,
            Err(e) => return Some(Err(e.into())),
        };
        if n2 > 0 {
            if let Err(e) = r2.read_sequence(&mut seq2) {
                return Some(Err(e.into()));
            }
        }
        match (n1, n2) {
            (0, 0) => None,
            (0, _) | (_, 0) => Some(Err(anyhow::anyhow!(
                "R1 and R2 FASTA files have unequal record counts (at pair {})",
                pair_num
            ))),
            _ => {
                let parse_raw = |def: &str| {
                    def.trim()
                        .trim_start_matches('>')
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_owned()
                };
                let raw1 = parse_raw(&def1);
                let raw2 = parse_raw(&def2);
                let name = match validate_pair_names(&raw1, &raw2, pair_num, "FASTA") {
                    Ok(n) => n,
                    Err(e) => return Some(Err(e)),
                };
                Some(Ok(Template {
                    name,
                    r1: std::mem::take(&mut seq1),
                    r2: Some(std::mem::take(&mut seq2)),
                }))
            }
        }
    })
}

fn iter_paired_fasta(
    path1: &std::path::Path,
    path2: &std::path::Path,
) -> anyhow::Result<impl Iterator<Item = anyhow::Result<Template>>> {
    let r1 = crate::open_fastx_reader(path1)
        .with_context(|| format!("failed to open R1 FASTA: {}", path1.display()))?;
    let r2 = crate::open_fastx_reader(path2)
        .with_context(|| format!("failed to open R2 FASTA: {}", path2.display()))?;
    Ok(iter_paired_fasta_from_readers(r1, r2))
}

/// Extract a UTF-8 read name from a `RecordBuf`, or error.
///
/// Errors when the name field is missing, empty, or the SAM placeholder `*`.
fn record_name(record: &noodles::sam::alignment::record_buf::RecordBuf) -> Result<String> {
    let bytes = record
        .name()
        .ok_or_else(|| anyhow::anyhow!("alignment record has no name"))?;
    if bytes.is_empty() || bytes == b"*" {
        anyhow::bail!("alignment record has empty or placeholder name");
    }
    std::str::from_utf8(bytes)
        .context("alignment record name is not valid UTF-8")
        .map(|s| s.to_owned())
}

/// Stream single-end alignment records as templates.
///
/// Secondary (0x100) and supplementary (0x800) records are always skipped;
/// `per_record` controls only whether paired-flag inspection happens.
///
/// When `per_record` is `false`, errors on any primary record with the paired
/// flag (0x1) set; callers should route to the query-grouped iterator
/// instead. When `per_record` is `true`, paired records are emitted as
/// independent single-end templates; their QNAME is suffixed with `/1` or
/// `/2` (per the FIRST_SEGMENT / LAST_SEGMENT flags) so each emitted record
/// has a unique name. Names that already end in `/1` or `/2` pass through
/// unchanged.
fn iter_single_end_alignment_templates<I>(
    records: I,
    per_record: bool,
) -> impl Iterator<Item = anyhow::Result<Template>>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
{
    records.filter_map(move |result| {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                return Some(Err(
                    anyhow::Error::from(e).context("failed to read alignment record")
                ))
            }
        };
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            return None;
        }
        if flags.is_segmented() && !per_record {
            let name = match record_name(&record) {
                Ok(n) => n,
                Err(e) => return Some(Err(e)),
            };
            return Some(Err(anyhow::anyhow!(
                "single-end mode: record {name:?} has the paired flag (0x1); \
                 queryname-sort the file first (samtools sort -n), then re-run \
                 (queryname-sorted paired SAM/BAM is handled automatically), \
                 or pass --per-record to emit each record independently"
            )));
        }
        let mut name = match record_name(&record) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        let already_suffixed = name.ends_with("/1") || name.ends_with("/2");
        if flags.is_segmented() && !already_suffixed {
            if flags.is_first_segment() {
                name.push_str("/1");
            } else if flags.is_last_segment() {
                name.push_str("/2");
            }
        }
        let mut seq = record.sequence().as_ref().to_vec();
        if flags.is_reverse_complemented() {
            reverse_complement(&mut seq);
        }
        Some(Ok(Template {
            name,
            r1: seq,
            r2: None,
        }))
    })
}

/// Stream query-grouped alignment records as templates.
///
/// Consecutive records sharing a QNAME form one group. Secondary (0x100) and
/// supplementary (0x800) records are skipped. R1 is assigned by first-segment
/// flag (0x40); R2 by last-segment flag (0x80).
fn iter_query_grouped_alignment_templates<I>(
    mut records: I,
) -> impl Iterator<Item = anyhow::Result<Template>>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
{
    let mut lookahead: Option<noodles::sam::alignment::record_buf::RecordBuf> = None;
    let mut done = false;

    std::iter::from_fn(move || {
        if done {
            return None;
        }

        // Find the first primary record to start this group (use lookahead if available)
        let start = if let Some(r) = lookahead.take() {
            r
        } else {
            loop {
                match records.next() {
                    None => {
                        done = true;
                        return None;
                    }
                    Some(Err(e)) => return Some(Err(e.into())),
                    Some(Ok(r)) => {
                        let f = r.flags();
                        if !f.is_secondary() && !f.is_supplementary() {
                            break r;
                        }
                    }
                }
            }
        };

        let group_name = match record_name(&start) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };

        let mut r1: Option<Vec<u8>> = None;
        let mut r2: Option<Vec<u8>> = None;

        // Assign start record to r1 or r2
        let flags = start.flags();
        // TODO: borrowed-Vec refactor; these per-record allocations could
        // share a reusable scratch buffer threaded through the iterator, at
        // the cost of an extra borrow lifetime in the public signature.
        let mut seq = start.sequence().as_ref().to_vec();
        if flags.is_reverse_complemented() {
            reverse_complement(&mut seq);
        }
        if flags.is_first_segment() {
            r1 = Some(seq);
        } else if flags.is_last_segment() {
            r2 = Some(seq);
        } else {
            r1 = Some(seq);
        }

        // Accumulate remaining records in this group
        loop {
            let record = loop {
                match records.next() {
                    None => {
                        done = true;
                        return match (r1.take(), r2.take()) {
                            (None, None) => None,
                            (Some(s), r2_opt) => Some(Ok(Template {
                                name: group_name,
                                r1: s,
                                r2: r2_opt,
                            })),
                            (None, Some(_)) => Some(Err(anyhow::anyhow!(
                                "R2 record found without R1 for template {:?}",
                                group_name
                            ))),
                        };
                    }
                    Some(Err(e)) => return Some(Err(e.into())),
                    Some(Ok(r)) => {
                        let f = r.flags();
                        if !f.is_secondary() && !f.is_supplementary() {
                            break r;
                        }
                    }
                }
            };

            let rec_name = match record_name(&record) {
                Ok(n) => n,
                Err(e) => return Some(Err(e)),
            };

            if rec_name != group_name {
                // New group; buffer this record and return the completed group
                lookahead = Some(record);
                return match (r1.take(), r2.take()) {
                    (None, None) => None, // shouldn't happen (start was primary)
                    (Some(s), r2_opt) => Some(Ok(Template {
                        name: group_name,
                        r1: s,
                        r2: r2_opt,
                    })),
                    (None, Some(_)) => Some(Err(anyhow::anyhow!(
                        "R2 record found without R1 for template {:?}",
                        group_name
                    ))),
                };
            }

            // Same group; assign to r1 or r2
            let flags = record.flags();
            let mut seq = record.sequence().as_ref().to_vec();
            if flags.is_reverse_complemented() {
                reverse_complement(&mut seq);
            }
            if flags.is_first_segment() {
                if r1.is_some() {
                    return Some(Err(anyhow::anyhow!(
                        "two primary R1 records for template {:?}",
                        rec_name
                    )));
                }
                r1 = Some(seq);
            } else if flags.is_last_segment() {
                if r2.is_some() {
                    return Some(Err(anyhow::anyhow!(
                        "two primary R2 records for template {:?}",
                        rec_name
                    )));
                }
                r2 = Some(seq);
            } else if r1.is_none() {
                r1 = Some(seq);
            } else if r2.is_none() {
                r2 = Some(seq);
            } else {
                return Some(Err(anyhow::anyhow!(
                    "more than two primary records for template {:?}",
                    rec_name
                )));
            }
        }
    })
}

/// Resolved input mode; determined once at startup.
///
/// `AutoFasta`/`AutoFastq` cover single-input FASTX where interleaved layout
/// is auto-detected at iterator time (or skipped if `--per-record` is set).
/// `PairedFasta`/`PairedFastq` are the explicit two-file paths.
enum InputMode {
    AutoFasta,
    AutoFastq,
    PairedFasta,
    PairedFastq,
    /// SAM/BAM/CRAM: extension-resolved path; open by AlignmentFormat downstream.
    Sam,
    /// Sniffer fired and recognised text FASTX. Carries the reader at byte 0.
    SniffedFastx {
        kind: crate::FastxKind,
        gzipped: bool,
        reader: std::io::BufReader<std::fs::File>,
    },
    /// Sniffer fired and recognised an alignment format. Carries the reader
    /// at byte 0; the alignment writer constructs the right noodles reader.
    SniffedAlignment {
        format: crate::SniffedFormat,
        gzipped: bool,
        reader: std::io::BufReader<std::fs::File>,
    },
}

/// Common alignment-iterator dispatch: route to query-grouped or single-end
/// iterator based on `query_grouped` (from header) and `per_record` (CLI).
fn dispatch_alignment<I>(
    record_iter: I,
    query_grouped: bool,
    per_record: bool,
    output: impl Write,
) -> Result<()>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
{
    if query_grouped && !per_record {
        write_fasta(iter_query_grouped_alignment_templates(record_iter), output)
    } else {
        write_fasta(
            iter_single_end_alignment_templates(record_iter, per_record),
            output,
        )
    }
}

/// Write SAM/BAM/CRAM templates as FASTA, streaming without full-file buffering.
///
/// Reads the header, detects query-grouped ordering from `@HD`, and dispatches
/// to the query-grouped iterator unless `--per-record` is set, otherwise to
/// the single-end iterator. CRAM is handled in-place so that
/// `reader.records(&header)` never needs to escape its stack frame as a
/// boxed iterator.
fn write_sam_prep(args: &PrepArgs, output: impl Write) -> Result<()> {
    use crate::AlignmentFormat;
    use noodles::sam::alignment::record_buf::RecordBuf;

    match AlignmentFormat::from_path(&args.input) {
        AlignmentFormat::Bam => {
            let mut reader = crate::open_bam_reader(&args.input)?;
            let header = reader.read_header().with_context(|| {
                format!(
                    "failed to read SAM/BAM/CRAM header: {}",
                    args.input.display()
                )
            })?;
            let query_grouped = crate::is_query_grouped(&header);
            let mut buf = RecordBuf::default();
            let record_iter =
                std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                    Ok(0) => None,
                    Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                    Err(e) => Some(Err(e)),
                });
            dispatch_alignment(record_iter, query_grouped, args.per_record, output)
        }
        AlignmentFormat::Cram => {
            let mut reader = crate::open_cram_reader(&args.input, args.cram_reference.as_deref())?;
            let header = reader.read_header().with_context(|| {
                format!(
                    "failed to read SAM/BAM/CRAM header: {}",
                    args.input.display()
                )
            })?;
            crate::require_cram_reference_if_mapped(&header, args.cram_reference.as_deref())?;
            let query_grouped = crate::is_query_grouped(&header);
            // `reader.records(&header)` borrows both `reader` and `header` with tied
            // lifetimes. Using it in-place (no boxing) is fine because neither borrow
            // escapes this match arm's stack frame.
            dispatch_alignment(
                reader.records(&header),
                query_grouped,
                args.per_record,
                output,
            )
        }
        AlignmentFormat::Sam => {
            use noodles::sam;
            let file = std::fs::File::open(&args.input)
                .with_context(|| format!("failed to open SAM: {}", args.input.display()))?;
            let mut reader = sam::io::Reader::new(std::io::BufReader::new(file));
            let header = reader.read_header().with_context(|| {
                format!(
                    "failed to read SAM/BAM/CRAM header: {}",
                    args.input.display()
                )
            })?;
            let query_grouped = crate::is_query_grouped(&header);
            let mut buf = RecordBuf::default();
            let record_iter =
                std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                    Ok(0) => None,
                    Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                    Err(e) => Some(Err(e)),
                });
            dispatch_alignment(record_iter, query_grouped, args.per_record, output)
        }
    }
}

fn resolve_mode(args: &PrepArgs) -> Result<InputMode> {
    // Fast path: extension is unambiguous FASTX.
    match crate::infer_format(&args.input) {
        crate::InferredFormat::Fastx(crate::FastxKind::Fasta) => {
            return Ok(if args.input2.is_some() {
                InputMode::PairedFasta
            } else {
                InputMode::AutoFasta
            });
        }
        crate::InferredFormat::Fastx(crate::FastxKind::Fastq) => {
            return Ok(if args.input2.is_some() {
                InputMode::PairedFastq
            } else {
                InputMode::AutoFastq
            });
        }
        crate::InferredFormat::Alignment(_) => {
            // Fall through: extension may be a real .sam/.bam/.cram (use
            // existing alignment dispatch) or missing/unknown (sniff).
        }
    }

    let ext = args
        .input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase());
    let unambiguous_alignment = matches!(ext.as_deref(), Some("sam") | Some("bam") | Some("cram"));
    if unambiguous_alignment {
        if args.input2.is_some() {
            anyhow::bail!(
                "-2/--input-2 is not valid for SAM/BAM/CRAM input; \
                 SAM/BAM/CRAM pairing is detected from the @HD header automatically"
            );
        }
        return Ok(InputMode::Sam);
    }

    // Sniff. The reader returned holds the bytes already peeked, so any
    // downstream noodles reader picks up at byte 0.
    let (format, gzipped, reader) = crate::sniff_input(&args.input)
        .with_context(|| format!("failed to open input: {}", args.input.display()))?;
    match format {
        crate::SniffedFormat::Fasta => {
            if args.input2.is_some() {
                anyhow::bail!("-2/--input-2 with FASTA stdin is not supported (need two files)");
            }
            Ok(InputMode::SniffedFastx {
                kind: crate::FastxKind::Fasta,
                gzipped,
                reader,
            })
        }
        crate::SniffedFormat::Fastq => {
            if args.input2.is_some() {
                anyhow::bail!("-2/--input-2 with FASTQ stdin is not supported (need two files)");
            }
            Ok(InputMode::SniffedFastx {
                kind: crate::FastxKind::Fastq,
                gzipped,
                reader,
            })
        }
        crate::SniffedFormat::Sam | crate::SniffedFormat::Bam | crate::SniffedFormat::Cram => {
            if args.input2.is_some() {
                anyhow::bail!(
                    "-2/--input-2 is not valid for SAM/BAM/CRAM input; \
                     SAM/BAM/CRAM pairing is detected from the @HD header automatically"
                );
            }
            Ok(InputMode::SniffedAlignment {
                format,
                gzipped,
                reader,
            })
        }
        crate::SniffedFormat::Unknown => anyhow::bail!(
            "could not infer format from input head bytes for {}; \
             supply a file with a known extension",
            args.input.display()
        ),
    }
}

/// `true` when the resolved input is a CRAM stream (extension-detected or
/// sniffed). Used to gate `--cram-reference`.
fn mode_uses_cram(mode: &InputMode, input: &std::path::Path) -> bool {
    match mode {
        InputMode::Sam => matches!(
            crate::AlignmentFormat::from_path(input),
            crate::AlignmentFormat::Cram
        ),
        InputMode::SniffedAlignment { format, .. } => {
            matches!(format, crate::SniffedFormat::Cram)
        }
        _ => false,
    }
}

/// Run the `prep` command.
pub fn run_prep(args: PrepArgs) -> Result<()> {
    let mode = resolve_mode(&args)?;

    if args.cram_reference.is_some() && !mode_uses_cram(&mode, &args.input) {
        anyhow::bail!(
            "--cram-reference is only valid for CRAM input; \
             remove the flag or pass a CRAM file"
        );
    }

    let output_file = std::fs::File::create(&args.output)
        .with_context(|| format!("failed to create output: {}", args.output.display()))?;
    let output = std::io::BufWriter::new(output_file);

    match mode {
        InputMode::AutoFastq => write_fasta(iter_fastq_auto(&args.input, args.per_record)?, output),
        InputMode::AutoFasta => write_fasta(iter_fasta_auto(&args.input, args.per_record)?, output),
        InputMode::PairedFastq => write_fasta(
            iter_paired_fastq(
                &args.input,
                args.input2.as_deref().expect("input2 set for PairedFastq"),
            )?,
            output,
        ),
        InputMode::PairedFasta => write_fasta(
            iter_paired_fasta(
                &args.input,
                args.input2.as_deref().expect("input2 set for PairedFasta"),
            )?,
            output,
        ),
        InputMode::Sam => write_sam_prep(&args, output),
        InputMode::SniffedFastx {
            kind: crate::FastxKind::Fastq,
            gzipped,
            reader,
        } => {
            let r = crate::into_text_bufread(reader, gzipped);
            write_fasta(iter_fastq_auto_from_reader(r, args.per_record), output)
        }
        InputMode::SniffedFastx {
            kind: crate::FastxKind::Fasta,
            gzipped,
            reader,
        } => {
            let r = crate::into_text_bufread(reader, gzipped);
            write_fasta(iter_fasta_auto_from_reader(r, args.per_record), output)
        }
        InputMode::SniffedAlignment {
            format,
            gzipped,
            reader,
        } => write_sniffed_alignment_prep(&args, format, gzipped, reader, output),
    }
}

/// Run the alignment-templates writer for SAM/BAM/CRAM that arrived via
/// `sniff_input`. The reader is already positioned at byte 0 of the raw
/// stream.
fn write_sniffed_alignment_prep(
    args: &PrepArgs,
    format: crate::SniffedFormat,
    gzipped: bool,
    reader: std::io::BufReader<std::fs::File>,
    output: impl Write,
) -> Result<()> {
    use noodles::sam::alignment::record_buf::RecordBuf;

    let header_ctx = || {
        format!(
            "failed to read SAM/BAM/CRAM header: {}",
            args.input.display()
        )
    };

    match format {
        crate::SniffedFormat::Sam => {
            use noodles::sam;
            let r = crate::into_text_bufread(reader, gzipped);
            let mut reader = sam::io::Reader::new(r);
            let header = reader.read_header().with_context(header_ctx)?;
            let query_grouped = crate::is_query_grouped(&header);
            let mut buf = RecordBuf::default();
            let record_iter =
                std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                    Ok(0) => None,
                    Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                    Err(e) => Some(Err(e)),
                });
            dispatch_alignment(record_iter, query_grouped, args.per_record, output)
        }
        crate::SniffedFormat::Bam => {
            use noodles::bam;
            use noodles::bgzf;
            // Pseudo-paths (e.g. /dev/stdin) keep the buffered sniff reader so
            // the bytes already peeked are not lost; real paths reopen via the
            // multithreaded BGZF reader for throughput.
            if crate::is_pseudo_path(&args.input) {
                let mut reader = bam::io::Reader::from(bgzf::io::Reader::new(reader));
                let header = reader.read_header().with_context(header_ctx)?;
                let query_grouped = crate::is_query_grouped(&header);
                let mut buf = RecordBuf::default();
                let record_iter =
                    std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                        Ok(0) => None,
                        Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                        Err(e) => Some(Err(e)),
                    });
                dispatch_alignment(record_iter, query_grouped, args.per_record, output)
            } else {
                drop(reader);
                let mut reader = crate::open_bam_reader(&args.input)?;
                let header = reader.read_header().with_context(header_ctx)?;
                let query_grouped = crate::is_query_grouped(&header);
                let mut buf = RecordBuf::default();
                let record_iter =
                    std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                        Ok(0) => None,
                        Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                        Err(e) => Some(Err(e)),
                    });
                dispatch_alignment(record_iter, query_grouped, args.per_record, output)
            }
        }
        crate::SniffedFormat::Cram => {
            use noodles::cram;
            if crate::is_pseudo_path(&args.input) {
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(crate::build_fasta_repo(
                        args.cram_reference.as_deref(),
                    )?)
                    .build_from_reader(reader);
                let header = reader.read_header().with_context(header_ctx)?;
                crate::require_cram_reference_if_mapped(&header, args.cram_reference.as_deref())?;
                let query_grouped = crate::is_query_grouped(&header);
                dispatch_alignment(
                    reader.records(&header),
                    query_grouped,
                    args.per_record,
                    output,
                )
            } else {
                drop(reader);
                let mut reader =
                    crate::open_cram_reader(&args.input, args.cram_reference.as_deref())?;
                let header = reader.read_header().with_context(header_ctx)?;
                crate::require_cram_reference_if_mapped(&header, args.cram_reference.as_deref())?;
                let query_grouped = crate::is_query_grouped(&header);
                dispatch_alignment(
                    reader.records(&header),
                    query_grouped,
                    args.per_record,
                    output,
                )
            }
        }
        _ => unreachable!("non-alignment SniffedFormat reached alignment writer"),
    }
}

/// Write templates as FASTA.
///
/// Single-end: `>name\nseq\n`. Paired: `>name\nR1NR2\n`. Sequences are
/// uppercased before write so Kraken (v1 and v2) sees only `A/C/G/T/N`;
/// soft-masked lowercase bases would otherwise be treated as ambiguous and
/// break minimizer chains.
pub(crate) fn write_fasta(
    templates: impl Iterator<Item = anyhow::Result<Template>>,
    mut out: impl Write,
) -> Result<()> {
    // Flush unconditionally, even when a template errors mid-stream, so a
    // partial output is finalized and a flush error is surfaced rather than
    // swallowed by `BufWriter`'s drop.
    let body: Result<()> = (|| {
        for result in templates {
            let mut tmpl = result?;
            writeln!(out, ">{}", tmpl.name).context("failed to write FASTA header")?;
            tmpl.r1.make_ascii_uppercase();
            out.write_all(&tmpl.r1)
                .context("failed to write sequence")?;
            if let Some(ref mut r2) = tmpl.r2 {
                r2.make_ascii_uppercase();
                out.write_all(b"N").context("failed to write N separator")?;
                out.write_all(r2).context("failed to write mate sequence")?;
            }
            writeln!(out).context("failed to write record newline")?;
        }
        Ok(())
    })();
    crate::finish_after(body, || out.flush().context("failed to flush output"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fastq_bytes(records: &[(&str, &str, &str)]) -> Vec<u8> {
        let mut buf = Vec::new();
        for (name, seq, qual) in records {
            writeln!(buf, "@{name}").unwrap();
            writeln!(buf, "{seq}").unwrap();
            writeln!(buf, "+").unwrap();
            writeln!(buf, "{qual}").unwrap();
        }
        buf
    }

    fn fasta_bytes(records: &[(&str, &str)]) -> Vec<u8> {
        let mut buf = Vec::new();
        for (name, seq) in records {
            writeln!(buf, ">{name}").unwrap();
            writeln!(buf, "{seq}").unwrap();
        }
        buf
    }

    #[test]
    fn test_is_query_grouped_empty_header_is_false() {
        let header = noodles::sam::Header::default();
        assert!(!crate::is_query_grouped(&header));
    }

    #[test]
    fn test_is_query_grouped_so_queryname_is_true() {
        use noodles::sam;

        // Parse a SAM header with SO:queryname
        let header: sam::Header = "@HD\tVN:1.6\tSO:queryname\n"
            .parse()
            .expect("failed to parse SAM header");
        assert!(
            crate::is_query_grouped(&header),
            "SO:queryname should be query-grouped"
        );
    }

    #[test]
    fn test_is_query_grouped_go_query_is_true() {
        use noodles::sam;

        // Parse a SAM header with GO:query
        let header: sam::Header = "@HD\tVN:1.6\tGO:query\n"
            .parse()
            .expect("failed to parse SAM header");
        assert!(
            crate::is_query_grouped(&header),
            "GO:query should be query-grouped"
        );
    }

    #[test]
    fn test_write_fasta_single_end() {
        let templates = vec![Ok(Template {
            name: "read1".into(),
            r1: b"ACGT".to_vec(),
            r2: None,
        })];
        let mut out = Vec::new();
        write_fasta(templates.into_iter(), &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">read1\nACGT\n");
    }

    #[test]
    fn test_write_fasta_paired() {
        let templates = vec![Ok(Template {
            name: "pair1".into(),
            r1: b"AAAA".to_vec(),
            r2: Some(b"TTTT".to_vec()),
        })];
        let mut out = Vec::new();
        write_fasta(templates.into_iter(), &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">pair1\nAAAANTTTT\n");
    }

    #[test]
    fn test_write_fasta_multiple() {
        let templates = vec![
            Ok(Template {
                name: "r1".into(),
                r1: b"ACGT".to_vec(),
                r2: None,
            }),
            Ok(Template {
                name: "r2".into(),
                r1: b"GGCC".to_vec(),
                r2: None,
            }),
        ];
        let mut out = Vec::new();
        write_fasta(templates.into_iter(), &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains(">r1\nACGT\n"), "got: {s}");
        assert!(s.contains(">r2\nGGCC\n"), "got: {s}");
    }

    #[test]
    fn test_collect_paired_fastq() {
        let r1 = fastq_bytes(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "IIII")]);
        let r2 = fastq_bytes(&[("p1", "TTTT", "JJJJ"), ("p2", "GGGG", "JJJJ")]);
        let templates = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].r1, b"AAAA");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
    }

    #[test]
    fn test_collect_paired_fastq_name_mismatch_is_error() {
        let r1 = fastq_bytes(&[("p1", "AAAA", "IIII")]);
        let r2 = fastq_bytes(&[("p2", "TTTT", "JJJJ")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_paired_fastq_unequal_lengths_is_error() {
        let r1 = fastq_bytes(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "IIII")]);
        let r2 = fastq_bytes(&[("p1", "TTTT", "JJJJ")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_paired_fastq_empty_name_is_error() {
        let r1 = b"@\nACGT\n+\nIIII\n";
        let r2 = fastq_bytes(&[("p1", "TTTT", "JJJJ")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1.as_ref())),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err(), "expected error for empty FASTQ name");
        assert!(result.unwrap_err().to_string().contains("empty name"));
    }

    #[test]
    fn test_collect_paired_fasta() {
        let r1 = fasta_bytes(&[("s1", "AAAA"), ("s2", "CCCC")]);
        let r2 = fasta_bytes(&[("s1", "TTTT"), ("s2", "GGGG")]);
        let templates = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[1].r2.as_deref(), Some(b"GGGG".as_ref()));
    }

    #[test]
    fn test_collect_paired_fasta_name_mismatch_is_error() {
        let r1 = fasta_bytes(&[("s1", "AAAA")]);
        let r2 = fasta_bytes(&[("s2", "TTTT")]);
        let result = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_paired_fasta_unequal_lengths_is_error() {
        let r1 = fasta_bytes(&[("s1", "AAAA"), ("s2", "CCCC")]);
        let r2 = fasta_bytes(&[("s1", "TTTT")]);
        let result = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_paired_fasta_empty_name_is_error() {
        let r1 = b">\nACGT\n";
        let r2 = fasta_bytes(&[("s1", "TTTT")]);
        let result = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1.as_ref())),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err(), "expected error for empty FASTA name");
        assert!(result.unwrap_err().to_string().contains("empty name"));
    }

    #[test]
    fn test_iter_fastq_auto_single_end_distinct_names() {
        let data = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "TTTT", "JJJJ")]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "read1");
        assert_eq!(templates[0].r1, b"ACGT");
        assert!(templates[0].r2.is_none());
        assert_eq!(templates[1].name, "read2");
    }

    #[test]
    fn test_iter_fastq_auto_slash1_detected_as_interleaved() {
        let data = fastq_bytes(&[
            ("pair1/1", "AAAA", "IIII"),
            ("pair1/2", "TTTT", "IIII"),
            ("pair2/1", "CCCC", "JJJJ"),
            ("pair2/2", "GGGG", "JJJJ"),
        ]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "pair1");
        assert_eq!(templates[0].r1, b"AAAA");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
        assert_eq!(templates[1].name, "pair2");
        assert_eq!(templates[1].r2.as_deref(), Some(b"GGGG".as_ref()));
    }

    #[test]
    fn test_iter_fastq_auto_slash2_first_is_error() {
        let data = fastq_bytes(&[("pair1/2", "TTTT", "IIII"), ("pair1/1", "AAAA", "IIII")]);
        let result =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("/2"), "got: {msg}");
        assert!(msg.contains("--per-record"), "got: {msg}");
    }

    #[test]
    fn test_iter_fastq_auto_matching_names_detected_as_interleaved() {
        // Casava 1.8+: identical first-token name, pair info in description.
        let data = fastq_bytes(&[
            ("HISEQ:1:1234 1:N:0:ATC", "AAAA", "IIII"),
            ("HISEQ:1:1234 2:N:0:ATC", "TTTT", "IIII"),
            ("HISEQ:1:5678 1:N:0:ATC", "CCCC", "JJJJ"),
            ("HISEQ:1:5678 2:N:0:ATC", "GGGG", "JJJJ"),
        ]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "HISEQ:1:1234");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
        assert_eq!(templates[1].name, "HISEQ:1:5678");
    }

    #[test]
    fn test_iter_fastq_auto_per_record_preserves_slash1_slash2() {
        // With per_record=true, /1 /2 records are emitted as single-end with suffixes preserved.
        let data = fastq_bytes(&[("pair1/1", "AAAA", "IIII"), ("pair1/2", "TTTT", "IIII")]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), true)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "pair1/1");
        assert!(templates[0].r2.is_none());
        assert_eq!(templates[1].name, "pair1/2");
    }

    #[test]
    fn test_iter_fastq_auto_per_record_does_not_error_on_slash2_first() {
        // /2-first error must not fire when per_record=true.
        let data = fastq_bytes(&[("pair1/2", "TTTT", "IIII")]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), true)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "pair1/2");
    }

    #[test]
    fn test_iter_fastq_auto_empty_input() {
        let templates = iter_fastq_auto_from_reader(
            std::io::BufReader::new(std::io::Cursor::new(b"".as_ref())),
            false,
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
        assert!(templates.is_empty());
    }

    #[test]
    fn test_iter_fastq_auto_single_record_is_single_end() {
        let data = fastq_bytes(&[("only", "ACGT", "IIII")]);
        let templates =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "only");
        assert!(templates[0].r2.is_none());
    }

    #[test]
    fn test_iter_fastq_auto_interleaved_odd_count_is_error() {
        let data = fastq_bytes(&[
            ("pair1/1", "AAAA", "IIII"),
            ("pair1/2", "TTTT", "IIII"),
            ("pair2/1", "CCCC", "JJJJ"),
        ]);
        let result =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_fastq_auto_interleaved_name_mismatch_is_error() {
        let data = fastq_bytes(&[("pair1/1", "AAAA", "IIII"), ("other/2", "TTTT", "IIII")]);
        let result =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("--per-record"), "got: {msg}");
    }

    #[test]
    fn test_iter_fasta_auto_single_end_distinct_names() {
        let data = fasta_bytes(&[("seq1", "ACGT"), ("seq2", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "seq1");
        assert!(templates[0].r2.is_none());
    }

    #[test]
    fn test_iter_fasta_auto_slash1_detected_as_interleaved() {
        let data = fasta_bytes(&[("p1/1", "AAAA"), ("p1/2", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "p1");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
    }

    #[test]
    fn test_iter_fasta_auto_matching_names_detected_as_interleaved() {
        let data = fasta_bytes(&[("read1 1:N", "AAAA"), ("read1 2:N", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "read1");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
    }

    #[test]
    fn test_iter_fasta_auto_slash2_first_is_error() {
        let data = fasta_bytes(&[("p1/2", "TTTT"), ("p1/1", "AAAA")]);
        let result =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("/2"), "got: {msg}");
        assert!(msg.contains("--per-record"), "got: {msg}");
    }

    #[test]
    fn test_iter_fasta_auto_per_record_short_circuits() {
        let data = fasta_bytes(&[("p1/1", "AAAA"), ("p1/2", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), true)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "p1/1");
        assert!(templates[0].r2.is_none());
    }

    #[test]
    fn test_iter_fasta_auto_empty_input() {
        let templates = iter_fasta_auto_from_reader(
            std::io::BufReader::new(std::io::Cursor::new(b"".as_ref())),
            false,
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
        assert!(templates.is_empty());
    }

    #[test]
    fn test_iter_fasta_auto_single_record_is_single_end() {
        let data = fasta_bytes(&[("only", "ACGT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "only");
    }

    #[test]
    fn test_iter_fasta_auto_per_record_does_not_error_on_slash2_first() {
        // /2-first error must not fire when per_record=true.
        let data = fasta_bytes(&[("p1/2", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), true)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "p1/2");
    }

    #[test]
    fn test_iter_fasta_auto_interleaved_odd_count_is_error() {
        let data = fasta_bytes(&[("p1/1", "AAAA"), ("p1/2", "TTTT"), ("p2/1", "CCCC")]);
        let result =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_fasta_auto_first_record_empty_name_errors() {
        // A FASTA record with an empty header (`>` followed by whitespace
        // only) must produce a clear error rather than a silently empty name.
        let data = b">\nACGT\n";
        let result =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        let err = result.unwrap_err();
        assert!(format!("{err:#}").contains("empty name"));
    }

    #[test]
    fn test_iter_fasta_auto_per_record_passes_through_pair_suffixes() {
        // With per_record=true, pair-suffix detection is bypassed: every
        // record (including those ending in /1 or /2) becomes a single-end
        // template with its raw name.
        let data = fasta_bytes(&[("read/1", "AAAA"), ("read/2", "TTTT")]);
        let templates =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), true)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "read/1");
        assert_eq!(templates[1].name, "read/2");
    }

    #[test]
    fn test_parse_fastq_name_non_utf8_errors() {
        // parse_fastq_name must surface a non-UTF-8-name error rather than
        // panic; the message includes a lossy rendering of the bad bytes.
        let err = parse_fastq_name(&[0xFF, 0xFE]).unwrap_err();
        assert!(format!("{err:#}").contains("non-UTF-8"));
    }

    #[test]
    fn test_iter_fasta_auto_interleaved_name_mismatch_is_error() {
        let data = fasta_bytes(&[("p1/1", "AAAA"), ("other/2", "TTTT")]);
        let result =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("--per-record"), "got: {msg}");
    }

    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record_buf::{RecordBuf, Sequence};

    fn make_sam_record(name: &str, flags: Flags, seq: &[u8]) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(name.as_bytes().into());
        *r.flags_mut() = flags;
        *r.sequence_mut() = Sequence::from(seq.to_vec());
        r
    }

    #[test]
    fn test_collect_single_end_sam_basic() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record("r1", Flags::default(), b"ACGT")),
            Ok(make_sam_record("r2", Flags::default(), b"TTTT")),
        ];
        let templates = iter_single_end_alignment_templates(records.into_iter(), false)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "r1");
        assert_eq!(templates[0].r1, b"ACGT");
    }

    #[test]
    fn test_collect_single_end_sam_skips_secondary() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record("r1", Flags::SECONDARY, b"ACGT")),
            Ok(make_sam_record("r2", Flags::default(), b"TTTT")),
        ];
        let templates = iter_single_end_alignment_templates(records.into_iter(), false)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "r2");
    }

    #[test]
    fn test_collect_single_end_sam_per_record_appends_pair_suffix() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "pair",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "pair",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"TTTT",
            )),
        ];
        let templates = iter_single_end_alignment_templates(records.into_iter(), true)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "pair/1");
        assert_eq!(templates[1].name, "pair/2");
    }

    #[test]
    fn test_collect_single_end_sam_per_record_does_not_double_suffix() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "pair/1",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "pair/2",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"TTTT",
            )),
        ];
        let templates = iter_single_end_alignment_templates(records.into_iter(), true)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "pair/1");
        assert_eq!(templates[1].name, "pair/2");
    }

    #[test]
    fn test_collect_single_end_sam_errors_on_paired_flag() {
        let records: Vec<std::io::Result<RecordBuf>> =
            vec![Ok(make_sam_record("r1", Flags::SEGMENTED, b"ACGT"))];
        let result = iter_single_end_alignment_templates(records.into_iter(), false)
            .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("r1"));
    }

    #[test]
    fn test_collect_single_end_sam_per_record_does_not_error_on_paired_flag() {
        use noodles::sam;
        let sam_text =
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\nread1\t1\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let mut reader =
            sam::io::Reader::new(std::io::BufReader::new(std::io::Cursor::new(sam_text)));
        let header = reader.read_header().unwrap();
        let mut buf = sam::alignment::record_buf::RecordBuf::default();
        let record_iter =
            std::iter::from_fn(move || match reader.read_record_buf(&header, &mut buf) {
                Ok(0) => None,
                Ok(_) => Some(Ok(std::mem::take(&mut buf))),
                Err(e) => Some(Err(e)),
            });
        let templates = iter_single_end_alignment_templates(record_iter, true)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "read1");
        assert!(templates[0].r2.is_none());
    }

    #[test]
    fn test_collect_query_grouped_sam_pair() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"TTTT",
            )),
        ];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].r1, b"AAAA");
        assert_eq!(templates[0].r2.as_deref(), Some(b"TTTT".as_ref()));
    }

    #[test]
    fn test_collect_query_grouped_sam_singleton() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![Ok(make_sam_record(
            "r1",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            b"ACGT",
        ))];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert!(templates[0].r2.is_none());
    }

    #[test]
    fn test_collect_query_grouped_sam_skips_supplementary_only_group() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record("r1", Flags::SUPPLEMENTARY, b"ACGT")),
            Ok(make_sam_record(
                "r2",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"TTTT",
            )),
        ];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].name, "r2");
    }

    #[test]
    fn test_collect_query_grouped_sam_r2_without_r1_is_error() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![Ok(make_sam_record(
            "r1",
            Flags::SEGMENTED | Flags::LAST_SEGMENT,
            b"ACGT",
        ))];
        let result = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_query_grouped_sam_multiple_pairs() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "p1",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "p1",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"TTTT",
            )),
            Ok(make_sam_record(
                "p2",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"CCCC",
            )),
            Ok(make_sam_record(
                "p2",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"GGGG",
            )),
        ];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[1].r1, b"CCCC");
    }

    #[test]
    fn test_collect_query_grouped_sam_no_segment_flags_both_records() {
        // Two primary records for the same QNAME with no segment flags set.
        // First goes to r1 (start-record fallback), second must go to r2, not
        // be silently dropped.
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record("r1", Flags::default(), b"AAAA")),
            Ok(make_sam_record("r1", Flags::default(), b"TTTT")),
        ];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].r1, b"AAAA");
        assert_eq!(
            templates[0].r2.as_deref(),
            Some(b"TTTT".as_ref()),
            "second no-flag record should become r2, not be silently dropped"
        );
    }

    #[test]
    fn test_prep_cram_single_end() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();

        let in_cram = dir.path().join("input.cram");
        let header = sam::Header::default();
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("r1".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
            w.try_finish(&header).unwrap();
        }

        let out_fa = dir.path().join("output.fa");

        super::run_prep(super::PrepArgs {
            input: in_cram,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();

        let result = std::fs::read_to_string(&out_fa).unwrap();
        assert_eq!(result, ">r1\nACGT\n");
    }

    /// Build an indexed FASTA and a CRAM whose header carries `@SQ` entries for
    /// it, so reading the CRAM requires the reference to decode slices.
    /// Returns (cram_path, fasta_path).
    fn write_mapped_cram_and_fasta(
        dir: &std::path::Path,
    ) -> (std::path::PathBuf, std::path::PathBuf) {
        use noodles::sam;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let fa_path = dir.join("ref.fa");
        let fai_path = dir.join("ref.fa.fai");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();
        std::fs::write(&fai_path, b"chr1\t4\t6\t4\t5\n").unwrap();

        let cram_path = dir.join("in.cram");
        let mut header = sam::Header::default();
        let len = NonZeroUsize::new(4).unwrap();
        header
            .reference_sequences_mut()
            .insert(b"chr1".as_ref().into(), Map::<ReferenceSequence>::new(len));
        let mut w = crate::open_cram_writer(&cram_path, Some(&fa_path)).unwrap();
        w.write_header(&header).unwrap();
        w.try_finish(&header).unwrap();
        (cram_path, fa_path)
    }

    /// Regression: a mapped CRAM (with `@SQ` in its header) used to drive
    /// noodles' decoder past the header without a reference and panic with
    /// "invalid slice reference sequence name". The fix bails with a clear
    /// error before any record is read.
    #[test]
    fn test_run_prep_mapped_cram_without_reference_errors_cleanly() {
        let dir = tempfile::TempDir::new().unwrap();
        let (cram_path, _fa) = write_mapped_cram_and_fasta(dir.path());
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: cram_path,
            input2: None,
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--cram-reference"), "got: {msg}");
        assert!(msg.contains("reference sequences"), "got: {msg}");
    }

    /// Same regression as above but exercising the pseudo-path branch
    /// (`/dev/fd/N`); this is the path the user's `samtools view ... | krak prep`
    /// pipeline takes. Before the fix the pseudo-path branch in
    /// `write_sniffed_alignment_prep` skipped the reference check and panicked.
    #[cfg(unix)]
    #[test]
    fn test_run_prep_mapped_cram_pseudo_path_without_reference_errors_cleanly() {
        use std::os::fd::AsRawFd;
        let dir = tempfile::TempDir::new().unwrap();
        let (cram_path, _fa) = write_mapped_cram_and_fasta(dir.path());

        let f = std::fs::File::open(&cram_path).unwrap();
        let fd = f.as_raw_fd();
        let pseudo = std::path::PathBuf::from(format!("/dev/fd/{fd}"));

        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: pseudo,
            input2: None,
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--cram-reference"), "got: {msg}");
    }

    /// Sequences must be uppercased on write so soft-masked / lowercase bases
    /// from CRAM-decoded reference reconstitution don't become ambiguous in
    /// Kraken's k-mer / minimizer scanner.
    #[test]
    fn test_write_fasta_uppercases_sequences() {
        let templates = vec![
            Ok(Template {
                name: "single".to_string(),
                r1: b"acgtACGTNn".to_vec(),
                r2: None,
            }),
            Ok(Template {
                name: "paired".to_string(),
                r1: b"aaaa".to_vec(),
                r2: Some(b"tttt".to_vec()),
            }),
        ];
        let mut buf = Vec::new();
        super::write_fasta(templates.into_iter(), &mut buf).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert_eq!(out, ">single\nACGTACGTNN\n>paired\nAAAANTTTT\n");
    }

    #[test]
    fn test_run_prep_input2_with_sam_extension_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("a.sam");
        std::fs::write(&in_path, b"@HD\tVN:1.6\n").unwrap();
        let in2_path = dir.path().join("b.fq");
        std::fs::write(&in2_path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: Some(in2_path),
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("-2/--input-2 is not valid for SAM/BAM/CRAM"));
    }

    #[test]
    fn test_run_prep_input2_with_sniffed_fasta_stdin_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("piped"); // extensionless → sniff
        std::fs::write(&in_path, b">r1\nACGT\n").unwrap();
        let in2_path = dir.path().join("b.fa");
        std::fs::write(&in2_path, b">r2\nTTTT\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: Some(in2_path),
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("-2/--input-2 with FASTA stdin"));
    }

    #[test]
    fn test_run_prep_input2_with_sniffed_fastq_stdin_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("piped"); // extensionless → sniff
        std::fs::write(&in_path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let in2_path = dir.path().join("b.fq");
        std::fs::write(&in2_path, b"@r2\nTTTT\n+\nJJJJ\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: Some(in2_path),
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("-2/--input-2 with FASTQ stdin"));
    }

    #[test]
    fn test_run_prep_input2_with_sniffed_alignment_stdin_errors() {
        use noodles::bam;
        use noodles::sam;
        let dir = tempfile::TempDir::new().unwrap();
        // Build a real BAM, then rename to extensionless so dispatch sniffs it.
        let bam_path = dir.path().join("real.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&bam_path).unwrap();
            w.write_header(&sam::Header::default()).unwrap();
        }
        let in_path = dir.path().join("piped");
        std::fs::rename(&bam_path, &in_path).unwrap();
        let in2_path = dir.path().join("b.fq");
        std::fs::write(&in2_path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: Some(in2_path),
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("-2/--input-2 is not valid for SAM/BAM/CRAM"));
    }

    #[test]
    fn test_run_prep_unknown_format_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("piped"); // extensionless
                                                // Bytes that are neither FASTA, FASTQ, BAM, nor CRAM.
        std::fs::write(&in_path, b"this is not a sequencing file\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: None,
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("could not infer format"));
    }

    #[test]
    fn test_run_prep_cram_reference_with_non_cram_input_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let ref_path = dir.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: in_path,
            input2: None,
            per_record: false,
            output: out_fa,
            cram_reference: Some(ref_path),
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("--cram-reference is only valid for CRAM"));
    }

    #[test]
    fn test_run_prep_extensionless_path_with_gzipped_fastq() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("piped"); // no extension
        {
            let f = std::fs::File::create(&in_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            enc.finish().unwrap();
        }

        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: in_path,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();
        assert_eq!(std::fs::read_to_string(&out_fa).unwrap(), ">r1\nACGT\n");
    }

    #[test]
    fn test_run_prep_extensionless_path_with_bam() {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("real.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&bam_path).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("r1".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
            // Drop writer to flush BGZF.
        }
        let in_path = dir.path().join("piped"); // no extension
        std::fs::rename(&bam_path, &in_path).unwrap();

        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: in_path,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();
        assert_eq!(std::fs::read_to_string(&out_fa).unwrap(), ">r1\nACGT\n");
    }

    #[test]
    fn test_run_prep_sam_extension_unambiguous_path() {
        // .sam-extension fast-path through run_prep -> write_sam_prep::Sam.
        let dir = tempfile::TempDir::new().unwrap();
        let in_sam = dir.path().join("input.sam");
        std::fs::write(
            &in_sam,
            b"@HD\tVN:1.6\n\
              @SQ\tSN:chr1\tLN:8\n\
              read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: in_sam,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();
        assert!(std::fs::read_to_string(&out_fa)
            .unwrap()
            .contains(">read1\nACGT\n"));
    }

    #[test]
    fn test_run_prep_bam_extension_unambiguous_path() {
        // Direct .bam-extension dispatch through run_prep -> write_sam_prep::Bam.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("input.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("rB".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: in_bam,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();
        assert_eq!(std::fs::read_to_string(&out_fa).unwrap(), ">rB\nACGT\n");
    }

    #[test]
    fn test_run_prep_paired_fasta_inputs() {
        // Direct dispatch through run_prep with two FASTA inputs (-1/-2).
        // Exercises iter_paired_fasta + iter_paired_fasta_from_readers and
        // verifies the FASTA writer emits both halves with /1 and /2.
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("r1.fa");
        let r2 = dir.path().join("r2.fa");
        std::fs::write(&r1, b">pair1\nACGT\n>pair2\nAAAA\n").unwrap();
        std::fs::write(&r2, b">pair1\nTGCA\n>pair2\nTTTT\n").unwrap();

        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: r1,
            input2: Some(r2),
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();

        // Kraken paired-end input format: one record per template, with the
        // R1 + 'N' + R2 sequences concatenated under the shared base name.
        let body = std::fs::read_to_string(&out_fa).unwrap();
        assert_eq!(body, ">pair1\nACGTNTGCA\n>pair2\nAAAANTTTT\n");
    }

    #[test]
    fn test_run_prep_paired_fasta_unequal_lengths_errors() {
        // R1 has more records than R2 → unequal-length error.
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("r1.fa");
        let r2 = dir.path().join("r2.fa");
        std::fs::write(&r1, b">p1\nACGT\n>p2\nAAAA\n").unwrap();
        std::fs::write(&r2, b">p1\nTGCA\n").unwrap();

        let out_fa = dir.path().join("out.fa");
        let err = super::run_prep(super::PrepArgs {
            input: r1,
            input2: Some(r2),
            per_record: false,
            output: out_fa,
            cram_reference: None,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("unequal record counts"));
    }

    #[test]
    fn test_run_prep_extensionless_path_with_cram() {
        // Exercise the sniffed-CRAM real-on-disk path in
        // write_sniffed_alignment_prep: when the input lacks an extension
        // we sniff CRAM and reopen via open_cram_reader so the
        // --cram-reference repository (if any) is wired in correctly.
        use noodles::sam;
        let dir = tempfile::TempDir::new().unwrap();
        let cram_path = dir.path().join("real.cram");
        {
            let mut w = crate::open_cram_writer(&cram_path, None).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            w.try_finish(&header).unwrap();
        }
        let in_path = dir.path().join("piped"); // no extension
        std::fs::rename(&cram_path, &in_path).unwrap();

        let out_fa = dir.path().join("out.fa");
        super::run_prep(super::PrepArgs {
            input: in_path,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();
        // Empty CRAM yields empty FASTA; the test exercises the dispatch
        // path; record content is verified elsewhere.
        assert_eq!(std::fs::read_to_string(&out_fa).unwrap(), "");
    }

    #[test]
    fn test_prep_gzipped_fastq_input() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        let in_fq_gz = dir.path().join("input.fq.gz");
        {
            let f = std::fs::File::create(&in_fq_gz).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
                .unwrap();
            enc.finish().unwrap();
        }

        let out_fa = dir.path().join("output.fa");
        super::run_prep(super::PrepArgs {
            input: in_fq_gz,
            input2: None,
            per_record: false,
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();

        let got = std::fs::read_to_string(&out_fa).unwrap();
        assert_eq!(got, ">r1\nACGT\n>r2\nTTTT\n");
    }

    #[test]
    fn test_prep_cram_paired() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();

        // Header with SO:queryname so paired mode is auto-detected
        let header: sam::Header = "@HD\tVN:1.6\tSO:queryname\n"
            .parse()
            .expect("failed to parse SAM header");

        let in_cram = dir.path().join("input.cram");
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();

            // R1: FIRST_SEGMENT flag (UNMAPPED required so noodles CRAM stores raw sequence bytes)
            let mut r1 = RecordBuf::default();
            *r1.name_mut() = Some("pair1".as_bytes().into());
            *r1.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::UNMAPPED;
            *r1.sequence_mut() = Sequence::from(b"AAAA".to_vec());
            *r1.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r1).unwrap();

            // R2: LAST_SEGMENT flag (UNMAPPED required so noodles CRAM stores raw sequence bytes)
            let mut r2 = RecordBuf::default();
            *r2.name_mut() = Some("pair1".as_bytes().into());
            *r2.flags_mut() = Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::UNMAPPED;
            *r2.sequence_mut() = Sequence::from(b"TTTT".to_vec());
            *r2.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r2).unwrap();

            w.try_finish(&header).unwrap();
        }

        let out_fa = dir.path().join("output.fa");

        super::run_prep(super::PrepArgs {
            input: in_cram,
            input2: None,
            per_record: false, // auto-detected from SO:queryname header
            output: out_fa.clone(),
            cram_reference: None,
        })
        .unwrap();

        let result = std::fs::read_to_string(&out_fa).unwrap();
        assert_eq!(result, ">pair1\nAAAANTTTT\n");
    }

    #[test]
    fn test_check_pair_suffixes_accepts_none_none() {
        assert!(check_pair_suffixes("read1", "read1").is_ok());
    }

    #[test]
    fn test_check_pair_suffixes_accepts_slash1_slash2() {
        assert!(check_pair_suffixes("read1/1", "read1/2").is_ok());
    }

    #[test]
    fn test_check_pair_suffixes_rejects_slash1_slash1() {
        assert!(check_pair_suffixes("pair1/1", "pair1/1").is_err());
    }

    #[test]
    fn test_check_pair_suffixes_rejects_slash2_slash1() {
        assert!(check_pair_suffixes("pair1/2", "pair1/1").is_err());
    }

    #[test]
    fn test_check_pair_suffixes_rejects_slash1_then_none() {
        assert!(check_pair_suffixes("pair1/1", "pair1").is_err());
    }

    #[test]
    fn test_iter_fastq_auto_rejects_duplicate_slash1() {
        let data = fastq_bytes(&[("pair1/1", "AAAA", "IIII"), ("pair1/1", "TTTT", "IIII")]);
        let result =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_paired_fastq_rejects_slash2_slash1() {
        let r1 = fastq_bytes(&[("pair1/2", "TTTT", "IIII")]);
        let r2 = fastq_bytes(&[("pair1/1", "AAAA", "IIII")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_paired_fastq_rejects_slash1_then_none() {
        let r1 = fastq_bytes(&[("pair1/1", "AAAA", "IIII")]);
        let r2 = fastq_bytes(&[("pair1", "TTTT", "IIII")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_paired_fastq_r1_empty_after_strip_errors() {
        // R1 raw name is "/1" → strips to "" → R1-empty error.
        let r1 = fastq_bytes(&[("/1", "ACGT", "IIII")]);
        let r2 = fastq_bytes(&[("/2", "TTTT", "JJJJ")]);
        let result = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        let err = result.unwrap_err();
        assert!(format!("{err:#}").contains("empty name"));
    }

    #[test]
    fn test_iter_fastq_auto_interleaved_r1_empty_after_strip_errors() {
        // First record name "/1" detected as interleaved; second record "/2".
        // After strip both names are empty → R1-empty error fires.
        let data = fastq_bytes(&[("/1", "ACGT", "IIII"), ("/2", "TTTT", "JJJJ")]);
        let result =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>();
        let err = result.unwrap_err();
        assert!(format!("{err:#}").contains("empty name"));
    }

    #[test]
    fn test_iter_fasta_auto_interleaved_odd_records_errors() {
        // First record ends with /1 → interleaved branch; no R2 follows.
        let data = fasta_bytes(&[("p1/1", "AAAA")]);
        let err =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap_err();
        assert!(format!("{err:#}").contains("interleaved FASTA has an odd number of records"));
    }

    #[test]
    fn test_iter_fastq_auto_interleaved_odd_records_errors() {
        let data = fastq_bytes(&[("p1/1", "ACGT", "IIII")]);
        let err =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap_err();
        assert!(format!("{err:#}").contains("interleaved FASTQ has an odd number of records"));
    }

    #[test]
    fn test_iter_paired_fastq_r2_empty_after_strip_errors() {
        // R1 raw "X/1" → "X"; R2 raw "/2" → "" → R2-empty error.
        let r1 = fastq_bytes(&[("X/1", "ACGT", "IIII")]);
        let r2 = fastq_bytes(&[("/2", "TTTT", "JJJJ")]);
        let err = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap_err();
        assert!(format!("{err:#}").contains("R2 FASTQ record has empty name"));
    }

    #[test]
    fn test_iter_paired_fasta_r1_empty_after_strip_errors() {
        let r1 = fasta_bytes(&[("/1", "AAAA")]);
        let r2 = fasta_bytes(&[("/2", "TTTT")]);
        let err = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap_err();
        assert!(format!("{err:#}").contains("R1 FASTA record has empty name"));
    }

    #[test]
    fn test_iter_paired_fasta_r2_empty_after_strip_errors() {
        let r1 = fasta_bytes(&[("X/1", "AAAA")]);
        let r2 = fasta_bytes(&[("/2", "TTTT")]);
        let err = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap_err();
        assert!(format!("{err:#}").contains("R2 FASTA record has empty name"));
    }

    #[test]
    fn test_iter_fasta_auto_interleaved_r1_empty_after_strip_errors() {
        let data = fasta_bytes(&[("/1", "AAAA"), ("/2", "TTTT")]);
        let err =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap_err();
        assert!(format!("{err:#}").contains("interleaved FASTA R1 has empty name"));
    }

    #[test]
    fn test_iter_fasta_auto_interleaved_r2_empty_after_strip_errors() {
        let data = fasta_bytes(&[("X/1", "AAAA"), ("/2", "TTTT")]);
        let err =
            iter_fasta_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap_err();
        assert!(format!("{err:#}").contains("interleaved FASTA R2 has empty name"));
    }

    #[test]
    fn test_iter_fastq_auto_interleaved_r2_empty_after_strip_errors() {
        let data = fastq_bytes(&[("X/1", "ACGT", "IIII"), ("/2", "TTTT", "JJJJ")]);
        let err =
            iter_fastq_auto_from_reader(std::io::BufReader::new(std::io::Cursor::new(data)), false)
                .collect::<anyhow::Result<Vec<_>>>()
                .unwrap_err();
        assert!(format!("{err:#}").contains("interleaved FASTQ R2 has empty name"));
    }

    #[test]
    fn test_iter_paired_fastq_r2_with_no_kraken_suffix_alongside_r1_no_suffix_works() {
        // (None, None) is a valid combination; same first-token Casava names.
        let r1 = fastq_bytes(&[("read1", "ACGT", "IIII"), ("read2", "AAAA", "IIII")]);
        let r2 = fastq_bytes(&[("read1", "TGCA", "JJJJ"), ("read2", "TTTT", "JJJJ")]);
        let templates = iter_paired_fastq_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
        assert_eq!(templates.len(), 2);
        assert_eq!(templates[0].name, "read1");
        assert_eq!(templates[1].name, "read2");
    }

    #[test]
    fn test_iter_paired_fasta_rejects_duplicate_slash1() {
        let r1 = fasta_bytes(&[("pair1/1", "AAAA")]);
        let r2 = fasta_bytes(&[("pair1/1", "TTTT")]);
        let result = iter_paired_fasta_from_readers(
            std::io::BufReader::new(std::io::Cursor::new(r1)),
            std::io::BufReader::new(std::io::Cursor::new(r2)),
        )
        .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
    }

    #[test]
    fn test_iter_query_grouped_sam_two_primary_r1_records_is_error() {
        // Both records carry FIRST_SEGMENT; should error with "two primary R1".
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"CCCC",
            )),
        ];
        let err = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap_err();
        assert!(format!("{err:#}").contains("two primary R1 records"));
    }

    #[test]
    fn test_iter_query_grouped_sam_two_primary_r2_records_is_error() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"AAAA",
            )),
            Ok(make_sam_record(
                "r1",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"CCCC",
            )),
        ];
        let err = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap_err();
        assert!(format!("{err:#}").contains("two primary R2 records"));
    }

    #[test]
    fn test_iter_query_grouped_sam_three_primary_records_is_error() {
        let records: Vec<std::io::Result<RecordBuf>> = vec![
            Ok(make_sam_record("r1", Flags::default(), b"AAAA")),
            Ok(make_sam_record("r1", Flags::default(), b"TTTT")),
            Ok(make_sam_record("r1", Flags::default(), b"CCCC")),
        ];
        let result = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>();
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("more than two"), "got: {msg}");
    }

    fn rc(input: &[u8]) -> Vec<u8> {
        let mut v = input.to_vec();
        reverse_complement(&mut v);
        v
    }

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(rc(b"ACGT"), b"ACGT");
        assert_eq!(rc(b"AAAA"), b"TTTT");
        assert_eq!(rc(b"ACGTN"), b"NACGT");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        // Soft-masked bases (lowercase) must reverse-complement preserving case.
        assert_eq!(rc(b"acgt"), b"acgt");
        assert_eq!(rc(b"aaaa"), b"tttt");
        assert_eq!(rc(b"cccc"), b"gggg");
    }

    #[test]
    fn test_reverse_complement_iupac_passthrough() {
        // Bytes outside ACGTNacgtn pass through unchanged (mirror order
        // reversed). IUPAC ambiguity codes are preserved as-is rather than
        // silently corrupted.
        assert_eq!(rc(b"ARYn"), b"nYRT");
        assert_eq!(rc(b"X"), b"X");
        assert_eq!(rc(b""), b"");
    }

    #[test]
    fn test_pair_suffix_detect_variants() {
        // Direct coverage of PairSuffix::detect for each branch.
        assert_eq!(PairSuffix::detect("read"), PairSuffix::None);
        assert_eq!(PairSuffix::detect("read/1"), PairSuffix::Slash1);
        assert_eq!(PairSuffix::detect("read/2"), PairSuffix::Slash2);
        assert_eq!(PairSuffix::detect(""), PairSuffix::None);
    }

    #[test]
    fn test_check_pair_suffixes_rejects_slash2_then_slash2() {
        // (Slash2, Slash2); both /2; must error like the (Slash1, Slash1)
        // case: paired files duplicated or out-of-order.
        let err = check_pair_suffixes("read/2", "read/2").unwrap_err();
        assert!(format!("{err:#}").contains("mismatched pair suffixes"));
    }

    #[test]
    fn test_check_pair_suffixes_rejects_none_then_slash1() {
        let err = check_pair_suffixes("read", "read/1").unwrap_err();
        assert!(format!("{err:#}").contains("mismatched pair suffixes"));
    }

    #[test]
    fn test_iter_query_grouped_reverse_complemented_flag_applies() {
        // Both R1 and R2 carry REVERSE_COMPLEMENTED. The query-grouped
        // template iterator must RC each segment independently; both the
        // start record (line ~793) and the same-group continuation
        // record (line ~858).
        let r1 = make_sam_record(
            "pair1",
            Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::REVERSE_COMPLEMENTED,
            b"AAAA",
        );
        let r2 = make_sam_record(
            "pair1",
            Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED,
            b"GGGG",
        );
        let records: Vec<std::io::Result<RecordBuf>> = vec![Ok(r1), Ok(r2)];
        let templates = iter_query_grouped_alignment_templates(records.into_iter())
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(templates.len(), 1);
        assert_eq!(templates[0].r1, b"TTTT", "R1 must be RC'd");
        assert_eq!(
            templates[0].r2.as_deref(),
            Some(b"CCCC".as_ref()),
            "R2 must be RC'd"
        );
    }

    #[test]
    fn test_iter_single_end_sam_reverse_complemented_flag_applies() {
        // Without REVERSE_COMPLEMENTED: sequence emitted as-is.
        let plain: Vec<std::io::Result<RecordBuf>> =
            vec![Ok(make_sam_record("r1", Flags::default(), b"AAAA"))];
        let plain_templates = iter_single_end_alignment_templates(plain.into_iter(), false)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(plain_templates[0].r1, b"AAAA");

        // With REVERSE_COMPLEMENTED: reverse-complemented before emit.
        let rc: Vec<std::io::Result<RecordBuf>> = vec![Ok(make_sam_record(
            "r1",
            Flags::REVERSE_COMPLEMENTED,
            b"AAAA",
        ))];
        let rc_templates = iter_single_end_alignment_templates(rc.into_iter(), false)
            .collect::<anyhow::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(rc_templates[0].r1, b"TTTT");
    }

    #[test]
    fn test_record_name_rejects_empty() {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(b"".as_ref().into());
        let result = record_name(&r);
        assert!(result.is_err());
    }

    #[test]
    fn test_record_name_rejects_placeholder() {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(b"*".as_ref().into());
        let result = record_name(&r);
        assert!(result.is_err());
    }
}
