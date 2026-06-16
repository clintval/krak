//! Revert aligned N-calls in SAM/BAM/CRAM to reference bases (`n2ref`).
//!
//! Requires a FASTA index (`.fai`) alongside the reference.

use std::io::BufReader;
use std::path::Path;

use anyhow::{bail, Context, Result};
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use crate::AlignmentFormat;

/// Arguments for the `n2ref` command.
pub struct N2RefArgs {
    /// Input SAM/BAM/CRAM file.
    pub input: std::path::PathBuf,
    /// Output SAM/BAM/CRAM file.
    pub output: std::path::PathBuf,
    /// Reference FASTA file (must have a `.fai` index alongside it).
    pub reference: std::path::PathBuf,
    /// Replacement quality score for converted N-calls (`None` = keep original).
    pub qual: Option<u8>,
    /// Number of bgzf compression worker threads for BAM output. Default 1
    /// (one compressor + one writer thread pipelined with the n2ref loop).
    /// Ignored for SAM (no compression) and CRAM (per-block codecs).
    pub threads: usize,
    /// bgzf compression level (0-9) for BAM output. Default 5. Ignored for
    /// SAM (no compression) and CRAM (per-block codecs).
    pub compression_level: u32,
}

/// Type alias for the indexed FASTA reader used to fetch reference bases.
type RefReader = noodles::fasta::io::IndexedReader<noodles::fasta::io::BufReader<std::fs::File>>;

/// Run the `n2ref` command.
///
/// Dispatches by file extension when unambiguous (`.sam`/`.bam`/`.cram`), and
/// otherwise sniffs the input head bytes via [`crate::sniff_input`]. Sniffing
/// is what makes `/dev/stdin` work for any of the three formats; without it,
/// `/dev/stdin` falls back to BAM (`AlignmentFormat::from_path`).
pub fn run_n2ref(args: N2RefArgs) -> Result<()> {
    let mut ref_reader = open_ref_reader(&args.reference)?;

    let ext = args
        .input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase());
    let unambiguous = matches!(ext.as_deref(), Some("sam") | Some("bam") | Some("cram"));

    if unambiguous {
        return match AlignmentFormat::from_path(&args.input) {
            AlignmentFormat::Bam => {
                let reader = crate::open_bam_reader(&args.input)?;
                n2ref_bam_with_reader(&args, reader, &mut ref_reader)
            }
            AlignmentFormat::Sam => {
                let file = std::fs::File::open(&args.input)
                    .with_context(|| format!("failed to open SAM: {}", args.input.display()))?;
                let reader = noodles::sam::io::Reader::new(BufReader::new(file));
                n2ref_sam_with_reader(&args, reader, &mut ref_reader)
            }
            AlignmentFormat::Cram => {
                let reader = crate::open_cram_reader(&args.input, Some(&args.reference))?;
                n2ref_cram_with_reader(&args, reader, &mut ref_reader)
            }
        };
    }

    let (sniffed, gzipped, peek_reader) = crate::sniff_input(&args.input)
        .with_context(|| format!("failed to open input: {}", args.input.display()))?;

    match sniffed {
        crate::SniffedFormat::Bam => {
            use noodles::bam;
            use noodles::bgzf;
            // BAM is always BGZF-framed; the BGZF reader wants the raw stream,
            // so we ignore the `gzipped` flag and never pre-decompress.
            let reader = bam::io::Reader::from(bgzf::io::Reader::new(peek_reader));
            n2ref_bam_with_reader(&args, reader, &mut ref_reader)
        }
        crate::SniffedFormat::Sam => {
            use noodles::sam;
            let r = crate::into_text_bufread(peek_reader, gzipped);
            let reader = sam::io::Reader::new(r);
            n2ref_sam_with_reader(&args, reader, &mut ref_reader)
        }
        crate::SniffedFormat::Cram => {
            // CRAM cannot be gzipped (sniff_input rejects gzipped CRAM). The
            // CRAM reader requires a reference repository to decode mapped
            // reads; wire in the same indexed FASTA we use for N-call
            // reversion.
            use noodles::cram;
            use noodles::fasta;
            use noodles::fasta::repository::adapters::IndexedReader as FastaIndexedAdapter;
            let adapter = fasta::io::indexed_reader::Builder::default()
                .build_from_path(&args.reference)
                .map(FastaIndexedAdapter::new)
                .with_context(|| {
                    format!(
                        "failed to open indexed (.fai present) FASTA: {}",
                        args.reference.display()
                    )
                })?;
            let repo = fasta::Repository::new(adapter);
            let reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_reader(peek_reader);
            n2ref_cram_with_reader(&args, reader, &mut ref_reader)
        }
        crate::SniffedFormat::Fasta | crate::SniffedFormat::Fastq => {
            bail!(
                "n2ref requires SAM/BAM/CRAM input; got FASTA/FASTQ for {}",
                args.input.display()
            )
        }
        crate::SniffedFormat::Unknown => bail!(
            "could not infer format from input head bytes for {}; \
             supply a file with a known extension",
            args.input.display()
        ),
    }
}

fn n2ref_bam_with_reader<R: std::io::Read>(
    args: &N2RefArgs,
    mut reader: noodles::bam::io::Reader<R>,
    ref_reader: &mut RefReader,
) -> Result<()> {
    use gzp::ZWriter as _;
    use noodles::bam;

    let header = reader.read_header().context("failed to read BAM header")?;

    let file = std::fs::File::create(&args.output)
        .with_context(|| format!("failed to create output BAM: {}", args.output.display()))?;
    let parz = gzp::par::compress::ParCompressBuilder::<gzp::deflate::Bgzf>::new()
        .num_threads(args.threads.max(1))
        .context("invalid --threads value for BAM bgzf writer")?
        .compression_level(gzp::Compression::new(args.compression_level))
        .from_writer(file);
    let mut writer = bam::io::Writer::from(parz);
    writer
        .write_header(&header)
        .context("failed to write BAM header")?;

    // Finalize the BGZF stream unconditionally, even if the loop errored: a
    // threaded gzp writer dropped without an explicit finish panics in its Drop
    // and masks the real error. `finish_after` preserves the loop's error.
    let body = process_records(
        reader.record_bufs(&header),
        &header,
        &mut writer,
        ref_reader,
        args.qual,
    );
    crate::finish_after(body, || {
        writer
            .into_inner()
            .finish()
            .map(|_| ())
            .map_err(|e| anyhow::anyhow!("failed to finish BAM BGZF stream: {e}"))
    })?;
    Ok(())
}

fn n2ref_sam_with_reader<R: std::io::BufRead>(
    args: &N2RefArgs,
    mut reader: noodles::sam::io::Reader<R>,
    ref_reader: &mut RefReader,
) -> Result<()> {
    use noodles::sam;

    let header = reader.read_header().context("failed to read SAM header")?;

    let out_file = std::fs::File::create(&args.output)
        .with_context(|| format!("failed to create output SAM: {}", args.output.display()))?;
    let mut writer = sam::io::Writer::new(std::io::BufWriter::new(out_file));
    writer
        .write_header(&header)
        .context("failed to write SAM header")?;

    // Flush unconditionally, even if the loop errored, so a partial SAM is
    // flushed (and a flush error surfaces) rather than left in the unflushed
    // buffer. The loop's error takes precedence.
    let body = process_records(
        reader.record_bufs(&header),
        &header,
        &mut writer,
        ref_reader,
        args.qual,
    );
    crate::finish_after(body, || {
        use std::io::Write as _;
        writer
            .into_inner()
            .flush()
            .context("failed to flush SAM writer")
    })?;
    Ok(())
}

fn n2ref_cram_with_reader<R: std::io::Read>(
    args: &N2RefArgs,
    mut reader: noodles::cram::io::Reader<R>,
    ref_reader: &mut RefReader,
) -> Result<()> {
    let header = reader.read_header().context("failed to read CRAM header")?;

    let mut writer = crate::open_cram_writer(&args.output, Some(&args.reference))?;
    writer
        .write_header(&header)
        .context("failed to write CRAM header")?;

    // Finish the CRAM writer unconditionally, even if the loop errored, so a
    // partial CRAM still gets its EOF container (and a finish error surfaces).
    // The loop's error takes precedence.
    let body = process_records(
        reader.records(&header),
        &header,
        &mut writer,
        ref_reader,
        args.qual,
    );
    crate::finish_after(body, || {
        writer
            .try_finish(&header)
            .context("failed to finish CRAM writer")
    })?;
    Ok(())
}

/// Open an indexed FASTA reader for fast random access to the reference.
fn open_ref_reader(path: &Path) -> Result<RefReader> {
    info!("Opening indexed reference: {}", path.display());
    noodles::fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| {
            format!(
                "failed to open indexed FASTA (.fai present): {}",
                path.display()
            )
        })
}

/// Compute the number of reference bases consumed by a record's CIGAR.
fn reference_span(record: &noodles::sam::alignment::record_buf::RecordBuf) -> usize {
    use noodles::sam::alignment::record::cigar::op::Kind;
    record
        .cigar()
        .as_ref()
        .iter()
        .filter_map(|op| match op.kind() {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => Some(op.len()),
            _ => None,
        })
        .sum()
}

/// Core record-processing loop shared across SAM, BAM, and CRAM.
///
/// Iterates `records`, reverts any aligned N-calls to their reference bases via
/// `ref_reader`, and writes each (possibly modified) record to `writer`.
/// The caller is responsible for reading/writing the header and finishing the writer.
fn process_records<W>(
    records: impl Iterator<Item = std::io::Result<noodles::sam::alignment::RecordBuf>>,
    header: &noodles::sam::Header,
    writer: &mut W,
    ref_reader: &mut RefReader,
    qual: Option<u8>,
) -> Result<()>
where
    W: AlignmentWrite,
{
    let mut converted = 0u64;
    let mut total = 0u64;

    for result in records {
        let mut record = result.context("failed to read record")?;
        total += 1;

        if !record.flags().is_unmapped() {
            let ref_seq = fetch_ref_seq(ref_reader, header, &record)?;
            converted += revert_n_calls(&mut record, &ref_seq, qual)?;
        }

        writer
            .write_alignment_record(header, &record)
            .context("failed to write record")?;
    }

    info!("Converted {converted} N-calls across {total} records.");
    Ok(())
}

/// Fetch the reference bases covering a record's aligned span from an indexed FASTA.
///
/// Returns the bases as an uppercase `Vec<u8>` starting at the record's alignment start
/// position, so index 0 of the returned slice corresponds to the first aligned reference base.
fn fetch_ref_seq(
    ref_reader: &mut RefReader,
    header: &noodles::sam::Header,
    record: &noodles::sam::alignment::record_buf::RecordBuf,
) -> Result<Vec<u8>> {
    use noodles::core::{Position, Region};

    let ref_name = resolve_ref_name(header, record)?;
    let Some(aln_start_pos) = record.alignment_start() else {
        bail!("mapped record has no alignment start");
    };
    let aln_start_1 = usize::from(aln_start_pos);
    let span = reference_span(record);
    if span == 0 {
        return Ok(Vec::new());
    }
    let aln_end_1 = aln_start_1 + span - 1;

    let start = Position::new(aln_start_1)
        .with_context(|| format!("invalid alignment start position: {aln_start_1}"))?;
    let end = Position::new(aln_end_1)
        .with_context(|| format!("invalid alignment end position: {aln_end_1}"))?;
    let region = Region::new(ref_name.as_bytes(), start..=end);

    let ref_record = ref_reader.query(&region).with_context(|| {
        format!("failed to query reference region {ref_name}:{aln_start_1}-{aln_end_1}")
    })?;

    Ok(ref_record
        .sequence()
        .as_ref()
        .iter()
        .map(|b: &u8| b.to_ascii_uppercase())
        .collect())
}

/// Get the reference sequence name for a record from the SAM header.
fn resolve_ref_name(
    header: &noodles::sam::Header,
    record: &noodles::sam::alignment::record_buf::RecordBuf,
) -> Result<String> {
    let Some(ref_id) = record.reference_sequence_id() else {
        bail!("mapped record has no reference sequence ID");
    };
    let (name, _) = header
        .reference_sequences()
        .get_index(ref_id)
        .ok_or_else(|| anyhow::anyhow!("reference sequence ID {ref_id} not in header"))?;
    Ok(std::str::from_utf8(name)
        .context("reference sequence name is not valid UTF-8")?
        .to_owned())
}

/// Revert N bases to their reference equivalents using the CIGAR alignment.
///
/// `ref_seq` must contain the reference bases starting at the record's alignment start
/// position (i.e. index 0 = first aligned reference base), as returned by [`fetch_ref_seq`].
///
/// Returns the number of bases converted.
fn revert_n_calls(
    record: &mut noodles::sam::alignment::record_buf::RecordBuf,
    ref_seq: &[u8],
    new_qual: Option<u8>,
) -> Result<u64> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    if ref_seq.is_empty() {
        return Ok(0);
    }
    // Fast path: secondaries with `SEQ=*` carry no read bases, so there's
    // nothing to revert and walking the CIGAR is wasted work.
    if record.sequence().as_ref().is_empty() {
        return Ok(0);
    }

    // Collect CIGAR ops to avoid borrow conflict while mutating the record.
    let cigar_ops: Vec<_> = record.cigar().as_ref().to_vec();
    let mut new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
    let mut new_quals: Vec<u8> = record.quality_scores().iter().collect();
    let mut read_pos = 0usize;
    let mut ref_pos = 0usize; // offset into ref_seq (which starts at alignment_start)
    let mut converted = 0u64;

    for op in &cigar_ops {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..op.len() {
                    if read_pos + i < new_seq.len()
                        && ref_pos + i < ref_seq.len()
                        && new_seq[read_pos + i].eq_ignore_ascii_case(&b'N')
                    {
                        new_seq[read_pos + i] = ref_seq[ref_pos + i];
                        if let Some(q) = new_qual {
                            if read_pos + i < new_quals.len() {
                                new_quals[read_pos + i] = q;
                            }
                        }
                        converted += 1;
                    }
                }
                read_pos += op.len();
                ref_pos += op.len();
            }
            Kind::Insertion | Kind::SoftClip => {
                read_pos += op.len();
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += op.len();
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if converted > 0 {
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        *record.sequence_mut() = Sequence::from(new_seq);
        if new_qual.is_some() {
            *record.quality_scores_mut() = QualityScores::from(new_quals);
        }
    }

    Ok(converted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::record::{
        cigar::{op::Kind, Op},
        Flags,
    };
    use noodles::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

    fn make_mapped_record(
        seq: &[u8],
        quals: &[u8],
        cigar_ops: Vec<Op>,
        ref_id: usize,
        aln_start: usize, // 1-based
    ) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.sequence_mut() = Sequence::from(seq.to_vec());
        *r.quality_scores_mut() = QualityScores::from(quals.to_vec());
        *r.cigar_mut() = Cigar::from(cigar_ops);
        *r.reference_sequence_id_mut() = Some(ref_id);
        *r.alignment_start_mut() = Position::new(aln_start);
        *r.flags_mut() = Flags::default();
        r
    }

    // In tests, ref_seq is the full chromosome starting at position 0.
    // All test records use aln_start=1 (1-based), so the region-local ref_seq
    // is identical to the full sequence, and ref_pos starts at 0 correctly.

    #[test]
    fn test_n_replaced_by_reference_base() {
        // seq: ANNT at positions 0-3, ref: ACGT
        let ref_seq = b"ACGTACGT";
        let seq = b"ANNT";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, ref_seq, None).unwrap();
        assert_eq!(converted, 2);
        let new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(&new_seq, b"ACGT");
        // Quality scores unchanged.
        let new_q: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(new_q, quals);
    }

    #[test]
    fn test_quality_replaced_when_qual_specified() {
        let ref_seq = b"AAAA";
        let seq = b"ANNA";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        revert_n_calls(&mut record, ref_seq, Some(5)).unwrap();
        let new_q: Vec<u8> = record.quality_scores().iter().collect();
        assert_eq!(new_q[0], 30); // unchanged (was A)
        assert_eq!(new_q[1], 5); // replaced
        assert_eq!(new_q[2], 5); // replaced
        assert_eq!(new_q[3], 30); // unchanged
    }

    #[test]
    fn test_no_n_no_conversion() {
        let ref_seq = b"ACGT";
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, ref_seq, None).unwrap();
        assert_eq!(converted, 0);
    }

    #[test]
    fn test_insertion_skips_ref_position() {
        // Read: AANNG (5 bases), CIGAR: 3M 2I
        // M positions: read[0..3] vs ref[0..3]="AAT"
        // I positions: read[3..5] not in ref -> Ns not replaced
        let ref_seq = b"AATXYZ";
        let seq = b"AANNG";
        let quals = vec![30u8; 5];
        let ops = vec![Op::new(Kind::Match, 3), Op::new(Kind::Insertion, 2)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, ref_seq, None).unwrap();
        assert_eq!(converted, 1); // read[2] N -> T
        let new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(new_seq[2], b'T');
        assert_eq!(new_seq[3], b'N'); // insertion, not replaced
    }

    #[test]
    fn test_soft_clip_skips_replacement() {
        // Read: NNACGT (6 bases), CIGAR: 2S4M
        // Soft-clipped Ns (positions 0-1) are not aligned -> not replaced.
        let ref_seq = b"ACGT";
        let seq = b"NNACGT";
        let quals = vec![30u8; 6];
        let ops = vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, ref_seq, None).unwrap();
        assert_eq!(converted, 0); // aligned ACGT has no Ns
        let new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(new_seq[0], b'N'); // soft-clipped, unchanged
        assert_eq!(new_seq[1], b'N'); // soft-clipped, unchanged
    }

    #[test]
    fn test_deletion_advances_ref_position() {
        // Read: ACT (3 bases), CIGAR: 1M 2D 2M
        // ref: A__GT (del skips ref[1..3], then maps GT)
        // seq[0]=A, ref[0]=A -> no change; seq[1..3] map to ref[3..5]
        let ref_seq = b"ANNGT";
        let seq = b"NNT";
        let quals = vec![30u8; 3];
        let ops = vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 2),
        ];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, ref_seq, None).unwrap();
        // seq[0] N -> ref[0] A; seq[1] N -> ref[3] G; seq[2] T not N
        assert_eq!(converted, 2);
        let new_seq: Vec<u8> = record.sequence().as_ref().to_vec();
        assert_eq!(new_seq[0], b'A');
        assert_eq!(new_seq[1], b'G');
        assert_eq!(new_seq[2], b'T');
    }

    #[test]
    fn test_empty_ref_seq_returns_zero() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, b"", None).unwrap();
        assert_eq!(converted, 0);
    }

    #[test]
    fn test_empty_seq_returns_zero_fast_path() {
        // Secondary records can carry SEQ=*; revert_n_calls must short-circuit
        // before walking the CIGAR, even with valid ref_seq present.
        use noodles::sam::alignment::record_buf::Sequence;
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(b"", &quals, ops, 0, 1);
        // Force empty sequence (constructor may have populated it).
        *record.sequence_mut() = Sequence::from(Vec::new());

        let converted = revert_n_calls(&mut record, b"ACGT", Some(40)).unwrap();
        assert_eq!(converted, 0, "empty SEQ must short-circuit");
        assert!(
            record.sequence().as_ref().is_empty(),
            "SEQ must remain empty"
        );
    }

    #[test]
    fn test_lowercase_n_is_replaced() {
        // The N detection must be ASCII-case-insensitive: lowercase 'n' is
        // a legal masked-base notation and should still be reverted.
        let seq = b"AnGT";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);

        let converted = revert_n_calls(&mut record, b"ACGT", None).unwrap();
        assert_eq!(converted, 1);
        let new_seq = record.sequence().as_ref();
        assert_eq!(new_seq[1], b'C', "lowercase n must be reverted to ref");
    }

    #[test]
    fn test_no_quality_with_explicit_qual_no_panic() {
        // When the input record has no quality scores (`*`) and --qual is
        // provided, conversion still succeeds (the qual is silently ignored
        // for the N positions); no panic on the empty new_quals indexing.
        use noodles::sam::alignment::record_buf::QualityScores;
        let seq = b"NCGT";
        let quals: Vec<u8> = vec![];
        let ops = vec![Op::new(Kind::Match, 4)];
        let mut record = make_mapped_record(seq, &quals, ops, 0, 1);
        // Force empty quality scores.
        *record.quality_scores_mut() = QualityScores::from(Vec::new());

        let converted = revert_n_calls(&mut record, b"ACGT", Some(40)).unwrap();
        assert_eq!(converted, 1);
        let new_seq = record.sequence().as_ref();
        assert_eq!(new_seq[0], b'A', "N must be replaced");
    }

    #[test]
    fn test_reference_span_match_only() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        let ops = vec![Op::new(Kind::Match, 4)];
        let record = make_mapped_record(seq, &quals, ops, 0, 1);
        assert_eq!(reference_span(&record), 4);
    }

    #[test]
    fn test_reference_span_with_indels() {
        let seq = b"ACGT";
        let quals = vec![30u8; 4];
        // 2M 1D 2M 1I 1M -> ref span = 2+1+2+0+1 = 6
        let ops = vec![
            Op::new(Kind::Match, 2),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, 2),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Match, 1),
        ];
        let record = make_mapped_record(seq, &quals, ops, 0, 1);
        assert_eq!(reference_span(&record), 6);
    }

    // Format detection is tested centrally in annotate::tests::test_format_detection.

    /// Integration test for CRAM processing.
    ///
    /// Uses `samtools` to generate a well-formed input CRAM and verifies that
    /// Write a minimal indexed FASTA at `dir/ref.fa` and return its path.
    fn write_ref_fa(dir: &std::path::Path) -> std::path::PathBuf {
        let fa_path = dir.join("ref.fa");
        let fai_path = dir.join("ref.fa.fai");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();
        // FAI: name  len  offset  bases_per_line  bytes_per_line
        std::fs::write(&fai_path, b"chr1\t4\t6\t4\t5\n").unwrap();
        fa_path
    }

    #[test]
    fn test_run_n2ref_bam_unambiguous_path() {
        // Direct dispatch through the .bam-extension branch in run_n2ref.
        // Verifies a single N-call is reverted to its reference base.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value;
        let _ = Value::Int32(0); // suppress unused-import warning if any

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());

        let in_path = dir.path().join("input.bam");
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_path).unwrap();
            w.write_header(&header).unwrap();
            let r = make_mapped_record(b"ANGT", &[30u8; 4], vec![Op::new(Kind::Match, 4)], 0, 1);
            w.write_alignment_record(&header, &r).unwrap();
        }

        let out_path = dir.path().join("out.bam");
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        // Read back and verify the N at position 1 was reverted to 'C'.
        let mut r = bam::io::reader::Builder.build_from_path(&out_path).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        let new_seq: Vec<u8> = recs[0].sequence().as_ref().to_vec();
        assert_eq!(&new_seq, b"ACGT");
    }

    /// Build a BAM with 32 mapped records (each containing an N to revert),
    /// then run `n2ref` writing to `out_name` with the given `threads` /
    /// `compression_level`. Returns (file size in bytes, count of records
    /// whose SEQ was successfully reverted to the reference `A` at position 1).
    fn run_n2ref_bam(
        tmpdir: &std::path::Path,
        out_name: &str,
        threads: usize,
        compression_level: u32,
    ) -> (u64, usize) {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;

        // 32-base reference of pure A's gives the compressor predictable
        // redundancy across records (so level 1 vs 9 differ measurably).
        let fa_path = tmpdir.join("ref.fa");
        let fai_path = tmpdir.join("ref.fa.fai");
        std::fs::write(&fa_path, b">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n").unwrap();
        std::fs::write(&fai_path, b"chr1\t32\t6\t32\t33\n").unwrap();

        let in_path = tmpdir.join("in.bam");
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:32\n".parse().unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_path).unwrap();
            w.write_header(&header).unwrap();
            for _ in 0..32 {
                // SEQ starts with N at pos 1, which n2ref must revert to 'A'.
                let r = make_mapped_record(
                    b"NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                    &[30u8; 32],
                    vec![Op::new(Kind::Match, 32)],
                    0,
                    1,
                );
                w.write_alignment_record(&header, &r).unwrap();
            }
        }

        let out_path = tmpdir.join(out_name);
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads,
            compression_level,
        })
        .unwrap();

        let size = std::fs::metadata(&out_path).unwrap().len();
        let mut r = bam::io::reader::Builder.build_from_path(&out_path).unwrap();
        let h = r.read_header().unwrap();
        let reverted = r
            .record_bufs(&h)
            .map(|res| res.unwrap())
            .filter(|rec| rec.sequence().as_ref().first().copied() == Some(b'A'))
            .count();
        (size, reverted)
    }

    #[test]
    fn test_n2ref_bam_higher_compression_level_yields_smaller_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let (size_low, reverted_low) = run_n2ref_bam(dir.path(), "low.bam", 1, 1);
        let (size_high, reverted_high) = run_n2ref_bam(dir.path(), "high.bam", 1, 9);
        assert_eq!(reverted_low, 32);
        assert_eq!(reverted_high, 32);
        assert!(
            size_high < size_low,
            "expected level 9 ({size_high} bytes) < level 1 ({size_low} bytes)"
        );
    }

    #[test]
    fn test_n2ref_bam_threads_one_and_many_round_trip_identically() {
        let dir = tempfile::TempDir::new().unwrap();
        let (_, reverted_serial) = run_n2ref_bam(dir.path(), "t1.bam", 1, 5);
        let (_, reverted_parallel) = run_n2ref_bam(dir.path(), "t4.bam", 4, 5);
        assert_eq!(reverted_serial, 32);
        assert_eq!(reverted_parallel, 32);
    }

    #[test]
    fn test_run_n2ref_sam_unambiguous_path() {
        // Direct dispatch through the .sam-extension branch in run_n2ref.
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());

        let in_path = dir.path().join("input.sam");
        std::fs::write(
            &in_path,
            b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNGT\tIIII\n",
        )
        .unwrap();

        let out_path = dir.path().join("out.sam");
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: Some(40),
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let body = std::fs::read_to_string(&out_path).unwrap();
        // Output SAM should contain ACGT (Ns reverted) and quality I (=40 phred+33)
        // at the converted positions.
        assert!(
            body.contains("\tACGT\t"),
            "expected reverted seq in: {body}"
        );
    }

    #[test]
    fn test_run_n2ref_unmapped_record_is_passthrough() {
        // Unmapped records have no alignment context; n2ref must skip them
        // without touching SEQ even if they contain Ns.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());

        let in_path = dir.path().join("input.bam");
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_path).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"unmapped".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"NNNN".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            *r.flags_mut() = Flags::UNMAPPED;
            w.write_alignment_record(&header, &r).unwrap();
        }

        let out_path = dir.path().join("out.bam");
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out_path).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(
            recs[0].sequence().as_ref().to_vec(),
            b"NNNN",
            "unmapped Ns must NOT be reverted"
        );
    }

    #[test]
    fn test_run_n2ref_cram_unambiguous_path_empty() {
        // Direct .cram-extension dispatch through run_n2ref. Use an empty CRAM
        // (header-only with no @SQ) to exercise the CRAM dispatch + reader
        // chain without requiring reference resolution.
        use noodles::sam;
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());
        let in_cram = dir.path().join("input.cram");
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            let header: sam::Header = "@HD\tVN:1.6\n".parse().unwrap();
            w.write_header(&header).unwrap();
            w.try_finish(&header).unwrap();
        }
        let out_path = dir.path().join("out.cram");
        run_n2ref(N2RefArgs {
            input: in_cram,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();
        // Output must be a valid CRAM that can be reopened.
        let mut reader = crate::open_cram_reader(&out_path, None).unwrap();
        reader.read_header().unwrap();
    }

    #[test]
    fn test_run_n2ref_extensionless_cram_via_sniff() {
        // Extensionless CRAM exercises the sniff-fallback Cram arm in
        // run_n2ref (lines ~93-106), wiring the FASTA repository through.
        use noodles::sam;
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());
        let cram_path = dir.path().join("real.cram");
        {
            let mut w = crate::open_cram_writer(&cram_path, None).unwrap();
            let header: sam::Header = "@HD\tVN:1.6\n".parse().unwrap();
            w.write_header(&header).unwrap();
            w.try_finish(&header).unwrap();
        }
        let in_path = dir.path().join("piped"); // no extension
        std::fs::rename(&cram_path, &in_path).unwrap();

        let out_path = dir.path().join("out.cram");
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();
        let mut reader = crate::open_cram_reader(&out_path, None).unwrap();
        reader.read_header().unwrap();
    }

    #[test]
    fn test_run_n2ref_extensionless_bam_via_sniff() {
        // Extensionless BAM exercises the sniff-fallback Bam arm in run_n2ref.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());

        let bam_path = dir.path().join("real.bam");
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&bam_path).unwrap();
            w.write_header(&header).unwrap();
            let r = make_mapped_record(b"ANGT", &[30u8; 4], vec![Op::new(Kind::Match, 4)], 0, 1);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let in_path = dir.path().join("piped"); // no extension
        std::fs::rename(&bam_path, &in_path).unwrap();

        let out_path = dir.path().join("out.bam");
        run_n2ref(N2RefArgs {
            input: in_path,
            output: out_path.clone(),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out_path).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs[0].sequence().as_ref().to_vec(), b"ACGT");
    }

    #[test]
    fn test_run_n2ref_fastx_input_errors() {
        // n2ref requires SAM/BAM/CRAM; FASTA/FASTQ input must produce a
        // clear error (not a panic in a downstream reader).
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());
        let bogus = dir.path().join("input"); // sniffs as FASTA via leading '>'
        std::fs::write(&bogus, b">read1\nACGT\n").unwrap();

        let err = run_n2ref(N2RefArgs {
            input: bogus,
            output: dir.path().join("out.sam"),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("requires SAM/BAM/CRAM"));
    }

    #[test]
    fn test_run_n2ref_unknown_format_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());
        let bogus = dir.path().join("input");
        std::fs::write(&bogus, b"\x00\x01garbage").unwrap();

        let err = run_n2ref(N2RefArgs {
            input: bogus,
            output: dir.path().join("out.sam"),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("could not infer format"));
    }

    #[test]
    fn test_run_n2ref_missing_input_file_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = write_ref_fa(dir.path());
        let err = run_n2ref(N2RefArgs {
            input: dir.path().join("does_not_exist.bam"),
            output: dir.path().join("out.bam"),
            reference: fa_path,
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("failed to open"),
            "expected file-open error, got: {msg}"
        );
    }

    #[test]
    fn test_resolve_ref_name_no_reference_sequence_id_errors() {
        // A record with no reference_sequence_id (unmapped equivalent) must
        // produce a clear error from resolve_ref_name.
        let header: noodles::sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        let mut record = RecordBuf::default();
        // alignment_start set but reference_sequence_id is None.
        *record.alignment_start_mut() = Position::new(1);
        let err = resolve_ref_name(&header, &record).unwrap_err();
        assert!(format!("{err:#}").contains("no reference sequence ID"));
    }

    #[test]
    fn test_fetch_ref_seq_no_alignment_start_errors() {
        // A mapped record with reference_sequence_id but no alignment_start
        // is malformed; fetch_ref_seq must fail fast.
        let header: noodles::sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        let mut record = RecordBuf::default();
        *record.reference_sequence_id_mut() = Some(0);

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = {
            let p = dir.path().join("ref.fa");
            std::fs::write(&p, b">chr1\nACGT\n").unwrap();
            std::fs::write(dir.path().join("ref.fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
            p
        };
        let mut ref_reader = open_ref_reader(&fa_path).unwrap();

        let err = fetch_ref_seq(&mut ref_reader, &header, &record).unwrap_err();
        assert!(format!("{err:#}").contains("no alignment start"));
    }

    #[test]
    fn test_fetch_ref_seq_zero_span_returns_empty() {
        // A CIGAR consisting only of soft-clips has zero reference_span;
        // fetch_ref_seq must short-circuit to an empty Vec without querying
        // the FASTA.
        let header: noodles::sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n".parse().unwrap();
        let record =
            make_mapped_record(b"ACGT", &[30u8; 4], vec![Op::new(Kind::SoftClip, 4)], 0, 1);

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = {
            let p = dir.path().join("ref.fa");
            std::fs::write(&p, b">chr1\nACGT\n").unwrap();
            std::fs::write(dir.path().join("ref.fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
            p
        };
        let mut ref_reader = open_ref_reader(&fa_path).unwrap();

        let result = fetch_ref_seq(&mut ref_reader, &header, &record).unwrap();
        assert!(result.is_empty());
    }

    /// `run_n2ref` replaces N-calls with reference bases.
    ///
    /// Skipped when `samtools` is not found on PATH.
    #[test]
    fn test_n2ref_cram() {
        use std::process::Command;

        // Skip if samtools is not available.
        if which::which("samtools").is_err() {
            eprintln!("SKIP: samtools not found; skipping test_n2ref_cram");
            return;
        }

        let dir = tempfile::TempDir::new().unwrap();

        // Write reference: >chr1\nACGT\n  (4 bases)
        let fa_path = dir.path().join("ref.fa");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();

        // Build index with samtools faidx.
        let status = Command::new("samtools")
            .args(["faidx", fa_path.to_str().unwrap()])
            .status()
            .expect("failed to run samtools faidx");
        assert!(status.success(), "samtools faidx failed");

        // Write a SAM file with one mapped record: NNGT at chr1:1, CIGAR=4M.
        let sam_path = dir.path().join("input.sam");
        std::fs::write(
            &sam_path,
            b"@HD\tVN:1.6\n\
              @SQ\tSN:chr1\tLN:4\n\
              r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNGT\tIIII\n",
        )
        .unwrap();

        // Convert to CRAM with samtools (produces a well-formed CRAM file).
        let in_cram = dir.path().join("input.cram");
        let status = Command::new("samtools")
            .args([
                "view",
                "-C",
                "-T",
                fa_path.to_str().unwrap(),
                "-o",
                in_cram.to_str().unwrap(),
                sam_path.to_str().unwrap(),
            ])
            .status()
            .expect("failed to run samtools view");
        assert!(status.success(), "samtools view failed");

        let out_cram = dir.path().join("output.cram");

        run_n2ref(N2RefArgs {
            input: in_cram,
            output: out_cram.clone(),
            reference: fa_path.clone(),
            qual: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        // Read output CRAM with samtools and verify N-calls were replaced.
        // samtools view decodes the CRAM and prints SAM-format lines; we check
        // that the sequence field is "ACGT" (N->A and N->C from the reference).
        let output = Command::new("samtools")
            .args([
                "view",
                "-T",
                fa_path.to_str().unwrap(),
                out_cram.to_str().unwrap(),
            ])
            .output()
            .expect("failed to run samtools view on output");
        assert!(output.status.success(), "samtools view of output failed");

        let stdout = String::from_utf8(output.stdout).unwrap();
        let seq_field = stdout
            .lines()
            .next()
            .and_then(|line| line.split('\t').nth(9))
            .unwrap_or("");
        assert_eq!(
            seq_field, "ACGT",
            "N-calls should be replaced with reference bases"
        );
    }
}
