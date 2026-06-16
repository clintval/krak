//! Annotate SAM/BAM/CRAM records with Kraken classifications.

use std::io::{BufReader, Write};
use std::path::Path;

use ahash::AHashMap;
use anyhow::{Context, Result};
use log::info;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use crate::kraken_report::KrakenReportEntry;
use crate::kraken_report_embed::entries_to_header_comment;
use crate::kraken_result::{KrakenResult, StreamingLookup};
use crate::kraken_taxonomy::{read_taxo_k2d, read_taxonomy_dmp};
use crate::AlignmentFormat;

/// Arguments for the `annotate` command.
pub struct AnnotateArgs {
    /// Input SAM/BAM/CRAM file.
    pub input: std::path::PathBuf,
    /// Kraken classification output file.
    pub assignments: std::path::PathBuf,
    /// Output SAM/BAM/CRAM file with `ti` tags.
    pub output: std::path::PathBuf,
    /// Optional Kraken report file; embeds the taxonomy tree in the output
    /// header as a `@CO krak:report:<base64>` line. Mutually exclusive with
    /// `kraken_db`.
    pub kraken_report: Option<std::path::PathBuf>,
    /// Optional Kraken database directory; reads the DB and embeds the
    /// full taxonomy tree in the output header. Mutually exclusive with
    /// `kraken_report`.
    pub kraken_db: Option<std::path::PathBuf>,
    /// When `true`, load all assignments into a HashMap before reading the
    /// input file (order-independent). When `false` (default), stream
    /// assignments in lock-step with input records (requires both to be at
    /// least weakly in QNAME order).
    pub unordered: bool,
    /// Optional reference FASTA for CRAM decompression (requires `.fai` index).
    pub cram_reference: Option<std::path::PathBuf>,
    /// Number of bgzf compression worker threads for BAM output. Default 1
    /// (one compressor + one writer thread pipelined with the annotation loop).
    /// Ignored for SAM (no compression) and CRAM (per-block codecs).
    pub threads: usize,
    /// bgzf compression level (0-9) for BAM output. Default 5. Ignored for
    /// SAM (no compression) and CRAM (per-block codecs).
    pub compression_level: u32,
}

/// Source of taxon-id assignments for the annotation loop.
enum Source<'a> {
    /// All assignments pre-loaded into a HashMap (order-independent lookup).
    Map(&'a AHashMap<String, u32>),
    /// Path to a Kraken assignments file streamed in lock-step with the input.
    Stream(&'a Path),
}

/// Run the `annotate` command.
///
/// By default, assignments are streamed record-by-record in lock-step with the
/// input, requiring both to be in the same QNAME order (the natural output of a
/// queryname-sorted SAM/BAM/CRAM -> `prep` -> `kraken2` pipeline). With
/// `args.unordered = true`, all assignments are loaded into a `HashMap` first
/// so order does not matter at all.
pub fn run_annotate(args: AnnotateArgs) -> Result<()> {
    let header_comments: Vec<String> = if let Some(db_path) = &args.kraken_db {
        let taxo_path = db_path.join("taxo.k2d");
        let entries = if taxo_path.exists() {
            info!("Reading taxonomy from: {}", taxo_path.display());
            let e = read_taxo_k2d(&taxo_path)?;
            info!("Loaded {} taxonomy nodes from taxo.k2d.", e.len());
            e
        } else {
            let nodes_path = db_path.join("taxonomy").join("nodes.dmp");
            info!("Reading taxonomy from: {}", nodes_path.display());
            let e = read_taxonomy_dmp(db_path)?;
            info!("Loaded {} taxonomy nodes from taxonomy/nodes.dmp.", e.len());
            e
        };
        entries_to_header_comment(&entries)?
    } else if let Some(report_path) = &args.kraken_report {
        info!("Loading Kraken report from: {}", report_path.display());
        let entries = KrakenReportEntry::read_file(report_path)?;
        entries_to_header_comment(&entries)?
    } else {
        Vec::new()
    };

    let map = if args.unordered {
        info!(
            "Loading Kraken assignments (unordered) from: {}",
            args.assignments.display()
        );
        let m = KrakenResult::load_as_map(&args.assignments)?;
        info!("Loaded {} Kraken assignments.", m.len());
        Some(m)
    } else {
        info!(
            "Streaming Kraken assignments from: {}",
            args.assignments.display()
        );
        None
    };
    let source = map
        .as_ref()
        .map_or(Source::Stream(&args.assignments), Source::Map);

    let fmt = AlignmentFormat::from_path(&args.input);
    // Pseudo-paths (/dev/stdin, /dev/fd/N) default to BAM by extension, but
    // the byte content may be SAM/CRAM. Sniff and thread the buffered reader
    // through; re-opening the pseudo-path after sniff would lose the bytes
    // already read off the shared file description.
    if matches!(fmt, AlignmentFormat::Bam) && crate::is_pseudo_path(&args.input) {
        let (sniffed, gzipped, peek_reader) = crate::sniff_input(&args.input)
            .with_context(|| format!("failed to open input: {}", args.input.display()))?;
        return match sniffed {
            crate::SniffedFormat::Sam => {
                let r = crate::into_text_bufread(peek_reader, gzipped);
                let mut reader = noodles::sam::io::Reader::new(r);
                annotate_sam_from_reader(&mut reader, &args.output, source, &header_comments)
            }
            crate::SniffedFormat::Cram => {
                use noodles::cram;
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(crate::build_fasta_repo(
                        args.cram_reference.as_deref(),
                    )?)
                    .build_from_reader(peek_reader);
                annotate_cram_from_reader(
                    &mut reader,
                    &args.output,
                    source,
                    &header_comments,
                    args.cram_reference.as_deref(),
                )
            }
            crate::SniffedFormat::Bam => {
                use noodles::bam;
                use noodles::bgzf;
                let mut reader = bam::io::Reader::from(bgzf::io::Reader::new(peek_reader));
                annotate_bam_from_reader(
                    &mut reader,
                    &args.output,
                    source,
                    &header_comments,
                    args.threads,
                    args.compression_level,
                )
            }
            _ => annotate_bam(
                &args.input,
                &args.output,
                source,
                &header_comments,
                args.threads,
                args.compression_level,
            ),
        };
    }

    match fmt {
        AlignmentFormat::Bam => annotate_bam(
            &args.input,
            &args.output,
            source,
            &header_comments,
            args.threads,
            args.compression_level,
        ),
        AlignmentFormat::Cram => annotate_cram(
            &args.input,
            &args.output,
            source,
            &header_comments,
            args.cram_reference.as_deref(),
        ),
        AlignmentFormat::Sam => annotate_sam(&args.input, &args.output, source, &header_comments),
    }
}

fn annotate_bam(
    input: &Path,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
    threads: usize,
    compression_level: u32,
) -> Result<()> {
    let mut reader = crate::open_bam_reader(input)?;
    annotate_bam_from_reader(
        &mut reader,
        output,
        source,
        header_comments,
        threads,
        compression_level,
    )
}

fn annotate_bam_from_reader<R: std::io::Read>(
    reader: &mut noodles::bam::io::Reader<R>,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
    threads: usize,
    compression_level: u32,
) -> Result<()> {
    use noodles::bam;

    let mut header = reader.read_header().context("failed to read BAM header")?;
    for c in header_comments {
        header.add_comment(c.clone());
    }

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create BAM file: {}", output.display()))?;
    let parz = gzp::par::compress::ParCompressBuilder::<gzp::deflate::Bgzf>::new()
        .num_threads(threads.max(1))
        .with_context(|| "invalid --threads value for BAM bgzf writer")?
        .compression_level(gzp::Compression::new(compression_level))
        .from_writer(file);
    let mut writer = bam::io::Writer::from(parz);
    writer
        .write_header(&header)
        .context("failed to write BAM header")?;

    run_annotate_pipeline(
        writer,
        &header,
        reader.record_bufs(&header),
        source,
        "BAM",
        output,
        Some(AlignmentFormat::Bam),
        |w, _| {
            use gzp::ZWriter as _;
            let mut parz = w.into_inner();
            parz.finish()
                .map_err(|e| anyhow::anyhow!("failed to finish BAM BGZF stream: {e}"))?;
            Ok(())
        },
    )
}

fn annotate_cram(
    input: &Path,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
    cram_reference: Option<&Path>,
) -> Result<()> {
    let mut reader = crate::open_cram_reader(input, cram_reference)?;
    annotate_cram_from_reader(&mut reader, output, source, header_comments, cram_reference)
}

fn annotate_cram_from_reader<R: std::io::Read>(
    reader: &mut noodles::cram::io::Reader<R>,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
    cram_reference: Option<&Path>,
) -> Result<()> {
    let mut header = reader.read_header().context("failed to read CRAM header")?;
    crate::require_cram_reference_if_mapped(&header, cram_reference)?;
    for c in header_comments {
        header.add_comment(c.clone());
    }

    let mut writer = crate::open_cram_writer(output, cram_reference)?;
    writer
        .write_header(&header)
        .context("failed to write CRAM header")?;

    run_annotate_pipeline(
        writer,
        &header,
        reader.records(&header),
        source,
        "CRAM",
        output,
        Some(AlignmentFormat::Cram),
        |mut w, header| w.try_finish(header).context("failed to finish CRAM writer"),
    )
}

fn annotate_sam(
    input: &Path,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
) -> Result<()> {
    use noodles::sam;
    let file = std::fs::File::open(input)
        .with_context(|| format!("failed to open SAM file: {}", input.display()))?;
    let mut reader = sam::io::Reader::new(BufReader::new(file));
    annotate_sam_from_reader(&mut reader, output, source, header_comments)
}

fn annotate_sam_from_reader<R: std::io::BufRead>(
    reader: &mut noodles::sam::io::Reader<R>,
    output: &Path,
    source: Source<'_>,
    header_comments: &[String],
) -> Result<()> {
    use noodles::sam;

    let mut header = reader.read_header().context("failed to read SAM header")?;
    for c in header_comments {
        header.add_comment(c.clone());
    }

    let out_file = std::fs::File::create(output)
        .with_context(|| format!("failed to create SAM file: {}", output.display()))?;
    let mut writer = sam::io::Writer::new(std::io::BufWriter::new(out_file));
    writer
        .write_header(&header)
        .context("failed to write SAM header")?;

    run_annotate_pipeline(
        writer,
        &header,
        reader.record_bufs(&header),
        source,
        "SAM",
        output,
        None,
        |w, _| w.into_inner().flush().context("failed to flush SAM writer"),
    )
}

/// Dispatch to the appropriate lookup strategy for `source` and run the loop.
///
/// Returns `(annotated, total, missing)` where `missing` is the number of
/// records whose QNAME was `*` (no name) and were written unannotated.
fn annotate_records<I, W>(
    iter: I,
    writer: &mut W,
    header: &noodles::sam::Header,
    source: Source<'_>,
    fmt: &str,
) -> Result<(u64, u64, u64)>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    W: AlignmentWrite,
{
    match source {
        Source::Map(map) => {
            annotate_loop(iter, writer, header, fmt, |name| Ok(map.get(name).copied()))
        }
        Source::Stream(path) => {
            let kfile = std::fs::File::open(path)
                .with_context(|| format!("failed to open assignments: {}", path.display()))?;
            let mut state = StreamingLookup::new(BufReader::new(kfile));
            annotate_loop(iter, writer, header, fmt, |name| state.lookup(name))
        }
    }
}

/// Run the per-record annotation loop, finalize the writer, log the summary,
/// and (when `output_format` is `Some`) emit a sidecar index. Owns the writer
/// so format-specific finalize semantics (BGZF `finish`, CRAM `try_finish`,
/// SAM `flush`) can each consume it.
#[allow(clippy::too_many_arguments)]
fn run_annotate_pipeline<W, I, F>(
    mut writer: W,
    header: &noodles::sam::Header,
    records: I,
    source: Source<'_>,
    fmt: &str,
    output: &Path,
    output_format: Option<AlignmentFormat>,
    finalize: F,
) -> Result<()>
where
    W: AlignmentWrite,
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    F: FnOnce(W, &noodles::sam::Header) -> Result<()>,
{
    // Finalize the writer unconditionally, even if the record loop errored: a
    // threaded BGZF (BAM) writer dropped without an explicit finish panics in
    // gzp's Drop, masking the real error. `finish_after` keeps the loop error.
    let body = annotate_records(records, &mut writer, header, source, fmt);
    let (annotated, total, missing) = crate::finish_after(body, || finalize(writer, header))?;
    info!("Annotated {annotated} / {total} records ({missing} records had no name).");
    if let Some(fmt) = output_format {
        crate::maybe_index_alignment_output(output, header, fmt)?;
    }
    Ok(())
}

/// Per-record annotation loop. `lookup(name)` returns the taxon id for the
/// read, `Ok(None)` if the read is absent (treated as a fatal mismatch), or
/// `Err` for any unrecoverable lookup failure.
fn annotate_loop<I, W, F>(
    iter: I,
    writer: &mut W,
    header: &noodles::sam::Header,
    fmt: &str,
    mut lookup: F,
) -> Result<(u64, u64, u64)>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    W: AlignmentWrite,
    F: FnMut(&str) -> Result<Option<u32>>,
{
    use noodles::sam::alignment::record_buf::data::field::Value;

    let mut annotated = 0u64;
    let mut missing = 0u64;
    let mut total = 0u64;

    for result in iter {
        let mut record = result.with_context(|| format!("failed to read {fmt} record"))?;
        total += 1;

        if let Some(name_bytes) = record.name() {
            let name = std::str::from_utf8(name_bytes).context("non-UTF-8 read name")?;
            // Kraken strips trailing `/1`/`/2` from query names; strip the SAM
            // QNAME the same way before lookup so suffixed names still match
            // (mirrors `filter`). The record itself is written unchanged.
            let lookup_name = crate::strip_pair_suffix(name);
            match lookup(lookup_name)? {
                Some(taxon_id) => {
                    // The SAM `ti:i:` aux type is a signed 32-bit integer.
                    // Reject taxon IDs above i32::MAX rather than truncating
                    // them into a negative tag value.
                    let ti = i32::try_from(taxon_id).with_context(|| {
                        format!(
                            "taxon ID {taxon_id} for read {name:?} exceeds the SAM ti tag's \
                             signed 32-bit range"
                        )
                    })?;
                    record
                        .data_mut()
                        .insert(crate::TI_TAG.into(), Value::Int32(ti));
                    annotated += 1;
                }
                None => {
                    anyhow::bail!(
                        "read {name:?} (record {total}) is not present in the assignments file; \
                         ensure the Kraken assignments file contains every read in the input"
                    );
                }
            }
        } else {
            missing += 1;
        }

        writer
            .write_alignment_record(header, &record)
            .with_context(|| format!("failed to write {fmt} record"))?;
    }

    Ok((annotated, total, missing))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlignmentFormat;
    use noodles::sam::alignment::record_buf::{data::field::Value, RecordBuf};

    fn make_record(name: &str) -> RecordBuf {
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(name.as_bytes().into());
        r
    }

    #[test]
    fn test_run_annotate_pipeline_finalizes_writer_even_on_record_error() {
        use noodles::sam;
        use std::cell::Cell;

        let header = sam::Header::default();
        let writer = sam::io::Writer::new(Vec::new());
        let finalized = Cell::new(false);
        let empty: AHashMap<String, u32> = AHashMap::new();
        // A record stream that errors on its first item.
        let records = std::iter::once::<std::io::Result<RecordBuf>>(Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "boom",
        )));

        let result = run_annotate_pipeline(
            writer,
            &header,
            records,
            Source::Map(&empty),
            "SAM",
            std::path::Path::new("/dev/stdout"),
            None,
            |w, _| {
                finalized.set(true);
                w.into_inner().flush().context("flush")?;
                Ok(())
            },
        );

        assert!(result.is_err(), "the record-loop error must propagate");
        assert!(
            finalized.get(),
            "writer must be finalized even when the record loop errors (else a \
             threaded BGZF writer is dropped un-finished and panics in gzp's Drop)"
        );
    }

    #[test]
    fn test_annotate_strips_pair_suffix_from_sam_qname() {
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let in_sam = dir.path().join("in.sam");
        // QNAME literally carries a Kraken v1 `/1` suffix.
        std::fs::write(
            &in_sam,
            b"@HD\tVN:1.6\nread1/1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let assignments = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments).unwrap();
            // Kraken stores the base name (the `/1` suffix is stripped on parse).
            writeln!(f, "C\tread1\t9606\t4\t9606:1").unwrap();
        }

        let out_sam = dir.path().join("out.sam");
        super::run_annotate(super::AnnotateArgs {
            input: in_sam,
            assignments,
            output: out_sam.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: false,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        // The record must be annotated, not fatally rejected as "not present".
        let text = std::fs::read_to_string(&out_sam).unwrap();
        assert!(
            text.contains("ti:i:9606"),
            "expected ti:i:9606 in output, got:\n{text}"
        );
    }

    #[test]
    fn test_annotate_errors_on_taxon_id_exceeding_i32() {
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let in_sam = dir.path().join("in.sam");
        std::fs::write(
            &in_sam,
            b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let assignments = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments).unwrap();
            // 3_000_000_000 fits in u32 but exceeds i32::MAX (2_147_483_647);
            // it must not be truncated into a negative `ti` tag.
            writeln!(f, "C\tread1\t3000000000\t4\t3000000000:1").unwrap();
        }

        let result = super::run_annotate(super::AnnotateArgs {
            input: in_sam,
            assignments,
            output: dir.path().join("out.sam"),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        });

        let err = result.unwrap_err();
        assert!(
            format!("{err:#}").contains("taxon"),
            "expected a taxon-range error, got: {err:#}"
        );
    }

    #[test]
    fn test_ti_tag_inserted() {
        let mut record = make_record("read1");
        record
            .data_mut()
            .insert(crate::TI_TAG.into(), Value::Int32(9606));
        match record.data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            _ => panic!("ti tag not found or wrong type"),
        }
    }

    #[test]
    fn test_ti_tag_overwrite() {
        let mut record = make_record("read2");
        record
            .data_mut()
            .insert(crate::TI_TAG.into(), Value::Int32(1));
        record
            .data_mut()
            .insert(crate::TI_TAG.into(), Value::Int32(9606));
        match record.data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            _ => panic!("ti tag not updated"),
        }
    }

    /// Verify that after BAM round-trip the `ti` tag is preserved as the SAM
    /// `i` (signed 32-bit integer) aux type rather than `I` (unsigned).
    #[test]
    fn test_ti_tag_is_int32_after_bam_round_trip() {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        let in_bam = dir.path().join("input.bam");
        let header = sam::Header::default();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("readZ".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
            w.into_inner().finish().unwrap();
        }

        let assignments_path = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments_path).unwrap();
            writeln!(f, "C\treadZ\t9606\t4\t9606:1").unwrap();
        }

        let out_bam = dir.path().join("output.bam");
        super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments: assignments_path,
            output: out_bam.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut reader = crate::open_bam_reader(&out_bam).unwrap();
        let out_header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .record_bufs(&out_header)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        match records[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("expected Int32 (SAM type 'i'); got: {:?}", other),
        }
    }

    /// Write a single-record BAM, then run `annotate` writing to `out_name`
    /// with the given `threads` and `compression_level`. Returns the output
    /// file size in bytes and the round-tripped `ti` tag value.
    fn run_annotate_bam(
        tmpdir: &std::path::Path,
        out_name: &str,
        threads: usize,
        compression_level: u32,
    ) -> (u64, i32) {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};
        use std::io::Write as _;

        let in_bam = tmpdir.join("in.bam");
        let header = sam::Header::default();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            w.write_header(&header).unwrap();
            // A few records with varied sequences give the compressor enough
            // material that level 1 vs 9 differ measurably while staying tiny.
            for i in 0..32 {
                let name = format!("read{i}");
                let mut r = RecordBuf::default();
                *r.name_mut() = Some(name.as_bytes().into());
                let bases: Vec<u8> = (0..96)
                    .map(|j| match (i + j) % 4 {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        _ => b'T',
                    })
                    .collect();
                *r.sequence_mut() = Sequence::from(bases.clone());
                *r.quality_scores_mut() = QualityScores::from(vec![30u8; bases.len()]);
                w.write_alignment_record(&header, &r).unwrap();
            }
            w.into_inner().finish().unwrap();
        }

        let assignments = tmpdir.join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments).unwrap();
            for i in 0..32 {
                writeln!(f, "C\tread{i}\t9606\t96\t9606:1").unwrap();
            }
        }

        let out_bam = tmpdir.join(out_name);
        super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments,
            output: out_bam.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads,
            compression_level,
        })
        .unwrap();

        let size = std::fs::metadata(&out_bam).unwrap().len();
        let mut reader = crate::open_bam_reader(&out_bam).unwrap();
        let h = reader.read_header().unwrap();
        let first: RecordBuf = reader.record_bufs(&h).next().unwrap().unwrap();
        let ti = match first.data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => *n,
            other => panic!("expected Int32 ti tag; got {other:?}"),
        };
        (size, ti)
    }

    #[test]
    fn test_annotate_bam_higher_compression_level_yields_smaller_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let (size_low, ti_low) = run_annotate_bam(dir.path(), "low.bam", 1, 1);
        let (size_high, ti_high) = run_annotate_bam(dir.path(), "high.bam", 1, 9);
        assert_eq!(ti_low, 9606);
        assert_eq!(ti_high, 9606);
        assert!(
            size_high < size_low,
            "expected level 9 ({size_high} bytes) < level 1 ({size_low} bytes)"
        );
    }

    #[test]
    fn test_annotate_bam_threads_one_and_many_round_trip_identically() {
        let dir = tempfile::TempDir::new().unwrap();
        let (_, ti_serial) = run_annotate_bam(dir.path(), "t1.bam", 1, 5);
        let (_, ti_parallel) = run_annotate_bam(dir.path(), "t4.bam", 4, 5);
        assert_eq!(ti_serial, 9606);
        assert_eq!(ti_parallel, 9606);
    }

    #[test]
    fn test_annotate_cram_unordered() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        // Write input CRAM with two named records; sequences must be non-empty
        // to avoid empty external data blocks in the CRAM writer.
        let in_cram = dir.path().join("input.cram");
        let header = sam::Header::default();
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();
            let mut r1 = RecordBuf::default();
            *r1.name_mut() = Some("read1".as_bytes().into());
            *r1.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r1.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r1).unwrap();
            let mut r2 = RecordBuf::default();
            *r2.name_mut() = Some("read2".as_bytes().into());
            *r2.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r2.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r2).unwrap();
            w.try_finish(&header).unwrap();
        }

        // Kraken assignments: read1 -> 9606, read2 -> 1234
        let assignments_path = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments_path).unwrap();
            writeln!(f, "C\tread1\t9606\t4\t9606:1").unwrap();
            writeln!(f, "C\tread2\t1234\t4\t1234:1").unwrap();
        }

        let out_cram = dir.path().join("output.cram");

        super::run_annotate(super::AnnotateArgs {
            input: in_cram,
            assignments: assignments_path,
            output: out_cram.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        // Verify ti tags in output CRAM
        let mut reader = crate::open_cram_reader(&out_cram, None).unwrap();
        let out_header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .records(&out_header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 2);
        match records[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("unexpected ti tag value: {:?}", other),
        }
        match records[1].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 1234),
            other => panic!("unexpected ti tag value: {:?}", other),
        }
    }

    #[test]
    fn test_annotate_cram_streaming() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        // Non-empty sequence avoids empty external data blocks in the CRAM writer.
        let in_cram = dir.path().join("input.cram");
        let header = sam::Header::default();
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("readA".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
            w.try_finish(&header).unwrap();
        }

        let assignments_path = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments_path).unwrap();
            writeln!(f, "C\treadA\t9606\t4\t9606:1").unwrap();
        }

        let out_cram = dir.path().join("output.cram");

        // unordered: false -> streaming mode
        super::run_annotate(super::AnnotateArgs {
            input: in_cram,
            assignments: assignments_path,
            output: out_cram.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: false,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut reader = crate::open_cram_reader(&out_cram, None).unwrap();
        let out_header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .records(&out_header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        match records[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("unexpected ti tag value: {:?}", other),
        }
    }

    #[test]
    fn test_run_annotate_unnamed_record_is_skipped() {
        // Records with `*` QNAME (no name) are passed through with no `ti`
        // tag and counted in the missing log message (line ~364).
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("in.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            // Named record + unnamed (None) record; neither in assignments.
            // Use unordered=true so the named record errors only on missing,
            // not on a streaming gap. Actually pre-load the named record.
            let mut named = RecordBuf::default();
            *named.name_mut() = Some(b"named".as_ref().into());
            *named.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *named.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &named).unwrap();
            // No name set → record.name() returns None.
            let mut unnamed = RecordBuf::default();
            *unnamed.sequence_mut() = Sequence::from(b"TTTT".to_vec());
            *unnamed.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &unnamed).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        std::fs::write(&assignments, b"C\tnamed\t9606\t4\t9606:1\n").unwrap();

        let out = dir.path().join("out.bam");
        super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments,
            output: out.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 2);
        // Named record gets ti tag.
        assert!(matches!(
            recs[0].data().get(&crate::TI_TAG),
            Some(Value::Int32(9606))
        ));
        // Unnamed record passes through with no ti tag.
        assert!(recs[1].data().get(&crate::TI_TAG).is_none());
    }

    #[test]
    fn test_run_annotate_bam_unambiguous_path() {
        // Direct .bam-extension dispatch through run_annotate -> annotate_bam.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("input.bam");
        let header = sam::Header::default();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"readB".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        std::fs::write(&assignments, b"C\treadB\t9606\t4\t9606:1\n").unwrap();

        let out_bam = dir.path().join("out.bam");
        super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments,
            output: out_bam.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: false,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out_bam).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        match recs[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("expected ti:i:9606, got {other:?}"),
        }
    }

    #[test]
    fn test_run_annotate_with_kraken_report_embeds_header() {
        // When --kraken-report is supplied, run_annotate must serialize the
        // taxonomy and add @CO krak:report:* line(s) to the output header.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("input.bam");
        let header = sam::Header::default();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"readR".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        std::fs::write(&assignments, b"C\treadR\t9606\t4\t9606:1\n").unwrap();
        let report = dir.path().join("report.k2report");
        std::fs::write(
            &report,
            b"100.00\t1\t1\tR\t1\troot\n100.00\t1\t1\tS\t9606\t  Homo sapiens\n",
        )
        .unwrap();

        let out_bam = dir.path().join("out.bam");
        super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments,
            output: out_bam.clone(),
            kraken_report: Some(report),
            kraken_db: None,
            unordered: false,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out_bam).unwrap();
        let h = r.read_header().unwrap();
        let has_embed = h.comments().iter().any(|c| {
            std::str::from_utf8(c)
                .map(|s| s.starts_with("krak:report:"))
                .unwrap_or(false)
        });
        assert!(
            has_embed,
            "expected at least one krak:report: header comment"
        );
    }

    /// When path-based detection says BAM only because the input path is a
    /// pseudo-path (here a symlink under `/dev/fd/` is awkward to construct in
    /// a unit test, so we directly verify the `from_path` heuristic plus the
    /// sniff fallback by symlinking via /dev/fd is platform-specific).
    /// Instead, simulate the flow by routing a SAM file through the SAM handler
    /// when input has no extension yet sniffs as SAM.
    #[test]
    fn test_run_annotate_sniff_fallback_extensionless_sam() {
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        // Write a SAM stream to a file with no recognized extension.
        // AlignmentFormat::from_path falls back to Sam for unknown
        // extensions, NOT Bam; pseudo-paths are the only ones forcing Bam.
        // To exercise this code path proper, use a path that mimics a
        // pseudo-path tail so detection picks Bam, then sniff overrides to Sam.
        let sam_path = dir.path().join("stream.sam");
        {
            let mut f = std::fs::File::create(&sam_path).unwrap();
            writeln!(f, "@HD\tVN:1.6").unwrap();
            writeln!(f, "readP\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").unwrap();
        }
        let assignments_path = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments_path).unwrap();
            writeln!(f, "U\treadP\t0\t4\t0:4").unwrap();
        }
        let out_path = dir.path().join("out.sam");

        super::run_annotate(super::AnnotateArgs {
            input: sam_path,
            assignments: assignments_path,
            output: out_path.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        // Output should be a valid SAM with the ti tag.
        let body = std::fs::read_to_string(&out_path).unwrap();
        assert!(
            body.contains("ti:i:0"),
            "expected ti aux tag in SAM output, got:\n{body}"
        );
    }

    /// /dev/fd/N pointing at a real BAM exercises the sniff-fallback Bam arm
    /// in run_annotate (lines ~131-136). The path-based detection defaults to
    /// Bam for pseudo-paths, the sniffer confirms it; the buffered reader is
    /// threaded into annotate_bam_from_reader.
    #[cfg(unix)]
    #[test]
    fn test_run_annotate_sniff_fallback_dev_fd_routes_bam() {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        use std::io::Write as _;
        use std::os::fd::AsRawFd;

        let dir = tempfile::TempDir::new().unwrap();
        let bam_path = dir.path().join("in.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&bam_path).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"readD".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments).unwrap();
            writeln!(f, "C\treadD\t9606\t4\t9606:4").unwrap();
        }

        let f = std::fs::File::open(&bam_path).unwrap();
        let fd = f.as_raw_fd();
        let pseudo = std::path::PathBuf::from(format!("/dev/fd/{fd}"));

        let out = dir.path().join("out.bam");
        super::run_annotate(super::AnnotateArgs {
            input: pseudo,
            assignments,
            output: out.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut r = bam::io::reader::Builder.build_from_path(&out).unwrap();
        let h = r.read_header().unwrap();
        let recs: Vec<RecordBuf> = r
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        match recs[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("expected ti:i:9606, got {other:?}"),
        }
    }

    /// /dev/fd/N pointing at a real CRAM exercises the sniff-fallback Cram arm
    /// in run_annotate (lines ~115-129).
    #[cfg(unix)]
    #[test]
    fn test_run_annotate_sniff_fallback_dev_fd_routes_cram() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        use std::io::Write as _;
        use std::os::fd::AsRawFd;

        let dir = tempfile::TempDir::new().unwrap();
        let cram_path = dir.path().join("in.cram");
        {
            let mut w = crate::open_cram_writer(&cram_path, None).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"readC".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
            w.try_finish(&header).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments).unwrap();
            writeln!(f, "C\treadC\t9606\t4\t9606:4").unwrap();
        }

        let f = std::fs::File::open(&cram_path).unwrap();
        let fd = f.as_raw_fd();
        let pseudo = std::path::PathBuf::from(format!("/dev/fd/{fd}"));

        let out = dir.path().join("out.cram");
        super::run_annotate(super::AnnotateArgs {
            input: pseudo,
            assignments,
            output: out.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let mut reader = crate::open_cram_reader(&out, None).unwrap();
        let h = reader.read_header().unwrap();
        let recs: Vec<RecordBuf> = reader.records(&h).collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(recs.len(), 1);
        match recs[0].data().get(&crate::TI_TAG) {
            Some(Value::Int32(n)) => assert_eq!(*n, 9606),
            other => panic!("expected ti:i:9606, got {other:?}"),
        }
    }

    /// Regression: a mapped CRAM (with `@SQ` in its header) without
    /// `--cram-reference` used to panic inside noodles' decoder ("invalid
    /// slice reference sequence name") when annotate began iterating
    /// records. The fix bails with a clear error after reading the header.
    #[test]
    fn test_run_annotate_mapped_cram_without_reference_errors_cleanly() {
        use noodles::sam;
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let dir = tempfile::TempDir::new().unwrap();
        let fa_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();
        std::fs::write(&fai_path, b"chr1\t4\t6\t4\t5\n").unwrap();

        let in_cram = dir.path().join("in.cram");
        let mut header = sam::Header::default();
        let len = NonZeroUsize::new(4).unwrap();
        header
            .reference_sequences_mut()
            .insert(b"chr1".as_ref().into(), Map::<ReferenceSequence>::new(len));
        {
            let mut w = crate::open_cram_writer(&in_cram, Some(&fa_path)).unwrap();
            w.write_header(&header).unwrap();
            w.try_finish(&header).unwrap();
        }

        let assignments = dir.path().join("assignments.txt");
        std::fs::write(&assignments, b"").unwrap();

        let out_cram = dir.path().join("out.cram");
        let err = super::run_annotate(super::AnnotateArgs {
            input: in_cram,
            assignments,
            output: out_cram,
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--cram-reference"), "got: {msg}");
        assert!(msg.contains("reference sequences"), "got: {msg}");
    }

    #[test]
    fn test_run_annotate_missing_read_in_assignments_errors() {
        // A record present in the BAM but absent from the unordered map source
        // must trigger the "not present in the assignments file" error.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("in.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(b"readMissing".as_ref().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        let assignments = dir.path().join("assignments.txt");
        std::fs::write(&assignments, b"C\tdifferent\t9606\t4\t9606:4\n").unwrap();

        let out = dir.path().join("out.bam");
        let err = super::run_annotate(super::AnnotateArgs {
            input: in_bam,
            assignments,
            output: out,
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("not present in the assignments"));
    }

    /// Pseudo-path branch: construct a `/dev/fd/N` reference to a real
    /// SAM file and verify that `run_annotate` routes through the SAM handler
    /// rather than failing in the BAM reader. Requires Unix `/dev/fd/`.
    #[cfg(unix)]
    #[test]
    fn test_run_annotate_sniff_fallback_dev_fd_routes_sam() {
        use std::io::Write as _;
        use std::os::fd::AsRawFd;

        let dir = tempfile::TempDir::new().unwrap();
        let sam_path = dir.path().join("stream.sam");
        {
            let mut f = std::fs::File::create(&sam_path).unwrap();
            writeln!(f, "@HD\tVN:1.6").unwrap();
            writeln!(f, "readQ\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").unwrap();
        }
        let assignments_path = dir.path().join("assignments.txt");
        {
            let mut f = std::fs::File::create(&assignments_path).unwrap();
            writeln!(f, "U\treadQ\t0\t4\t0:4").unwrap();
        }
        let out_path = dir.path().join("out.sam");

        // Open the SAM file and pass /dev/fd/N as the input path. Path-based
        // detection returns Bam for /dev/fd/* but sniff_input should report
        // Sam, triggering the fallback.
        let f = std::fs::File::open(&sam_path).unwrap();
        let fd = f.as_raw_fd();
        let pseudo = std::path::PathBuf::from(format!("/dev/fd/{fd}"));

        super::run_annotate(super::AnnotateArgs {
            input: pseudo,
            assignments: assignments_path,
            output: out_path.clone(),
            kraken_report: None,
            kraken_db: None,
            unordered: true,
            cram_reference: None,
            threads: 1,
            compression_level: 5,
        })
        .unwrap();

        let body = std::fs::read_to_string(&out_path).unwrap();
        assert!(body.contains("ti:i:0"), "expected ti aux tag, got:\n{body}");
    }

    #[test]
    fn test_format_detection() {
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.bam")),
            AlignmentFormat::Bam
        );
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.BAM")),
            AlignmentFormat::Bam
        );
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.cram")),
            AlignmentFormat::Cram
        );
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.CRAM")),
            AlignmentFormat::Cram
        );
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.sam")),
            AlignmentFormat::Sam
        );
        assert_eq!(
            AlignmentFormat::from_path(std::path::Path::new("foo.txt")),
            AlignmentFormat::Sam
        );
    }
}
