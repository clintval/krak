//! Filter FASTX/SAM/BAM/CRAM records by Kraken classifications.

use std::io::{BufReader, Write};

use ahash::{AHashMap, AHashSet};
use std::path::Path;

use anyhow::{Context, Result};
use log::{info, warn};
use noodles::sam::alignment::io::Write as AlignmentWrite;
use serde::Serialize;

use crate::kraken_report::KrakenTaxonomyTree;
use crate::kraken_result::{KrakenResult, StreamingLookup};
use crate::{AlignmentFormat, FastxKind, InferredFormat};

/// Expand the initial taxon ID set with descendants and/or taxon 0.
///
/// - If `include_descendants` is true and `tree` is provided, every descendant
///   of each target taxon is inserted.
/// - If `include_unclassified` is true, taxon ID `0` is inserted.
///
/// When a `tree` is provided, every input target (other than the synthetic 0
/// taxon) must exist in the tree; an unknown ID is an error. Without a tree,
/// targets are passed through unchecked.
fn build_taxon_set(
    base: &AHashSet<u32>,
    tree: Option<&KrakenTaxonomyTree>,
    include_descendants: bool,
    include_unclassified: bool,
) -> Result<AHashSet<u32>> {
    let mut set = base.clone();
    if let Some(tree) = tree {
        for taxon in base.iter().copied() {
            if taxon == 0 {
                continue;
            }
            if !tree.contains(taxon) {
                anyhow::bail!(
                    "taxon id {taxon} from --taxon-id is not present in the taxonomy tree"
                );
            }
        }
        if include_descendants {
            for taxon in base.iter().copied() {
                set.extend(tree.descendants_of(taxon));
            }
        }
    }
    if include_unclassified {
        set.insert(0);
    }
    Ok(set)
}

/// Arguments for the `filter` command.
#[derive(Debug, Clone)]
pub struct FilterArgs {
    /// Input SAM/BAM/CRAM file annotated with `ti` tags.
    pub input: std::path::PathBuf,
    /// Output SAM/BAM/CRAM for passing records.
    pub output: std::path::PathBuf,
    /// Kraken report file; fallback when no report is embedded in the SAM/BAM/CRAM header.
    pub kraken_report: Option<std::path::PathBuf>,
    /// TSV metrics output file. If `None`, metrics are only logged.
    pub metrics: Option<std::path::PathBuf>,
    /// Target taxon IDs to keep.
    pub taxon_ids: AHashSet<u32>,
    /// Optional output for rejected records.
    pub rejects: Option<std::path::PathBuf>,
    /// Keep reads assigned to ancestors of target taxa.
    pub allow_ancestors: bool,
    /// Maximum edit distance from reference for off-taxa rescue.
    /// `None` means edit-distance rescue is not requested.
    pub rescue_max_edit_distance: Option<u32>,
    /// Maximum number of indel events for off-taxa rescue.
    pub rescue_max_indels: Option<u32>,
    /// Maximum length of any single indel for off-taxa rescue.
    pub rescue_max_indel_length: Option<u32>,
    /// Reduce edit distance limit by 1 per `rescue_n_adjustment` Ns in a read.
    pub rescue_n_adjustment: Option<u32>,
    /// Process each alignment record independently, bypassing template grouping.
    /// Required when input is not query-grouped/queryname-sorted.
    pub per_record: bool,
    /// Kraken2 per-read classification file. Required with FASTA/FASTQ input.
    pub classifications: Option<std::path::PathBuf>,
    /// Expand target taxa to include all descendants in the clade.
    /// Requires --kraken-report (or an embedded report in the SAM/BAM/CRAM header).
    pub include_descendants: bool,
    /// Also keep reads classified as taxon 0 (unclassified).
    pub include_unclassified: bool,
    /// Optional reference FASTA for CRAM decompression (requires `.fai` index).
    pub cram_reference: Option<std::path::PathBuf>,
    /// When `true`, records that lack a `ti` tag (no Kraken assignment) are kept.
    /// When `false` (default), such records are rejected.
    /// Only applies to SAM/BAM/CRAM; FASTX uses --classifications and unmapped
    /// reads cannot be 'unannotated'.
    pub keep_unannotated: bool,
    /// When `true`, load all assignments into a HashMap before reading the
    /// FASTA/FASTQ input (order-independent). When `false` (default), stream
    /// assignments in lock-step with input records using a self-sizing
    /// lookahead buffer (handles modest disorder such as multi-threaded
    /// Kraken v1 output). Only applies to FASTA/FASTQ input; SAM/BAM/CRAM
    /// input reads taxon IDs from `ti` tags.
    pub unordered: bool,
    /// Number of bgzf compression worker threads for `.gz` outputs.
    /// At `1` (default), one compressor + one writer thread pipeline with the
    /// main filter loop. Higher values fan compression out across more workers.
    /// SAM and CRAM outputs ignore this value.
    pub threads: usize,
    /// bgzf compression level (0-9) for `.gz` outputs. Default 5.
    pub compression_level: u32,
}

/// MD tag.
const MD_TAG: [u8; 2] = [b'M', b'D'];

/// Counts of edits from reference for a single record.
#[derive(Debug, Default)]
struct EditCounts {
    /// Number of substitutions (Ns excluded).
    substitutions: i64,
    /// Signed lengths of indel events: positive = insertion, negative = deletion.
    indels: Vec<i32>,
}

impl EditCounts {
    /// Total edit events (substitutions + number of indel events).
    fn total(&self) -> i64 {
        self.substitutions + self.indels.len() as i64
    }

    /// Number of indel events.
    fn indel_count(&self) -> usize {
        self.indels.len()
    }

    /// Maximum absolute indel length, or 0 if no indels.
    fn max_indel_length(&self) -> u32 {
        self.indels
            .iter()
            .map(|l| l.unsigned_abs())
            .max()
            .unwrap_or(0)
    }
}

/// Metrics for a single filter pass.
#[derive(Debug, Default, Serialize)]
pub struct TaxaFilterMetric {
    /// `true` if counts are per-template, `false` if per-record.
    pub template: bool,
    /// Records/templates that matched target taxon IDs exactly.
    pub on_taxa: u64,
    /// Records/templates rescued as ancestors of target taxa.
    pub rescued_ancestors: u64,
    /// Records/templates rescued by edit-distance proximity.
    /// `None` for FASTA/FASTQ mode (no MD/CIGAR available).
    pub rescued_variants: Option<u64>,
    /// Records/templates kept solely because they lacked a `ti` tag and
    /// `--keep-unannotated` was active.
    /// `None` for FASTA/FASTQ mode (no SAM tags).
    pub unannotated_but_kept: Option<u64>,
    /// Total kept.
    pub num_kept: u64,
    /// Total rejected.
    pub num_filtered: u64,
    /// Total processed.
    pub total: u64,
    /// Fraction rejected.
    pub frac_removed: f64,
    /// `--rescue-max-edit-distance` parameter used.
    pub rescue_max_edit_distance: Option<u32>,
    /// `--rescue-max-indels` parameter used.
    pub rescue_max_indels: Option<u32>,
    /// `--rescue-max-indel-length` parameter used.
    pub rescue_max_indel_length: Option<u32>,
    /// `--rescue-n-adjustment` parameter used.
    pub rescue_n_adjustment: Option<u32>,
    /// Records/templates missing the `ti` tag.
    /// `None` for FASTA/FASTQ mode (no SAM tags).
    pub missing_ti_tag: Option<u64>,
    /// Reads whose name was absent from the `--classifications` map.
    /// `None` for SAM/BAM/CRAM mode.
    pub missing_classifications: Option<u64>,
    /// Whether `--keep-unannotated` was active, suppressing the "no ti tag" warning.
    /// Not a metric column; excluded from TSV output.
    #[serde(skip)]
    pub keep_unannotated: bool,
}

impl TaxaFilterMetric {
    fn finalize(&mut self) {
        self.num_kept = self.on_taxa
            + self.rescued_ancestors
            + self.rescued_variants.unwrap_or(0)
            + self.unannotated_but_kept.unwrap_or(0);
        self.num_filtered = self.total - self.num_kept;
        self.frac_removed = if self.total == 0 {
            0.0
        } else {
            self.num_filtered as f64 / self.total as f64
        };
    }
}

/// Run the `filter` command.
pub fn run_filter(args: FilterArgs) -> Result<()> {
    if matches!(args.rescue_n_adjustment, Some(0)) {
        anyhow::bail!("--rescue-n-adjustment must be >= 1, or unset");
    }
    let metrics = run_filter_dispatch(&args)?;

    if let Some(path) = &args.metrics {
        write_metrics(path, &metrics)?;
    }
    log_metrics(&metrics);
    Ok(())
}

/// Dispatch the filter command to the right per-format handler. When the
/// path's extension is ambiguous (no `.fa`/`.fq`/`.sam`/`.bam`/`.cram` and
/// peers), defer to the byte sniffer.
fn run_filter_dispatch(args: &FilterArgs) -> Result<Vec<TaxaFilterMetric>> {
    // Fast path: extension says FASTX.
    if matches!(crate::infer_format(&args.input), InferredFormat::Fastx(_)) {
        if args.classifications.is_none() {
            anyhow::bail!("FASTA/FASTQ input requires --classifications (-c)");
        }
        return Ok(vec![filter_fastx_runner(args)?]);
    }

    let classifications_for_alignment_err = || {
        anyhow::anyhow!(
            "--classifications (-c) is not valid for SAM/BAM/CRAM input; \
             taxon IDs are read from the `ti` tag"
        )
    };

    // Extension says Sam (or unknown). Use the existing alignment dispatch
    // when the extension is unambiguous (.sam/.bam/.cram).
    let ext = args
        .input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase());
    let unambiguous_alignment = matches!(ext.as_deref(), Some("sam") | Some("bam") | Some("cram"));
    if unambiguous_alignment {
        if args.classifications.is_some() {
            return Err(classifications_for_alignment_err());
        }
        return match AlignmentFormat::from_path(&args.input) {
            AlignmentFormat::Bam => filter_bam(args),
            AlignmentFormat::Cram => filter_cram(args),
            AlignmentFormat::Sam => filter_sam(args),
        };
    }

    // Sniff. The reader returned holds the bytes already peeked.
    let (sniffed, gzipped, reader) = crate::sniff_input(&args.input)
        .with_context(|| format!("failed to open input: {}", args.input.display()))?;
    if matches!(
        sniffed,
        crate::SniffedFormat::Sam | crate::SniffedFormat::Bam | crate::SniffedFormat::Cram
    ) && args.classifications.is_some()
    {
        return Err(classifications_for_alignment_err());
    }
    match sniffed {
        crate::SniffedFormat::Fasta | crate::SniffedFormat::Fastq => {
            if args.classifications.is_none() {
                anyhow::bail!("FASTA/FASTQ input requires --classifications (-c)");
            }
            let r = crate::into_text_bufread(reader, gzipped);
            Ok(vec![filter_fastx_from_reader(args, sniffed, r)?])
        }
        crate::SniffedFormat::Sam => {
            use noodles::sam;
            let r = crate::into_text_bufread(reader, gzipped);
            let reader = sam::io::Reader::new(r);
            filter_sam_with_reader(args, reader)
        }
        crate::SniffedFormat::Bam => {
            use noodles::bam;
            use noodles::bgzf;
            let reader = bam::io::Reader::from(bgzf::io::Reader::new(reader));
            filter_bam_with_reader(args, reader)
        }
        crate::SniffedFormat::Cram => {
            use noodles::cram;
            let reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(crate::build_fasta_repo(
                    args.cram_reference.as_deref(),
                )?)
                .build_from_reader(reader);
            filter_cram_with_reader(args, reader)
        }
        crate::SniffedFormat::Unknown => anyhow::bail!(
            "could not infer format from input head bytes for {}; \
             supply a file with a known extension",
            args.input.display()
        ),
    }
}

/// Path-based FASTX entry point: builds the taxonomy tree, expanded taxon
/// set, filter params, and chooses whether to load assignments into a map
/// (`--unordered`) or stream them with a lookahead buffer (default), then
/// dispatches to FASTQ or FASTA based on the file extension.
fn filter_fastx_runner(args: &FilterArgs) -> Result<TaxaFilterMetric> {
    let classifications_path = args
        .classifications
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("FASTX mode requires --classifications (-c)"))?;

    let tree: Option<KrakenTaxonomyTree> = args
        .kraken_report
        .as_deref()
        .map(KrakenTaxonomyTree::from_file)
        .transpose()?;
    validate_tree_requirements(&tree, args)?;
    let expanded = build_taxon_set(
        &args.taxon_ids,
        tree.as_ref(),
        args.include_descendants,
        args.include_unclassified,
    )?;
    let params = FilterParams {
        taxon_ids: &expanded,
        tree: tree.as_ref(),
        allow_ancestors: args.allow_ancestors,
        rescue_max_edit_distance: args.rescue_max_edit_distance,
        rescue_max_indels: args.rescue_max_indels,
        rescue_max_indel_length: args.rescue_max_indel_length,
        rescue_n_adjustment: args.rescue_n_adjustment,
        keep_unannotated: args.keep_unannotated,
    };

    let map = load_kraken_map_if_unordered(args, classifications_path)?;
    let source = make_fastx_source(map.as_ref(), classifications_path);

    match crate::infer_format(&args.input) {
        InferredFormat::Fastx(FastxKind::Fastq) => filter_fastq_path(args, &params, source),
        _ => filter_fasta_path(args, &params, source),
    }
}

/// Sniffed-FASTX entry point: build params + tree like the path-based runner,
/// but feed the source and reader through `filter_*_from_reader`.
fn filter_fastx_from_reader(
    args: &FilterArgs,
    sniffed: crate::SniffedFormat,
    reader: Box<dyn std::io::BufRead>,
) -> Result<TaxaFilterMetric> {
    let classifications_path = args
        .classifications
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("FASTX mode requires --classifications (-c)"))?;

    let tree: Option<KrakenTaxonomyTree> = args
        .kraken_report
        .as_deref()
        .map(KrakenTaxonomyTree::from_file)
        .transpose()?;
    validate_tree_requirements(&tree, args)?;
    let expanded = build_taxon_set(
        &args.taxon_ids,
        tree.as_ref(),
        args.include_descendants,
        args.include_unclassified,
    )?;
    let params = FilterParams {
        taxon_ids: &expanded,
        tree: tree.as_ref(),
        allow_ancestors: args.allow_ancestors,
        rescue_max_edit_distance: args.rescue_max_edit_distance,
        rescue_max_indels: args.rescue_max_indels,
        rescue_max_indel_length: args.rescue_max_indel_length,
        rescue_n_adjustment: args.rescue_n_adjustment,
        keep_unannotated: args.keep_unannotated,
    };

    let map = load_kraken_map_if_unordered(args, classifications_path)?;
    let source = make_fastx_source(map.as_ref(), classifications_path);

    match sniffed {
        crate::SniffedFormat::Fastq => filter_fastq_from_reader(args, &params, source, reader),
        crate::SniffedFormat::Fasta => filter_fasta_from_reader(args, &params, source, reader),
        _ => unreachable!("filter_fastx_from_reader called with non-FASTX format"),
    }
}

/// When `--unordered` is set, eagerly load the full assignments file into a
/// `HashMap` so later lookups are order-independent. When omitted, return
/// `None` and the caller will stream with a lookahead buffer.
fn load_kraken_map_if_unordered(
    args: &FilterArgs,
    classifications_path: &Path,
) -> Result<Option<AHashMap<String, u32>>> {
    if args.unordered {
        info!(
            "Loading Kraken assignments (unordered) from: {}",
            classifications_path.display()
        );
        let m = KrakenResult::load_as_map(classifications_path)?;
        info!("Loaded {} Kraken assignments.", m.len());
        Ok(Some(m))
    } else {
        info!(
            "Streaming Kraken assignments from: {}",
            classifications_path.display()
        );
        Ok(None)
    }
}

fn make_fastx_source<'a>(
    map: Option<&'a AHashMap<String, u32>>,
    path: &'a Path,
) -> FastxSource<'a> {
    map.map_or(FastxSource::Stream(path), FastxSource::Map)
}

/// Bail if tree-dependent flags are set but no tree is available.
fn validate_tree_requirements(tree: &Option<KrakenTaxonomyTree>, args: &FilterArgs) -> Result<()> {
    if tree.is_some() {
        return Ok(());
    }
    for (flag, set) in [
        ("--allow-ancestors", args.allow_ancestors),
        ("--include-descendants", args.include_descendants),
    ] {
        if set {
            anyhow::bail!(
                "{flag} requires a taxonomy tree; pass --kraken-report \
                 (for SAM/BAM/CRAM, embedding via `krak annotate` also works)"
            );
        }
    }
    Ok(())
}

/// Build a taxonomy tree from an already-parsed SAM header (embedded report)
/// or from a report file on disk, whichever is available.
fn tree_from_header_or_file(
    header: &noodles::sam::Header,
    args: &FilterArgs,
) -> Result<Option<KrakenTaxonomyTree>> {
    let embedded = crate::kraken_report_embed::entries_from_header(header)?;
    if let Some(entries) = embedded {
        return Ok(Some(KrakenTaxonomyTree::from_entries(&entries)?));
    }
    args.kraken_report
        .as_deref()
        .map(KrakenTaxonomyTree::from_file)
        .transpose()
}

struct FilterParams<'a> {
    taxon_ids: &'a AHashSet<u32>,
    tree: Option<&'a KrakenTaxonomyTree>,
    allow_ancestors: bool,
    rescue_max_edit_distance: Option<u32>,
    rescue_max_indels: Option<u32>,
    rescue_max_indel_length: Option<u32>,
    rescue_n_adjustment: Option<u32>,
    keep_unannotated: bool,
}

/// Source of taxon-id assignments for the FASTX filter loop.
///
/// Mirrors the `Source` enum used by `crate::annotate` for SAM/BAM/CRAM:
/// `Map` is `--unordered` (everything pre-loaded), `Stream` is the default
/// (advances the assignments file in lock-step with the FASTX input, with a
/// lookahead buffer that absorbs modest disorder such as multi-threaded
/// Kraken v1 work-unit interleaving).
enum FastxSource<'a> {
    /// All assignments pre-loaded into a HashMap (order-independent lookup).
    Map(&'a AHashMap<String, u32>),
    /// Path to a Kraken assignments file streamed in lock-step with the input.
    Stream(&'a Path),
}

/// Per-loop lookup state. Holds either a borrowed map (the unordered path)
/// or an owned streaming reader, dispatching `lookup` accordingly.
enum LookupState<'a> {
    Map(&'a AHashMap<String, u32>),
    Stream(StreamingLookup<BufReader<std::fs::File>>),
}

impl<'a> LookupState<'a> {
    /// Build a `LookupState` from a `FastxSource`, opening the assignments
    /// file in streaming mode.
    fn from_source(source: FastxSource<'a>) -> Result<Self> {
        match source {
            FastxSource::Map(m) => Ok(Self::Map(m)),
            FastxSource::Stream(path) => {
                let file = std::fs::File::open(path).with_context(|| {
                    format!("failed to open Kraken assignments file: {}", path.display())
                })?;
                Ok(Self::Stream(StreamingLookup::new(BufReader::new(file))))
            }
        }
    }

    /// Resolve `name` to a taxon id. Returns `Ok(None)` when the name is
    /// absent from the assignments (matching the existing FASTX semantics:
    /// missing reads count toward `missing_classifications` and are treated as
    /// taxon 0 by the caller).
    fn lookup(&mut self, name: &str) -> Result<Option<u32>> {
        match self {
            Self::Map(m) => Ok(m.get(name).copied()),
            Self::Stream(s) => s.lookup(name),
        }
    }

    /// Drain any kraken assignments that were never consumed against an input
    /// record. The map-based source has no notion of unconsumed entries, so it
    /// returns 0; the stream-based source counts buffered entries plus the
    /// remaining stream tail.
    fn count_unconsumed(&mut self) -> usize {
        match self {
            Self::Map(_) => 0,
            Self::Stream(s) => s.count_unconsumed(),
        }
    }
}

/// Check taxon-based classification only (no edit-distance rescue).
/// Returns `Some(decision)` if the taxon matches on-taxa or ancestor criteria,
/// `None` if neither condition is met.
fn classify_by_taxon(taxon_id: u32, params: &FilterParams<'_>) -> Option<RecordDecision> {
    if params.taxon_ids.contains(&taxon_id) {
        return Some(RecordDecision::OnTaxa);
    }
    if params.allow_ancestors {
        if let Some(tree) = params.tree {
            if tree.is_ancestor_of_any(taxon_id, params.taxon_ids) {
                return Some(RecordDecision::AncestorRescue);
            }
        }
    }
    None
}

/// Classify a single FASTX record using taxon ID alone (no edit rescue).
fn classify_fastx(taxon_id: u32, params: &FilterParams<'_>) -> RecordDecision {
    classify_by_taxon(taxon_id, params).unwrap_or(RecordDecision::Reject)
}

/// FASTX output sink. For `.gz` outputs the file is encoded as BGZF
/// (block-gzip) via `gzp::ParCompress` with the libdeflate backend; the
/// resulting file is multi-member gzip and is read transparently by every
/// standard gzip reader (`zcat`, `gunzip`, `flate2`, samtools, etc.). Plain
/// outputs are buffered file writes with no compression.
enum FastxSink {
    Plain(std::io::BufWriter<std::fs::File>),
    Bgzf(std::io::BufWriter<Box<dyn gzp::ZWriter<std::fs::File>>>),
}

impl std::io::Write for FastxSink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Bgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Bgzf(w) => w.flush(),
        }
    }
}

impl FastxSink {
    fn create(path: &Path, threads: usize, compression_level: u32) -> std::io::Result<Self> {
        let is_gz = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);
        let file = std::fs::File::create(path)?;
        if is_gz {
            Ok(Self::Bgzf(std::io::BufWriter::new(build_bgzf_writer(
                file,
                threads,
                compression_level,
            )?)))
        } else {
            Ok(Self::Plain(std::io::BufWriter::new(file)))
        }
    }

    /// Finalize the writer chain. For BGZF outputs this flushes the BufWriter
    /// and then calls `ZWriter::finish` so the BGZF EOF block is emitted.
    fn finalize(self) -> std::io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Bgzf(w) => {
                let mut zw = w.into_inner().map_err(|e| e.into_error())?;
                zw.finish().map_err(std::io::Error::other)?;
                Ok(())
            }
        }
    }
}

/// Build a libdeflate-backed BGZF writer with `threads` compression workers
/// and the given level (0-9). At `threads = 1`, one compressor and one writer
/// thread pipeline with the caller.
fn build_bgzf_writer(
    file: std::fs::File,
    threads: usize,
    level: u32,
) -> std::io::Result<Box<dyn gzp::ZWriter<std::fs::File>>> {
    let parz = gzp::par::compress::ParCompressBuilder::<gzp::deflate::Bgzf>::new()
        .num_threads(threads.max(1))
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?
        .compression_level(gzp::Compression::new(level))
        .from_writer(file);
    Ok(Box::new(parz))
}

/// Strip the FASTQ description suffix (after first whitespace), then strip any
/// trailing `/1` or `/2` Kraken-style pair suffix.  Returns the base read name.
fn fastq_base_name(name_bytes: &[u8]) -> Result<String> {
    let raw = std::str::from_utf8(name_bytes).context("non-UTF-8 FASTQ read name")?;
    let no_desc = raw.split_whitespace().next().unwrap_or(raw);
    Ok(crate::strip_pair_suffix(no_desc).to_owned())
}

fn filter_fastq_path(
    args: &FilterArgs,
    params: &FilterParams<'_>,
    source: FastxSource<'_>,
) -> Result<TaxaFilterMetric> {
    let reader_in = crate::open_fastx_reader(&args.input)
        .with_context(|| format!("failed to open FASTQ: {}", args.input.display()))?;
    filter_fastq_from_reader(args, params, source, reader_in)
}

fn filter_fastq_from_reader(
    args: &FilterArgs,
    params: &FilterParams<'_>,
    source: FastxSource<'_>,
    reader_in: Box<dyn std::io::BufRead>,
) -> Result<TaxaFilterMetric> {
    use noodles::fastq;

    let mut lookup = LookupState::from_source(source)?;

    let mut reader = fastq::io::Reader::new(reader_in);

    let out_sink = FastxSink::create(&args.output, args.threads, args.compression_level)
        .with_context(|| format!("failed to create output FASTQ: {}", args.output.display()))?;
    let mut out = fastq::io::Writer::new(out_sink);

    let mut reject_out: Option<fastq::io::Writer<FastxSink>> = args
        .rejects
        .as_deref()
        .map(|p| {
            FastxSink::create(p, args.threads, args.compression_level)
                .with_context(|| format!("failed to create rejects FASTQ: {}", p.display()))
                .map(fastq::io::Writer::new)
        })
        .transpose()?;

    let mut metric = TaxaFilterMetric {
        missing_classifications: Some(0),
        ..Default::default()
    };

    let mut first = fastq::Record::default();
    if reader.read_record(&mut first)? == 0 {
        metric.finalize();
        return Ok(metric);
    }
    let first_raw_name = {
        let raw = std::str::from_utf8(first.name()).context("non-UTF-8 FASTQ read name")?;
        raw.split_whitespace().next().unwrap_or(raw).to_owned()
    };
    let interleaved = !args.per_record && first_raw_name.ends_with("/1");

    // Validate: if first record ends with /2 (wrong order), and NOT per-record, raise error.
    if !args.per_record && first_raw_name.ends_with("/2") {
        anyhow::bail!(
            "interleaved FASTQ: first record '{first_raw_name}' ends with /2; \
             interleaved input must begin with a /1 record.\n\
             Use --per-record to process records independently."
        );
    }

    // `pending` holds the already-read first record; subsequent iterations read
    // fresh records from `reader`.
    let mut pending: Option<fastq::Record> = Some(first);

    loop {
        let r1 = if let Some(r) = pending.take() {
            r
        } else {
            let mut r = fastq::Record::default();
            match reader.read_record(&mut r)? {
                0 => break,
                _ => r,
            }
        };

        let r2: Option<fastq::Record> = if interleaved {
            let mut r = fastq::Record::default();
            if reader.read_record(&mut r)? == 0 {
                let base = fastq_base_name(r1.name())?;
                anyhow::bail!(
                    "interleaved FASTQ: '{base}/1' has no following /2 record (truncated file?).\n\
                     Use --per-record to process records independently."
                );
            }
            // Validate pairing
            let r1_base = fastq_base_name(r1.name())?;
            let r2_base = fastq_base_name(r.name())?;
            if r1_base != r2_base {
                anyhow::bail!(
                    "interleaved FASTQ: '{r1_base}/1' is not immediately followed by '{r1_base}/2' \
                     (got '{r2_base}').\n\
                     Use --per-record to process records independently."
                );
            }
            Some(r)
        } else {
            None
        };

        let rec_count = 1 + r2.is_some() as u64;
        metric.total += rec_count;

        let name = fastq_base_name(r1.name())?;
        let taxon_id = match lookup.lookup(&name)? {
            Some(id) => id,
            None => {
                if let Some(ref mut v) = metric.missing_classifications {
                    *v += 1;
                }
                0
            }
        };

        match classify_fastx(taxon_id, params) {
            RecordDecision::OnTaxa => {
                metric.on_taxa += rec_count;
                out.write_record(&r1).context("failed to write FASTQ")?;
                if let Some(ref r) = r2 {
                    out.write_record(r).context("failed to write FASTQ")?;
                }
            }
            RecordDecision::AncestorRescue => {
                metric.rescued_ancestors += rec_count;
                out.write_record(&r1).context("failed to write FASTQ")?;
                if let Some(ref r) = r2 {
                    out.write_record(r).context("failed to write FASTQ")?;
                }
            }
            _ => {
                if let Some(ref mut rout) = reject_out {
                    rout.write_record(&r1)
                        .context("failed to write reject FASTQ")?;
                    if let Some(ref r) = r2 {
                        rout.write_record(r)
                            .context("failed to write reject FASTQ")?;
                    }
                }
            }
        }
    }

    out.into_inner()
        .finalize()
        .context("failed to finalize FASTQ output")?;
    if let Some(rout) = reject_out {
        rout.into_inner()
            .finalize()
            .context("failed to finalize reject FASTQ output")?;
    }
    let leftover = lookup.count_unconsumed();
    if leftover > 0 {
        warn!(
            "{leftover} Kraken assignments remained after the input ran out; \
             the --classifications file appears to extend past the FASTQ input."
        );
    }
    metric.finalize();
    Ok(metric)
}

fn filter_fasta_path(
    args: &FilterArgs,
    params: &FilterParams<'_>,
    source: FastxSource<'_>,
) -> Result<TaxaFilterMetric> {
    let reader_in = crate::open_fastx_reader(&args.input)
        .with_context(|| format!("failed to open FASTA: {}", args.input.display()))?;
    filter_fasta_from_reader(args, params, source, reader_in)
}

fn filter_fasta_from_reader(
    args: &FilterArgs,
    params: &FilterParams<'_>,
    source: FastxSource<'_>,
    reader_in: Box<dyn std::io::BufRead>,
) -> Result<TaxaFilterMetric> {
    use noodles::fasta;

    let mut lookup = LookupState::from_source(source)?;
    let mut reader = fasta::io::Reader::new(reader_in);

    let out_sink = FastxSink::create(&args.output, args.threads, args.compression_level)
        .with_context(|| format!("failed to create output FASTA: {}", args.output.display()))?;
    let mut out = fasta::io::Writer::new(out_sink);

    let mut reject_out: Option<fasta::io::Writer<FastxSink>> = args
        .rejects
        .as_deref()
        .map(|p| {
            FastxSink::create(p, args.threads, args.compression_level)
                .with_context(|| format!("failed to create rejects FASTA: {}", p.display()))
                .map(fasta::io::Writer::new)
        })
        .transpose()?;

    let mut metric = TaxaFilterMetric {
        missing_classifications: Some(0),
        ..Default::default()
    };

    for result in reader.records() {
        let record = result.context("failed to read FASTA record")?;
        let name = std::str::from_utf8(record.name()).context("non-UTF-8 FASTA record name")?;
        let name = crate::strip_pair_suffix(name);

        metric.total += 1;
        let taxon_id_opt = lookup.lookup(name)?;
        if taxon_id_opt.is_none() {
            if let Some(ref mut v) = metric.missing_classifications {
                *v += 1;
            }
        }
        let taxon_id = taxon_id_opt.unwrap_or(0);

        match classify_fastx(taxon_id, params) {
            RecordDecision::OnTaxa => {
                metric.on_taxa += 1;
                out.write_record(&record)
                    .context("failed to write FASTA record")?;
            }
            RecordDecision::AncestorRescue => {
                metric.rescued_ancestors += 1;
                out.write_record(&record)
                    .context("failed to write FASTA record")?;
            }
            _ => {
                if let Some(ref mut rout) = reject_out {
                    rout.write_record(&record)
                        .context("failed to write reject FASTA record")?;
                }
            }
        }
    }

    out.into_inner()
        .finalize()
        .context("failed to finalize FASTA output")?;
    if let Some(rout) = reject_out {
        rout.into_inner()
            .finalize()
            .context("failed to finalize reject FASTA output")?;
    }
    let leftover = lookup.count_unconsumed();
    if leftover > 0 {
        warn!(
            "{leftover} Kraken assignments remained after the input ran out; \
             the --classifications file appears to extend past the FASTA input."
        );
    }
    metric.finalize();
    Ok(metric)
}

fn filter_bam(args: &FilterArgs) -> Result<Vec<TaxaFilterMetric>> {
    let reader = crate::open_bam_reader(&args.input)?;
    filter_bam_with_reader(args, reader)
}

fn filter_bam_with_reader<R: std::io::Read>(
    args: &FilterArgs,
    mut reader: noodles::bam::io::Reader<R>,
) -> Result<Vec<TaxaFilterMetric>> {
    use noodles::bam;

    let header = reader.read_header().context("failed to read BAM header")?;
    let (tree, expanded) = prepare_filter_state(&header, args)?;
    let params = make_params(&expanded, &tree, args);

    let mut writer = bam::io::writer::Builder
        .build_from_path(&args.output)
        .with_context(|| format!("failed to create output BAM: {}", args.output.display()))?;
    writer
        .write_header(&header)
        .context("failed to write header")?;

    let mut reject_writer: Option<bam::io::Writer<_>> = args
        .rejects
        .as_deref()
        .map(|p| {
            bam::io::writer::Builder
                .build_from_path(p)
                .with_context(|| format!("failed to create rejects BAM: {}", p.display()))
                .and_then(|mut w| {
                    w.write_header(&header)?;
                    Ok(w)
                })
        })
        .transpose()?;

    let metrics = dispatch_filter_records(
        args,
        &params,
        &header,
        reader.record_bufs(&header),
        &mut writer,
        &mut reject_writer,
        "samtools sort -n input.bam",
    )?;

    writer
        .into_inner()
        .finish()
        .context("failed to finish BAM BGZF stream")?;
    if let Some(rw) = reject_writer {
        rw.into_inner()
            .finish()
            .context("failed to finish rejects BAM BGZF stream")?;
    }
    crate::maybe_index_alignment_output(&args.output, &header, AlignmentFormat::Bam)?;
    if let Some(rejects_path) = args.rejects.as_deref() {
        crate::maybe_index_alignment_output(rejects_path, &header, AlignmentFormat::Bam)?;
    }
    Ok(metrics)
}

fn filter_sam(args: &FilterArgs) -> Result<Vec<TaxaFilterMetric>> {
    use noodles::sam;
    let file = std::fs::File::open(&args.input)
        .with_context(|| format!("failed to open SAM: {}", args.input.display()))?;
    let reader = sam::io::Reader::new(BufReader::new(file));
    filter_sam_with_reader(args, reader)
}

fn filter_sam_with_reader<R: std::io::BufRead>(
    args: &FilterArgs,
    mut reader: noodles::sam::io::Reader<R>,
) -> Result<Vec<TaxaFilterMetric>> {
    use noodles::sam;

    let header = reader.read_header().context("failed to read SAM header")?;
    let (tree, expanded) = prepare_filter_state(&header, args)?;
    let params = make_params(&expanded, &tree, args);

    let out_file = std::fs::File::create(&args.output)
        .with_context(|| format!("failed to create output SAM: {}", args.output.display()))?;
    let mut writer = sam::io::Writer::new(std::io::BufWriter::new(out_file));
    writer
        .write_header(&header)
        .context("failed to write header")?;

    let mut reject_writer: Option<sam::io::Writer<_>> = args
        .rejects
        .as_deref()
        .map(|p| {
            std::fs::File::create(p)
                .with_context(|| format!("failed to create rejects SAM: {}", p.display()))
                .and_then(|f| {
                    let mut w = sam::io::Writer::new(std::io::BufWriter::new(f));
                    w.write_header(&header)
                        .context("failed to write rejects header")?;
                    Ok(w)
                })
        })
        .transpose()?;

    let metrics = dispatch_filter_records(
        args,
        &params,
        &header,
        reader.record_bufs(&header),
        &mut writer,
        &mut reject_writer,
        "samtools sort -n input.sam",
    )?;

    writer
        .into_inner()
        .flush()
        .context("failed to flush SAM writer")?;
    if let Some(rw) = reject_writer {
        rw.into_inner()
            .flush()
            .context("failed to flush rejects SAM writer")?;
    }
    Ok(metrics)
}

fn filter_cram(args: &FilterArgs) -> Result<Vec<TaxaFilterMetric>> {
    let reader = crate::open_cram_reader(&args.input, args.cram_reference.as_deref())?;
    filter_cram_with_reader(args, reader)
}

fn filter_cram_with_reader<R: std::io::Read>(
    args: &FilterArgs,
    mut reader: noodles::cram::io::Reader<R>,
) -> Result<Vec<TaxaFilterMetric>> {
    let header = reader.read_header().context("failed to read CRAM header")?;
    crate::require_cram_reference_if_mapped(&header, args.cram_reference.as_deref())?;
    let (tree, expanded) = prepare_filter_state(&header, args)?;
    let params = make_params(&expanded, &tree, args);

    let mut writer = crate::open_cram_writer(&args.output, args.cram_reference.as_deref())?;
    writer
        .write_header(&header)
        .context("failed to write CRAM header")?;

    let mut reject_writer = args
        .rejects
        .as_deref()
        .map(|p| {
            crate::open_cram_writer(p, args.cram_reference.as_deref()).and_then(|mut w| {
                w.write_header(&header)
                    .context("failed to write rejects CRAM header")?;
                Ok(w)
            })
        })
        .transpose()?;

    let result = dispatch_filter_records(
        args,
        &params,
        &header,
        reader.records(&header),
        &mut writer,
        &mut reject_writer,
        "samtools sort -n -T ref.fa input.cram",
    );

    writer
        .try_finish(&header)
        .context("failed to finish CRAM writer")?;
    if let Some(rw) = reject_writer.as_mut() {
        rw.try_finish(&header)
            .context("failed to finish CRAM rejects writer")?;
    }
    crate::maybe_index_alignment_output(&args.output, &header, AlignmentFormat::Cram)?;
    if let Some(rejects_path) = args.rejects.as_deref() {
        crate::maybe_index_alignment_output(rejects_path, &header, AlignmentFormat::Cram)?;
    }
    result
}

fn make_params<'a>(
    expanded: &'a AHashSet<u32>,
    tree: &'a Option<KrakenTaxonomyTree>,
    args: &FilterArgs,
) -> FilterParams<'a> {
    FilterParams {
        taxon_ids: expanded,
        tree: tree.as_ref(),
        allow_ancestors: args.allow_ancestors,
        rescue_max_edit_distance: args.rescue_max_edit_distance,
        rescue_max_indels: args.rescue_max_indels,
        rescue_max_indel_length: args.rescue_max_indel_length,
        rescue_n_adjustment: args.rescue_n_adjustment,
        keep_unannotated: args.keep_unannotated,
    }
}

/// Build the taxonomy tree (from header-embedded report or `--kraken-report`),
/// validate it, and expand the user's target taxon set with descendants and/or
/// taxon 0. Shared across the BAM/SAM/CRAM dispatch entrypoints.
fn prepare_filter_state(
    header: &noodles::sam::Header,
    args: &FilterArgs,
) -> Result<(Option<KrakenTaxonomyTree>, AHashSet<u32>)> {
    let tree = tree_from_header_or_file(header, args)?;
    validate_tree_requirements(&tree, args)?;
    let expanded = build_taxon_set(
        &args.taxon_ids,
        tree.as_ref(),
        args.include_descendants,
        args.include_unclassified,
    )?;
    Ok((tree, expanded))
}

/// Build the keep/reject writer closures and run either per-record or
/// per-template filtering against `records`. `sort_hint` is the
/// format-specific suggestion that appears in the queryname-sort error
/// (e.g. `samtools sort -n input.bam`).
fn dispatch_filter_records<I, W, Wrej>(
    args: &FilterArgs,
    params: &FilterParams<'_>,
    header: &noodles::sam::Header,
    records: I,
    writer: &mut W,
    reject_writer: &mut Option<Wrej>,
    sort_hint: &str,
) -> Result<Vec<TaxaFilterMetric>>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    W: AlignmentWrite,
    Wrej: AlignmentWrite,
{
    use noodles::sam::alignment::record_buf::RecordBuf;
    let keep = |r: &RecordBuf| {
        writer
            .write_alignment_record(header, r)
            .context("write failed")
    };
    let rej = |r: &RecordBuf| {
        reject_writer.as_mut().map_or(Ok(()), |rw| {
            rw.write_alignment_record(header, r)
                .context("write reject failed")
        })
    };

    if args.per_record {
        Ok(vec![filter_by_record(records, params, keep, rej)?])
    } else if crate::is_query_grouped(header) {
        let (rec, tmpl) = filter_by_template(records, params, keep, rej)?;
        Ok(vec![rec, tmpl])
    } else {
        anyhow::bail!(
            "input is not queryname-sorted or query-grouped; \
             krak filter evaluates read pairs as a unit by default.\n\
             Sort with:  {sort_hint}\n\
             Override:   --per-record  (processes each alignment line independently)"
        )
    }
}

/// Stream-filter records one at a time; no buffering; O(1) memory per record.
fn filter_by_record<I, W, R>(
    records: I,
    params: &FilterParams<'_>,
    mut write_keep: W,
    mut write_reject: R,
) -> Result<TaxaFilterMetric>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    W: FnMut(&noodles::sam::alignment::record_buf::RecordBuf) -> Result<()>,
    R: FnMut(&noodles::sam::alignment::record_buf::RecordBuf) -> Result<()>,
{
    let mut metric = TaxaFilterMetric {
        rescued_variants: Some(0),
        unannotated_but_kept: Some(0),
        missing_ti_tag: Some(0),
        missing_classifications: None,
        rescue_max_edit_distance: params.rescue_max_edit_distance,
        rescue_max_indels: params.rescue_max_indels,
        rescue_max_indel_length: params.rescue_max_indel_length,
        rescue_n_adjustment: params.rescue_n_adjustment,
        keep_unannotated: params.keep_unannotated,
        ..Default::default()
    };

    for result in records {
        let record = result.context("failed to read alignment record")?;
        metric.total += 1;
        let taxon_id = read_ti_tag(&record);

        if taxon_id.is_none() {
            if let Some(ref mut v) = metric.missing_ti_tag {
                *v += 1;
            }
        }

        match classify_record(&record, taxon_id, params) {
            RecordDecision::OnTaxa => {
                metric.on_taxa += 1;
                write_keep(&record)?;
            }
            RecordDecision::AncestorRescue => {
                metric.rescued_ancestors += 1;
                write_keep(&record)?;
            }
            RecordDecision::EditRescue => {
                if let Some(ref mut v) = metric.rescued_variants {
                    *v += 1;
                }
                write_keep(&record)?;
            }
            RecordDecision::UnannotatedKeep => {
                if let Some(ref mut v) = metric.unannotated_but_kept {
                    *v += 1;
                }
                write_keep(&record)?;
            }
            RecordDecision::Reject => {
                write_reject(&record)?;
            }
        }
    }

    metric.finalize();
    Ok(metric)
}

/// Filter by template (stream consecutive records grouped by query name).
fn filter_by_template<I, W, R>(
    records: I,
    params: &FilterParams<'_>,
    mut write_keep: W,
    mut write_reject: R,
) -> Result<(TaxaFilterMetric, TaxaFilterMetric)>
where
    I: Iterator<Item = std::io::Result<noodles::sam::alignment::record_buf::RecordBuf>>,
    W: FnMut(&noodles::sam::alignment::record_buf::RecordBuf) -> Result<()>,
    R: FnMut(&noodles::sam::alignment::record_buf::RecordBuf) -> Result<()>,
{
    use noodles::sam::alignment::record_buf::RecordBuf;

    let mut rec_metric = TaxaFilterMetric {
        rescued_variants: Some(0),
        unannotated_but_kept: Some(0),
        missing_ti_tag: Some(0),
        missing_classifications: None,
        rescue_max_edit_distance: params.rescue_max_edit_distance,
        rescue_max_indels: params.rescue_max_indels,
        rescue_max_indel_length: params.rescue_max_indel_length,
        rescue_n_adjustment: params.rescue_n_adjustment,
        keep_unannotated: params.keep_unannotated,
        ..Default::default()
    };
    let mut tmpl_metric = TaxaFilterMetric {
        template: true,
        rescued_variants: Some(0),
        unannotated_but_kept: Some(0),
        missing_ti_tag: Some(0),
        missing_classifications: None,
        rescue_max_edit_distance: params.rescue_max_edit_distance,
        rescue_max_indels: params.rescue_max_indels,
        rescue_max_indel_length: params.rescue_max_indel_length,
        rescue_n_adjustment: params.rescue_n_adjustment,
        keep_unannotated: params.keep_unannotated,
        ..Default::default()
    };

    let mut current_name: Option<String> = None;
    let mut current_template: Vec<RecordBuf> = Vec::new();

    // Helper: process one complete template and update metrics.
    let process_template = |template: &[RecordBuf],
                            rec_m: &mut TaxaFilterMetric,
                            tmpl_m: &mut TaxaFilterMetric,
                            write_keep: &mut W,
                            write_reject: &mut R|
     -> Result<()> {
        tmpl_m.total += 1;
        rec_m.total += template.len() as u64;

        let taxon_ids: Vec<Option<u32>> = template.iter().map(read_ti_tag).collect();
        let missing = taxon_ids.iter().filter(|t| t.is_none()).count() as u64;
        if missing > 0 {
            if let Some(ref mut v) = tmpl_m.missing_ti_tag {
                *v += 1;
            }
        }
        if let Some(ref mut v) = rec_m.missing_ti_tag {
            *v += missing;
        }

        let template_refs: Vec<&RecordBuf> = template.iter().collect();
        match classify_template(&template_refs, &taxon_ids, params) {
            RecordDecision::OnTaxa => {
                tmpl_m.on_taxa += 1;
                rec_m.on_taxa += template.len() as u64;
                for r in template {
                    write_keep(r)?;
                }
            }
            RecordDecision::AncestorRescue => {
                tmpl_m.rescued_ancestors += 1;
                rec_m.rescued_ancestors += template.len() as u64;
                for r in template {
                    write_keep(r)?;
                }
            }
            RecordDecision::EditRescue => {
                if let Some(ref mut v) = tmpl_m.rescued_variants {
                    *v += 1;
                }
                if let Some(ref mut v) = rec_m.rescued_variants {
                    *v += template.len() as u64;
                }
                for r in template {
                    write_keep(r)?;
                }
            }
            RecordDecision::UnannotatedKeep => {
                if let Some(ref mut v) = tmpl_m.unannotated_but_kept {
                    *v += 1;
                }
                if let Some(ref mut v) = rec_m.unannotated_but_kept {
                    *v += template.len() as u64;
                }
                for r in template {
                    write_keep(r)?;
                }
            }
            RecordDecision::Reject => {
                for r in template {
                    write_reject(r)?;
                }
            }
        }
        Ok(())
    };

    for result in records {
        let record = result.context("failed to read alignment record")?;
        let name_bytes = record.name().ok_or_else(|| {
            anyhow::anyhow!(
                "record has no name; template mode requires all records to have a QNAME"
            )
        })?;
        let name = std::str::from_utf8(name_bytes)
            .context("non-UTF-8 record name")?
            .to_owned();

        if Some(name.as_str()) != current_name.as_deref() {
            // Flush the previous template before starting the next one.
            if !current_template.is_empty() {
                process_template(
                    &current_template,
                    &mut rec_metric,
                    &mut tmpl_metric,
                    &mut write_keep,
                    &mut write_reject,
                )?;
                current_template.clear();
            }
            current_name = Some(name);
        }
        current_template.push(record);
    }
    // Flush the last template.
    if !current_template.is_empty() {
        process_template(
            &current_template,
            &mut rec_metric,
            &mut tmpl_metric,
            &mut write_keep,
            &mut write_reject,
        )?;
    }

    rec_metric.finalize();
    tmpl_metric.finalize();
    Ok((rec_metric, tmpl_metric))
}

#[derive(Debug, PartialEq, Eq)]
enum RecordDecision {
    OnTaxa,
    AncestorRescue,
    EditRescue,
    UnannotatedKeep,
    Reject,
}

fn classify_record(
    record: &noodles::sam::alignment::record_buf::RecordBuf,
    taxon_id: Option<u32>,
    params: &FilterParams<'_>,
) -> RecordDecision {
    let Some(tid) = taxon_id else {
        if params.keep_unannotated {
            return RecordDecision::UnannotatedKeep;
        }
        return RecordDecision::Reject;
    };
    // Has ti tag: apply cascade.
    if let Some(d) = classify_by_taxon(tid, params) {
        return d;
    }
    if within_edit_limits(record, params) {
        return RecordDecision::EditRescue;
    }
    RecordDecision::Reject
}

fn classify_template(
    template: &[&noodles::sam::alignment::record_buf::RecordBuf],
    taxon_ids: &[Option<u32>],
    params: &FilterParams<'_>,
) -> RecordDecision {
    // Taxon decisions are made on PRIMARY records only; supplementary and
    // secondary alignments share their primary's QNAME and are typically
    // re-annotated independently, so their `ti` tag may disagree. Letting a
    // supplementary on-taxa taxon rescue an off-taxa primary would silently
    // keep records the user did not intend to retain.
    let primary_taxon_ids: Vec<Option<u32>> = template
        .iter()
        .zip(taxon_ids.iter())
        .filter_map(|(r, tid)| {
            (!r.flags().is_supplementary() && !r.flags().is_secondary()).then_some(*tid)
        })
        .collect();

    // On-taxa: any primary read with a taxon in the target set.
    if primary_taxon_ids
        .iter()
        .flatten()
        .any(|&tid| params.taxon_ids.contains(&tid))
    {
        return RecordDecision::OnTaxa;
    }

    // Ancestor rescue: any primary read whose taxon is an ancestor of a target.
    if primary_taxon_ids
        .iter()
        .flatten()
        .filter_map(|&tid| classify_by_taxon(tid, params))
        .any(|d| matches!(d, RecordDecision::AncestorRescue))
    {
        return RecordDecision::AncestorRescue;
    }

    // Edit rescue: ALL primary reads must pass. Only applies when at least one
    // primary has a taxon annotation; a fully-unannotated template must not be
    // rescued here (consistent with classify_record, which never edit-rescues
    // an unannotated record).
    let has_any_taxon = primary_taxon_ids.iter().any(|t| t.is_some());
    let primary: Vec<_> = template
        .iter()
        .filter(|r| !r.flags().is_supplementary() && !r.flags().is_secondary())
        .collect();

    if has_any_taxon && !primary.is_empty() && primary.iter().all(|r| within_edit_limits(r, params))
    {
        return RecordDecision::EditRescue;
    }

    // Unannotated: only applies when ALL primary reads lack a taxon
    // annotation. A partially-annotated template (some primaries have a ti
    // tag, even if off-taxa) is not unannotated; it is annotated but failed
    // all classification criteria, so it should be Rejected.
    if primary_taxon_ids.iter().all(|t| t.is_none()) && params.keep_unannotated {
        return RecordDecision::UnannotatedKeep;
    }

    RecordDecision::Reject
}

/// Read the `ti` tag from a record, returning `None` if absent or not an integer.
pub(crate) fn read_ti_tag(record: &noodles::sam::alignment::record_buf::RecordBuf) -> Option<u32> {
    use noodles::sam::alignment::record_buf::data::field::Value;
    let n: i64 = match record.data().get(&crate::TI_TAG)? {
        Value::Int8(n) => i64::from(*n),
        Value::UInt8(n) => i64::from(*n),
        Value::Int16(n) => i64::from(*n),
        Value::UInt16(n) => i64::from(*n),
        Value::Int32(n) => i64::from(*n),
        Value::UInt32(n) => i64::from(*n),
        _ => return None,
    };
    u32::try_from(n).ok()
}

/// Returns `true` if the record is within all edit distance limits.
///
/// Returns `false` if the record has no MD tag or is unmapped.
fn within_edit_limits(
    record: &noodles::sam::alignment::record_buf::RecordBuf,
    params: &FilterParams<'_>,
) -> bool {
    // Rescue is opt-in: if no rescue knob is set, edit-rescue is disabled
    // entirely (even a 0-edit record is not rescued).
    if params.rescue_max_edit_distance.is_none()
        && params.rescue_max_indels.is_none()
        && params.rescue_max_indel_length.is_none()
        && params.rescue_n_adjustment.is_none()
    {
        return false;
    }
    if record.flags().is_unmapped() {
        return false;
    }
    let Some(edits) = count_edits(record) else {
        return false; // no MD tag
    };

    // An unset --rescue-max-edit-distance means "no edit-count ceiling"; only
    // the indel/length knobs (if set) gate the record. A `Some(n)` value caps
    // total edits at `n`, optionally reduced by the N-adjustment.
    if let Some(base) = params.rescue_max_edit_distance {
        let max_edits = match params.rescue_n_adjustment {
            Some(n) if n > 0 => {
                let n_count = count_aligned_ns(record) as u32;
                let adjustment = n_count / n;
                base.saturating_sub(adjustment)
            }
            _ => base,
        };
        if edits.total() > max_edits as i64 {
            return false;
        }
    }
    if let Some(max_ind) = params.rescue_max_indels {
        if edits.indel_count() > max_ind as usize {
            return false;
        }
    }
    if let Some(max_len) = params.rescue_max_indel_length {
        if edits.max_indel_length() > max_len {
            return false;
        }
    }
    true
}

/// Count N bases that fall at aligned (M/=/X) CIGAR positions.
///
/// Ns at soft-clip (S) or insertion (I) positions are excluded because the
/// MD tag only records mismatches at reference-consuming aligned positions;
/// subtracting those Ns would incorrectly reduce the substitution count.
fn count_aligned_ns(record: &noodles::sam::alignment::record_buf::RecordBuf) -> i64 {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let seq = record.sequence().as_ref();
    let mut read_pos: usize = 0;
    let mut count: i64 = 0;
    for op in record.cigar().as_ref().iter() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in read_pos..read_pos + len {
                    if seq
                        .get(i)
                        .map(|&b| b.eq_ignore_ascii_case(&b'N'))
                        .unwrap_or(false)
                    {
                        count += 1;
                    }
                }
                read_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => read_pos += len,
            _ => {} // Deletion, Skip, HardClip, Pad do not consume read bases
        }
    }
    count
}

/// Count edits from the reference using MD tag (substitutions) and CIGAR (indels).
///
/// N bases at aligned positions are only subtracted from the substitution
/// count when the MD tag actually recorded a mismatch at that position (the
/// read had an N where the reference had A/C/G/T). When the *reference* base
/// itself is N there is no MD entry; all aligned bases match the N reference;
/// so subtracting the read's Ns would over-correct. Walking MD and CIGAR
/// jointly avoids that pitfall.
fn count_edits(record: &noodles::sam::alignment::record_buf::RecordBuf) -> Option<EditCounts> {
    use noodles::sam::alignment::record::cigar::{op::Kind, Op};
    use noodles::sam::alignment::record_buf::data::field::Value;

    let md_bytes: &[u8] = match record.data().get(&MD_TAG)? {
        Value::String(s) => s.as_ref(),
        _ => return None,
    };

    // Walk CIGAR + MD jointly. CIGAR drives the read-position cursor; MD drives
    // the reference-aware substitution decisions. For each MD substitution we
    // check the read base at the corresponding aligned position: if it is N we
    // skip the substitution (the read N is "free" against an A/C/G/T reference);
    // otherwise it counts. Reference-is-N positions produce no MD substitution
    // entry, so they're naturally absent from the count.
    let seq_bytes = record.sequence().as_ref();
    let mut indels: Vec<i32> = Vec::new();
    let mut substitutions: i64 = 0;

    // CIGAR iteration state.
    let mut cigar_iter = record.cigar().as_ref().iter();
    let mut current_op: Option<(Kind, usize)> = cigar_iter.next().map(|o| (o.kind(), o.len()));
    let mut read_pos: usize = 0;

    // Advance through CIGAR until we land on a read+ref-consuming aligned base
    // (M/=/X). Returns Some(read_pos) for an aligned base, or None at EOF of
    // the CIGAR. Soft-clips and insertions advance read_pos but don't surface
    // an aligned base; deletions/skips/etc don't advance read_pos.
    let consume_aligned_base = |cigar_iter: &mut std::slice::Iter<Op>,
                                current_op: &mut Option<(Kind, usize)>,
                                read_pos: &mut usize|
     -> Option<usize> {
        loop {
            match current_op {
                None => return None,
                Some((kind, remaining)) => {
                    if *remaining == 0 {
                        *current_op = cigar_iter.next().map(|o| (o.kind(), o.len()));
                        continue;
                    }
                    match *kind {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            let pos = *read_pos;
                            *read_pos += 1;
                            *remaining -= 1;
                            return Some(pos);
                        }
                        Kind::Insertion | Kind::SoftClip => {
                            *read_pos += *remaining;
                            *remaining = 0;
                        }
                        _ => {
                            *remaining = 0;
                        }
                    }
                }
            }
        }
    };

    // Walk MD bytes.
    let mut i = 0;
    while i < md_bytes.len() {
        let b = md_bytes[i];
        match b {
            b'0'..=b'9' => {
                let mut n: usize = 0;
                while i < md_bytes.len() && md_bytes[i].is_ascii_digit() {
                    n = n * 10 + (md_bytes[i] - b'0') as usize;
                    i += 1;
                }
                // Advance the CIGAR cursor over `n` aligned bases.
                for _ in 0..n {
                    if consume_aligned_base(&mut cigar_iter, &mut current_op, &mut read_pos)
                        .is_none()
                    {
                        break;
                    }
                }
            }
            b'^' => {
                // Deletion run: skip every alphabetic byte that follows.
                i += 1;
                while i < md_bytes.len() && md_bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
            }
            b'A'..=b'Z' | b'a'..=b'z' => {
                // Substitution at the next aligned read position.
                if let Some(pos) =
                    consume_aligned_base(&mut cigar_iter, &mut current_op, &mut read_pos)
                {
                    let read_base = seq_bytes.get(pos).copied().unwrap_or(b'N');
                    if !read_base.eq_ignore_ascii_case(&b'N') {
                        substitutions += 1;
                    }
                } else {
                    // Defensive: MD called for a substitution but no aligned
                    // base remains. Count it (matches prior behavior of
                    // counting MD-letters even without sequence support).
                    substitutions += 1;
                }
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }

    // Collect indel events from CIGAR.
    for op in record.cigar().as_ref().iter() {
        match op.kind() {
            Kind::Insertion => indels.push(op.len() as i32),
            Kind::Deletion => indels.push(-(op.len() as i32)),
            _ => {}
        }
    }

    Some(EditCounts {
        substitutions,
        indels,
    })
}

fn write_metrics(path: &Path, metrics: &[TaxaFilterMetric]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to create metrics file: {}", path.display()))?;
    for m in metrics {
        wtr.serialize(m)
            .context("failed to serialize metrics row")?;
    }
    wtr.flush().context("failed to flush metrics file")?;
    Ok(())
}

fn log_metrics(metrics: &[TaxaFilterMetric]) {
    for m in metrics {
        let unit = if m.template { "templates" } else { "records" };
        info!("On-taxa {unit}:          {}", m.on_taxa);
        info!("Ancestor-rescued {unit}: {}", m.rescued_ancestors);
        if let Some(rv) = m.rescued_variants {
            info!("Edit-rescued {unit}:     {rv}");
        }
        if let Some(ubk) = m.unannotated_but_kept {
            info!("Unannotated-kept {unit}: {ubk}");
        }
        info!("Kept {unit}:             {}", m.num_kept);
        info!("Rejected {unit}:         {}", m.num_filtered);
        info!("Fraction {unit} removed: {:.4}", m.frac_removed);
        if let Some(mtt) = m.missing_ti_tag {
            if mtt > 0 && !m.keep_unannotated {
                warn!(
                    "{mtt} {unit} missing `ti` tag and were not annotated; \
                     run `krak annotate` first."
                );
            }
        }
        if let Some(mko) = m.missing_classifications {
            if mko > 0 {
                warn!("{mko} {unit} not found in --classifications.");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::{op::Kind, Op};
    use noodles::sam::alignment::record_buf::{data::field::Value, Cigar, RecordBuf, Sequence};

    fn record_with_ti(taxon_id: u32) -> RecordBuf {
        let mut r = RecordBuf::default();
        r.data_mut()
            .insert(crate::TI_TAG.into(), Value::UInt32(taxon_id));
        r
    }

    fn record_with_ti_and_md(taxon_id: u32, md: &str, ops: Vec<Op>, seq: &[u8]) -> RecordBuf {
        use noodles::sam::alignment::record::Flags;
        let mut r = RecordBuf::default();
        r.data_mut()
            .insert(crate::TI_TAG.into(), Value::UInt32(taxon_id));
        r.data_mut()
            .insert(MD_TAG.into(), Value::String(md.as_bytes().into()));
        *r.sequence_mut() = Sequence::from(seq.to_vec());
        *r.cigar_mut() = Cigar::from(ops);
        *r.flags_mut() = Flags::empty(); // clear UNMAPPED so within_edit_limits processes it
        r
    }

    fn no_tree_params(taxon_ids: &AHashSet<u32>) -> FilterParams<'_> {
        FilterParams {
            taxon_ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: None,
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        }
    }

    #[test]
    fn test_on_taxa_record() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let record = record_with_ti(9606);
        let params = no_tree_params(&ids);
        assert_eq!(
            classify_record(&record, read_ti_tag(&record), &params),
            RecordDecision::OnTaxa
        );
    }

    #[test]
    fn test_reject_off_taxa_no_rescue() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let record = record_with_ti(1234);
        let params = no_tree_params(&ids);
        assert_eq!(
            classify_record(&record, read_ti_tag(&record), &params),
            RecordDecision::Reject
        );
    }

    #[test]
    fn test_edit_rescue_within_distance() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        // 1 SNP in MD, 10M cigar, no Ns in read.
        let record =
            record_with_ti_and_md(1234, "5A4", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.total(), 1);

        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert_eq!(
            classify_record(&record, read_ti_tag(&record), &params),
            RecordDecision::EditRescue
        );
    }

    #[test]
    fn test_count_edits_snps_only() {
        let record =
            record_with_ti_and_md(0, "5A2G2", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 2);
        assert!(edits.indels.is_empty());
        assert_eq!(edits.total(), 2);
    }

    #[test]
    fn test_count_edits_n_not_counted() {
        // MD says position 4 is a substitution; read base at position 4 is N
        // (a free pass against any A/C/G/T reference) -> 0 substitutions.
        let record = record_with_ti_and_md(0, "4A5", vec![Op::new(Kind::Match, 10)], b"ACGTNACGTA");
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 0);
    }

    #[test]
    fn test_count_edits_insertion() {
        let record = record_with_ti_and_md(
            0,
            "8",
            vec![
                Op::new(Kind::Match, 5),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 3),
            ],
            b"ACGTACCGTA",
        );
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 0);
        assert_eq!(edits.indels, vec![2i32]);
        assert_eq!(edits.total(), 1);
    }

    #[test]
    fn test_count_edits_deletion() {
        let record = record_with_ti_and_md(
            0,
            "5^AC5",
            vec![
                Op::new(Kind::Match, 5),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 5),
            ],
            b"ACGTAACGTA",
        );
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 0);
        assert_eq!(edits.indels, vec![-2i32]);
        assert_eq!(edits.total(), 1);
    }

    #[test]
    fn test_count_edits_n_in_softclip_not_subtracted() {
        // 5S + 10M; 1 MD substitution; 5 Ns in soft-clip, 0 Ns in aligned region.
        // Soft-clip Ns must NOT reduce the substitution count.
        let record = record_with_ti_and_md(
            0,
            "5A4",
            vec![Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 10)],
            b"NNNNNACGTAC",
        );
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 1);
    }

    #[test]
    fn test_count_edits_n_in_insertion_not_subtracted() {
        // 5M + 2I (NN) + 3M; 1 MD substitution in the 5M stretch; insertion Ns ignored.
        let record = record_with_ti_and_md(
            0,
            "4A0",
            vec![
                Op::new(Kind::Match, 5),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 3),
            ],
            b"ACGTANNACG",
        );
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 1);
    }

    #[test]
    fn test_count_edits_mixed_aligned_and_softclip_ns() {
        // 3S + 7M; 2 MD substitutions; 3 Ns in soft-clip AND 2 Ns at aligned positions.
        // Only the 2 aligned Ns should be subtracted -> net 0.
        //
        // Sequence NNNANCGNAT: soft-clip NNN, then aligned A·N·C·G·N·A·T (7M).
        // MD "1A2T2": 1 match, ref=A sub (aligned pos 1, read=N), 2 matches,
        //             ref=T sub (aligned pos 4, read=N), 2 matches; totals 7 ✓
        let record = record_with_ti_and_md(
            0,
            "1A2T2",
            vec![Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 7)],
            b"NNNANCGNAT",
        );
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 0); // 2 MD subs - 2 aligned Ns = 0
    }

    #[test]
    fn test_n_adjustment_softclip_ns_do_not_reduce_limit() {
        // 3S + 7M; 3 Ns are all in soft-clip, 0 Ns at aligned positions; 2 real substitutions.
        // rescue_n_adjustment=1, rescue_max_edit_distance=1 -> adjusted limit = 1 - 0 = 1.
        // The 2 real subs exceed the limit, so the record should be REJECTED.
        //
        // Under the old (all-Ns) logic: subs = max(0, 2−3) = 0; max_edits = max(0,1−3) = 0;
        // 0 ≤ 0 -> wrongly rescued.  New logic correctly rejects.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: Some(1),
            keep_unannotated: false,
        };
        // Sequence NNNAAGAGAT: soft-clip NNN, then aligned A·A·G·A·G·A·T (7M).
        // MD "2A1T2": 2 matches, ref=A sub (read=G), 1 match, ref=T sub (read=G), 2 matches
        //            ; 2 real subs, no aligned Ns.
        let record = record_with_ti_and_md(
            1234,
            "2A1T2",
            vec![Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, 7)],
            b"NNNAAGAGAT",
        );
        assert!(!within_edit_limits(&record, &params));
    }

    #[test]
    fn test_count_edits_no_md_tag_returns_none() {
        let record = RecordBuf::default();
        assert!(count_edits(&record).is_none());
    }

    #[test]
    fn test_read_ti_tag_int32() {
        let mut r = RecordBuf::default();
        r.data_mut()
            .insert(crate::TI_TAG.into(), Value::Int32(9606));
        assert_eq!(read_ti_tag(&r), Some(9606));
    }

    #[test]
    fn test_read_ti_tag_missing() {
        let r = RecordBuf::default();
        assert_eq!(read_ti_tag(&r), None);
    }

    #[test]
    fn test_read_ti_tag_all_integer_widths() {
        // Real SAM/BAM emitters can store the ti tag at any integer width.
        // Ensure all eight Value::Int* / UInt* variants decode correctly so
        // krak filter can consume tags written by any conforming tool.
        for v in [
            Value::Int8(7),
            Value::UInt8(7),
            Value::Int16(7),
            Value::UInt16(7),
            Value::Int32(7),
            Value::UInt32(7),
        ] {
            let mut r = RecordBuf::default();
            r.data_mut().insert(crate::TI_TAG.into(), v);
            assert_eq!(read_ti_tag(&r), Some(7));
        }
    }

    #[test]
    fn test_read_ti_tag_negative_signed_returns_none() {
        // Negative values are not valid taxon IDs; treat them as absent
        // rather than wrap to a huge u32 silently.
        for v in [Value::Int8(-1), Value::Int16(-1), Value::Int32(-1)] {
            let mut r = RecordBuf::default();
            r.data_mut().insert(crate::TI_TAG.into(), v);
            assert_eq!(read_ti_tag(&r), None);
        }
    }

    #[test]
    fn test_read_ti_tag_string_value_returns_none() {
        // A `ti:Z:9606` tag (string-typed) is unsupported; Kraken always
        // emits integer taxon IDs. We must NOT silently parse the string;
        // return None so the record is treated as unannotated.
        let mut r = RecordBuf::default();
        r.data_mut()
            .insert(crate::TI_TAG.into(), Value::String(b"9606".to_vec().into()));
        assert_eq!(read_ti_tag(&r), None);
    }

    #[test]
    fn test_filter_by_record_counts() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        let records = vec![
            record_with_ti(9606),
            record_with_ti(1234),
            record_with_ti(9606),
        ];
        let mut kept = 0u32;
        let mut rejected = 0u32;
        let metric = filter_by_record(
            records.into_iter().map(Ok),
            &params,
            |_| {
                kept += 1;
                Ok(())
            },
            |_| {
                rejected += 1;
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(metric.on_taxa, 2);
        assert_eq!(metric.num_kept, 2);
        assert_eq!(metric.num_filtered, 1);
        assert_eq!(metric.total, 3);
        assert_eq!(kept, 2);
        assert_eq!(rejected, 1);
    }

    #[test]
    fn test_filter_by_record_missing_ti_counted() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        let records = vec![RecordBuf::default(), record_with_ti(9606)];
        let metric =
            filter_by_record(records.into_iter().map(Ok), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.missing_ti_tag, Some(1));
        assert_eq!(metric.total, 2);
    }

    #[test]
    fn test_filter_metrics_fraction() {
        let ids: AHashSet<u32> = [1u32].into_iter().collect();
        let params = no_tree_params(&ids);
        let records = vec![record_with_ti(2), record_with_ti(2)];
        let metric =
            filter_by_record(records.into_iter().map(Ok), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.frac_removed, 1.0);
        assert_eq!(metric.num_filtered, 2);
    }

    #[test]
    fn test_max_indel_count_limit() {
        // 2 indel events, but rescue_max_indels = 1 -> reject.
        let record = record_with_ti_and_md(
            1234,
            "3^A3^C3",
            vec![
                Op::new(Kind::Match, 3),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 3),
                Op::new(Kind::Deletion, 1),
                Op::new(Kind::Match, 3),
            ],
            b"ACGTACGTA",
        );
        // For this test we test count_edits directly (unmapped record skips within_edit_limits).
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.indel_count(), 2);
    }

    #[test]
    fn test_count_edits_aligned_ns_only_free_at_substitution_position() {
        // Under the joint MD/CIGAR walk, aligned Ns are "free" only at the
        // exact positions MD calls a substitution. Read has 4 leading Ns and
        // an MD substitution at position 5 (an 'A' read base) -> subs = 1.
        let record =
            record_with_ti_and_md(1234, "5A4", vec![Op::new(Kind::Match, 10)], b"NNNNACGTAC");
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.substitutions, 1);
    }

    #[test]
    fn test_filter_by_template_on_taxa() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        // Two reads in same template: one on-taxa, one off-taxa.
        let mut r1 = record_with_ti(9606);
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 = record_with_ti(1234);
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let records = vec![r1, r2];
        let mut kept = 0usize;
        let (rec_metric, tmpl_metric) = filter_by_template(
            records.into_iter().map(Ok),
            &params,
            |_| {
                kept += 1;
                Ok(())
            },
            |_| Ok(()),
        )
        .unwrap();

        assert_eq!(tmpl_metric.on_taxa, 1);
        assert_eq!(tmpl_metric.total, 1);
        assert_eq!(rec_metric.on_taxa, 2); // both reads in kept template counted as on-taxa
        assert_eq!(kept, 2);
    }

    #[test]
    fn test_classify_template_ignores_supplementary_taxon_for_on_taxa() {
        // B3: a template whose primary record is off-taxa but whose
        // supplementary alignment carries an on-taxa ti tag must NOT be
        // classified as OnTaxa; supplementary/secondary taxon IDs are
        // ignored for classification, consistent with the edit-rescue path.
        use noodles::sam::alignment::record::Flags;
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);

        let mut primary = record_with_ti(1234); // off-taxa
        *primary.name_mut() = Some(b"read1".as_ref().into());
        *primary.flags_mut() = Flags::empty();

        let mut supp = record_with_ti(9606); // on-taxa, but supplementary
        *supp.name_mut() = Some(b"read1".as_ref().into());
        *supp.flags_mut() = Flags::SUPPLEMENTARY;

        let taxon_ids: Vec<Option<u32>> = vec![Some(1234), Some(9606)];
        let template_refs: Vec<&RecordBuf> = vec![&primary, &supp];

        let decision = classify_template(&template_refs, &taxon_ids, &params);
        assert_eq!(
            decision,
            RecordDecision::Reject,
            "supplementary on-taxa taxon must not rescue an off-taxa primary"
        );
    }

    #[test]
    fn test_classify_template_ignores_secondary_taxon_for_on_taxa() {
        // Same as above but with SECONDARY flag instead of SUPPLEMENTARY.
        use noodles::sam::alignment::record::Flags;
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);

        let mut primary = record_with_ti(1234);
        *primary.name_mut() = Some(b"read1".as_ref().into());
        *primary.flags_mut() = Flags::empty();

        let mut secondary = record_with_ti(9606);
        *secondary.name_mut() = Some(b"read1".as_ref().into());
        *secondary.flags_mut() = Flags::SECONDARY;

        let taxon_ids: Vec<Option<u32>> = vec![Some(1234), Some(9606)];
        let template_refs: Vec<&RecordBuf> = vec![&primary, &secondary];

        let decision = classify_template(&template_refs, &taxon_ids, &params);
        assert_eq!(
            decision,
            RecordDecision::Reject,
            "secondary on-taxa taxon must not rescue an off-taxa primary"
        );
    }

    #[test]
    fn test_within_edit_limits_no_rescue_requested_returns_false_even_for_zero_edits() {
        // B4: with default --max-edit-distance 0 and no other rescue knobs set,
        // edit rescue must be disabled entirely. A record with 0 edits must
        // NOT be rescued just because it satisfies the trivially-zero default.
        let record =
            record_with_ti_and_md(1234, "10", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        let edits = count_edits(&record).unwrap();
        assert_eq!(edits.total(), 0, "fixture must have zero edits");

        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);

        assert!(
            !within_edit_limits(&record, &params),
            "0-edit record must not be rescued when no rescue knob is set"
        );
    }

    #[test]
    fn test_within_edit_limits_explicit_zero_distance_with_max_indels_allows_rescue() {
        // Opting in via --max-indels makes the rescue path active, and a
        // zero-edit, zero-indel record then passes.
        let record =
            record_with_ti_and_md(1234, "10", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(0),
            rescue_max_indels: Some(0),
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert!(
            within_edit_limits(&record, &params),
            "0-edit record must be rescued when user opts in via --rescue-max-indels"
        );
    }

    #[test]
    fn test_fastq_base_name_strips_slash1() {
        assert_eq!(fastq_base_name(b"read1/1").unwrap(), "read1");
    }

    #[test]
    fn test_fastq_base_name_strips_slash2() {
        assert_eq!(fastq_base_name(b"read1/2").unwrap(), "read1");
    }

    #[test]
    fn test_fastq_base_name_strips_description() {
        // Illumina-style "read 1:N:0:ATCG"; keep only the first token.
        assert_eq!(fastq_base_name(b"read1 1:N:0:ATCG").unwrap(), "read1");
        assert_eq!(fastq_base_name(b"read1/1 desc").unwrap(), "read1");
    }

    #[test]
    fn test_fastq_base_name_no_suffix() {
        assert_eq!(fastq_base_name(b"read1").unwrap(), "read1");
    }

    #[test]
    fn test_fastq_base_name_non_utf8_errors() {
        let err = fastq_base_name(&[0xFF, 0xFE]).unwrap_err();
        assert!(format!("{err:#}").contains("non-UTF-8"));
    }

    #[test]
    fn test_classify_by_taxon_ancestor_rescue() {
        // tree: 1 -> 9606. Target = 9606, allow_ancestors=true.
        // A record annotated with taxon 1 (root) is an ancestor of 9606.
        use crate::kraken_report::KrakenReportEntry;
        let entries = vec![
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "R".to_owned(),
                taxon_id: 1,
                name: "root".to_owned(),
                indent: 0,
            },
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "S".to_owned(),
                taxon_id: 9606,
                name: "Homo sapiens".to_owned(),
                indent: 2,
            },
        ];
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: Some(&tree),
            allow_ancestors: true,
            rescue_max_edit_distance: None,
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert_eq!(
            classify_by_taxon(1, &params),
            Some(RecordDecision::AncestorRescue)
        );
        // Without --allow-ancestors, no rescue.
        let no_anc = FilterParams {
            allow_ancestors: false,
            ..params
        };
        assert_eq!(classify_by_taxon(1, &no_anc), None);
    }

    #[test]
    fn test_within_edit_limits_unmapped_returns_false() {
        // Unmapped records have no MD tag and no alignment context;
        // within_edit_limits must reject them outright (line ~1712).
        use noodles::sam::alignment::record::Flags;
        let mut record =
            record_with_ti_and_md(0, "5A4", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        *record.flags_mut() = Flags::UNMAPPED;
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert!(!within_edit_limits(&record, &params));
    }

    #[test]
    fn test_within_edit_limits_no_md_tag_returns_false() {
        // A mapped record without an MD tag cannot be edit-rescued
        // (line ~1715).
        use noodles::sam::alignment::record::Flags;
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::empty();
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert!(!within_edit_limits(&r, &params));
    }

    #[test]
    fn test_within_edit_limits_max_indels_rejects() {
        // CIGAR carries 2 insertions; with --max-indels 1 the record must
        // be rejected from edit-rescue (line ~1732).
        let record = record_with_ti_and_md(
            0,
            "5",
            vec![
                Op::new(Kind::Match, 5),
                Op::new(Kind::Insertion, 1),
                Op::new(Kind::Match, 0),
                Op::new(Kind::Insertion, 1),
            ],
            b"ACGTAACGTA",
        );
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100),
            rescue_max_indels: Some(1),
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert!(!within_edit_limits(&record, &params));
    }

    #[test]
    fn test_within_edit_limits_max_indel_length_rejects() {
        // Single 5bp insertion exceeds --max-indel-length 2 (line ~1737).
        let record = record_with_ti_and_md(
            0,
            "5",
            vec![Op::new(Kind::Match, 5), Op::new(Kind::Insertion, 5)],
            b"ACGTACCCCC",
        );
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100),
            rescue_max_indels: None,
            rescue_max_indel_length: Some(2),
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        assert!(!within_edit_limits(&record, &params));
    }

    #[test]
    fn test_within_edit_limits_n_adjustment_lowers_threshold() {
        // 4 aligned Ns + rescue_n_adjustment=2 reduces threshold by 2.
        // 1 MD substitution remains a single "edit"; with original
        // threshold 1 -> still passes; with threshold 1 - 2 = saturating
        // 0 -> fails.
        let record = record_with_ti_and_md(0, "5A4", vec![Op::new(Kind::Match, 10)], b"NNNNAACGTA");
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params_pass = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: Some(8),
            keep_unannotated: false,
        };
        assert!(within_edit_limits(&record, &params_pass));

        let params_fail = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: Some(2),
            keep_unannotated: false,
        };
        assert!(!within_edit_limits(&record, &params_fail));
    }

    #[test]
    fn test_count_aligned_ns_skips_softclip_and_insertion() {
        // 5S + 5M with 5 N at the start (soft-clipped) and 1 N inside the
        // aligned region. Only the aligned N counts.
        let record = record_with_ti_and_md(
            0,
            "5",
            vec![Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, 5)],
            b"NNNNNAANCC",
        );
        assert_eq!(count_aligned_ns(&record), 1);
    }

    #[test]
    fn test_filter_by_record_ancestor_rescue() {
        // Per-record: an off-taxa annotation that is an ancestor of the
        // target tree gets routed to AncestorRescue and counted in
        // rescued_ancestors.
        use crate::kraken_report::KrakenReportEntry;
        let entries = vec![
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "R".to_owned(),
                taxon_id: 1,
                name: "root".to_owned(),
                indent: 0,
            },
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "S".to_owned(),
                taxon_id: 9606,
                name: "Homo sapiens".to_owned(),
                indent: 2,
            },
        ];
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: Some(&tree),
            allow_ancestors: true,
            rescue_max_edit_distance: None,
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        let r = record_with_ti(1); // root, ancestor of 9606
        let metric =
            filter_by_record(vec![Ok(r)].into_iter(), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.rescued_ancestors, 1);
        assert_eq!(metric.on_taxa, 0);
    }

    #[test]
    fn test_filter_by_record_edit_rescue_metric() {
        // Per-record: off-taxa, but within edit-distance threshold → counted
        // in rescued_variants.
        let r = record_with_ti_and_md(1234, "5A4", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };
        let metric =
            filter_by_record(vec![Ok(r)].into_iter(), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.rescued_variants, Some(1));
    }

    #[test]
    fn test_filter_by_record_unannotated_keep_metric() {
        // Per-record: no `ti` tag at all + --keep-unannotated → kept,
        // counted in unannotated_but_kept.
        let r = RecordBuf::default();
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let mut params = no_tree_params(&ids);
        params.keep_unannotated = true;
        let metric =
            filter_by_record(vec![Ok(r)].into_iter(), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.unannotated_but_kept, Some(1));
        assert_eq!(metric.missing_ti_tag, Some(1));
    }

    #[test]
    fn test_filter_by_record_reject_routes_to_reject_sink() {
        // Per-record: off-taxa, no rescue → reject sink.
        let r = record_with_ti(1234);
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = no_tree_params(&ids);
        let mut kept = 0usize;
        let mut rejected = 0usize;
        let metric = filter_by_record(
            vec![Ok(r)].into_iter(),
            &params,
            |_| {
                kept += 1;
                Ok(())
            },
            |_| {
                rejected += 1;
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(kept, 0);
        assert_eq!(rejected, 1);
        assert_eq!(metric.on_taxa, 0);
    }

    #[test]
    fn test_filter_by_template_ancestor_rescue() {
        // Template-level ancestor rescue + per-record metric counting.
        use crate::kraken_report::KrakenReportEntry;
        let entries = vec![
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "R".to_owned(),
                taxon_id: 1,
                name: "root".to_owned(),
                indent: 0,
            },
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "S".to_owned(),
                taxon_id: 9606,
                name: "Homo sapiens".to_owned(),
                indent: 2,
            },
        ];
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: Some(&tree),
            allow_ancestors: true,
            rescue_max_edit_distance: None,
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };

        let mut r1 = record_with_ti(1); // ancestor
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 = record_with_ti(1);
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let (rec, tmpl) = filter_by_template(
            vec![r1, r2].into_iter().map(Ok),
            &params,
            |_| Ok(()),
            |_| Ok(()),
        )
        .unwrap();
        assert_eq!(tmpl.rescued_ancestors, 1);
        assert_eq!(rec.rescued_ancestors, 2);
    }

    #[test]
    fn test_filter_by_template_edit_rescue() {
        // Template-level edit rescue: both primary records pass within-edit-limits.
        let mut r1 =
            record_with_ti_and_md(1234, "5A4", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 =
            record_with_ti_and_md(1234, "5A4", vec![Op::new(Kind::Match, 10)], b"ACGTACGTAC");
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let ids: AHashSet<u32> = [9606].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(1),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };

        let (rec, tmpl) = filter_by_template(
            vec![r1, r2].into_iter().map(Ok),
            &params,
            |_| Ok(()),
            |_| Ok(()),
        )
        .unwrap();
        assert_eq!(tmpl.rescued_variants, Some(1));
        assert_eq!(rec.rescued_variants, Some(2));
    }

    #[test]
    fn test_filter_by_template_writes_rejects_to_separate_sink() {
        // Off-taxa template with no rescue knobs → Reject; reject closure
        // must receive the records, the keep closure must NOT.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);

        let mut r1 = record_with_ti(1234);
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 = record_with_ti(5678);
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let records = vec![r1, r2];
        let mut kept = 0usize;
        let mut rejected = 0usize;
        let (rec_metric, tmpl_metric) = filter_by_template(
            records.into_iter().map(Ok),
            &params,
            |_| {
                kept += 1;
                Ok(())
            },
            |_| {
                rejected += 1;
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(kept, 0, "no records may be kept");
        assert_eq!(rejected, 2, "both records must hit the reject sink");
        assert_eq!(tmpl_metric.total, 1);
        assert_eq!(rec_metric.total, 2);
        assert_eq!(tmpl_metric.on_taxa, 0);
    }

    #[test]
    fn test_filter_by_template_unannotated_keep() {
        // A template whose primaries lack `ti` tags entirely is kept when
        // --keep-unannotated is set; counts must land in unannotated_but_kept.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let mut params = no_tree_params(&ids);
        params.keep_unannotated = true;

        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let records = vec![r1, r2];
        let mut kept = 0usize;
        let (rec_metric, tmpl_metric) = filter_by_template(
            records.into_iter().map(Ok),
            &params,
            |_| {
                kept += 1;
                Ok(())
            },
            |_| Ok(()),
        )
        .unwrap();

        assert_eq!(
            kept, 2,
            "both records must be kept under --keep-unannotated"
        );
        assert_eq!(tmpl_metric.unannotated_but_kept, Some(1));
        assert_eq!(rec_metric.unannotated_but_kept, Some(2));
        assert_eq!(tmpl_metric.missing_ti_tag, Some(1));
    }

    #[test]
    fn test_classify_template_fully_unannotated_does_not_edit_rescue() {
        // A template where NO record has a ti tag must not be classified as
        // EditRescue even if all primary reads pass edit limits. This must be
        // consistent with classify_record, which never edit-rescues a record
        // that has no taxon annotation.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100), // extremely lenient; any record would pass
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: true,
        };

        // Build two records with no ti tag at all.
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some(b"read1".as_ref().into());
        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some(b"read1".as_ref().into());

        let taxon_ids: Vec<Option<u32>> = vec![None, None];
        let template_refs: Vec<&RecordBuf> = vec![&r1, &r2];

        let decision = classify_template(&template_refs, &taxon_ids, &params);
        assert_eq!(
            decision,
            RecordDecision::UnannotatedKeep,
            "fully unannotated template must be UnannotatedKeep, not EditRescue"
        );
    }

    #[test]
    fn test_classify_template_fully_unannotated_rejects_when_keep_unannotated_false() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: Some(100),
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: false,
        };

        let r1 = RecordBuf::default();
        let r2 = RecordBuf::default();
        let taxon_ids: Vec<Option<u32>> = vec![None, None];
        let template_refs: Vec<&RecordBuf> = vec![&r1, &r2];

        let decision = classify_template(&template_refs, &taxon_ids, &params);
        assert_eq!(
            decision,
            RecordDecision::Reject,
            "fully unannotated template with keep_unannotated=false must be Reject"
        );
    }

    #[test]
    fn test_classify_template_partially_annotated_off_taxa_is_not_unannotated_keep() {
        // R1 has a ti tag pointing to an off-taxa taxon; R2 has no ti tag.
        // The template is NOT fully unannotated, so UnannotatedKeep must NOT
        // fire even when keep_unannotated=true. The correct outcome is Reject.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = FilterParams {
            taxon_ids: &ids,
            tree: None,
            allow_ancestors: false,
            rescue_max_edit_distance: None, // no edit rescue
            rescue_max_indels: None,
            rescue_max_indel_length: None,
            rescue_n_adjustment: None,
            keep_unannotated: true,
        };

        let r1 = RecordBuf::default();
        let r2 = RecordBuf::default();
        // R1 is annotated with an off-taxa taxon; R2 has no annotation.
        let taxon_ids: Vec<Option<u32>> = vec![Some(1234), None];
        let template_refs: Vec<&RecordBuf> = vec![&r1, &r2];

        let decision = classify_template(&template_refs, &taxon_ids, &params);
        assert_eq!(
            decision,
            RecordDecision::Reject,
            "partially-annotated template (some ti tags present) must be Reject, not UnannotatedKeep"
        );
    }

    #[test]
    fn test_metrics_fastx_optional_fields_are_none() {
        // A FASTX-mode metric initializes with None for SAM-specific fields.
        let metric = TaxaFilterMetric::default();
        assert!(metric.rescued_variants.is_none());
        assert!(metric.missing_ti_tag.is_none());
        assert!(metric.missing_classifications.is_none());
    }

    #[test]
    fn test_metrics_sam_mode_optional_fields_are_some() {
        // SAM filter_by_record initializes with Some(0).
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        let metric = filter_by_record(std::iter::empty(), &params, |_| Ok(()), |_| Ok(())).unwrap();
        assert_eq!(metric.rescued_variants, Some(0));
        assert_eq!(metric.missing_ti_tag, Some(0));
        assert!(metric.missing_classifications.is_none());
    }

    #[test]
    fn test_build_taxon_set_include_descendants() {
        use crate::kraken_report::{KrakenReportEntry, KrakenTaxonomyTree};
        let entries = vec![
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "R".to_owned(),
                taxon_id: 1,
                name: "root".to_owned(),
                indent: 0,
            },
            KrakenReportEntry {
                pct_fragments: 0.0,
                num_fragments_clade: 0,
                num_fragments_direct: 0,
                rank_code: "S".to_owned(),
                taxon_id: 9606,
                name: "human".to_owned(),
                indent: 2,
            },
        ];
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        let initial: AHashSet<u32> = [1u32].into_iter().collect();
        let result = build_taxon_set(&initial, Some(&tree), true, false).unwrap();
        assert!(result.contains(&1), "original target kept");
        assert!(result.contains(&9606), "descendant added");
        assert!(!result.contains(&0), "unclassified not inserted");
    }

    #[test]
    fn test_build_taxon_set_include_unclassified() {
        let initial: AHashSet<u32> = [9606u32].into_iter().collect();
        let result = build_taxon_set(&initial, None, false, true).unwrap();
        assert!(result.contains(&0), "taxon 0 inserted");
        assert!(result.contains(&9606), "original target kept");
    }

    #[test]
    fn test_build_taxon_set_no_flags() {
        let initial: AHashSet<u32> = [9606u32].into_iter().collect();
        let result = build_taxon_set(&initial, None, false, false).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result.contains(&9606));
    }

    #[test]
    fn test_build_taxon_set_unknown_taxon_errors() {
        use crate::kraken_report::{KrakenReportEntry, KrakenTaxonomyTree};
        let entries = vec![KrakenReportEntry {
            pct_fragments: 0.0,
            num_fragments_clade: 0,
            num_fragments_direct: 0,
            rank_code: "R".to_owned(),
            taxon_id: 1,
            name: "root".to_owned(),
            indent: 0,
        }];
        let tree = KrakenTaxonomyTree::from_entries(&entries).unwrap();
        let initial: AHashSet<u32> = [9999u32].into_iter().collect();
        let err = build_taxon_set(&initial, Some(&tree), false, false).unwrap_err();
        assert!(
            err.to_string().contains("9999"),
            "error mentions missing id: {err}"
        );
    }

    #[test]
    fn test_classify_fastx_on_taxa() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        assert_eq!(classify_fastx(9606, &params), RecordDecision::OnTaxa);
    }

    #[test]
    fn test_classify_fastx_reject() {
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        assert_eq!(classify_fastx(1234, &params), RecordDecision::Reject);
    }

    #[test]
    fn test_classify_fastx_unclassified_in_set() {
        // After build_taxon_set inserts 0, taxon_id 0 -> OnTaxa.
        let ids: AHashSet<u32> = [0u32].into_iter().collect();
        let params = no_tree_params(&ids);
        assert_eq!(classify_fastx(0, &params), RecordDecision::OnTaxa);
    }

    #[test]
    fn test_classify_fastx_missing_from_map_becomes_taxon_0_and_rejects() {
        // A read absent from kraken_map gets taxon_id 0; if 0 not in set -> Reject.
        let ids: AHashSet<u32> = [9606u32].into_iter().collect();
        let params = no_tree_params(&ids);
        // Simulate lookup miss by using taxon_id 0 directly.
        assert_eq!(classify_fastx(0, &params), RecordDecision::Reject);
    }

    fn make_ti_record(name: &str, taxon_id: u32) -> noodles::sam::alignment::record_buf::RecordBuf {
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };
        let mut r = RecordBuf::default();
        *r.name_mut() = Some(name.as_bytes().into());
        *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
        *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
        r.data_mut()
            .insert(crate::TI_TAG.into(), Value::UInt32(taxon_id));
        r
    }

    /// Regression: a mapped CRAM (with `@SQ` in its header) without
    /// `--cram-reference` used to panic deep inside noodles' decoder
    /// ("invalid slice reference sequence name") when filter tried to read
    /// records. The fix bails with a clear error after reading the header.
    #[test]
    fn test_run_filter_mapped_cram_without_reference_errors_cleanly() {
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

        let out_cram = dir.path().join("out.cram");
        let err = run_filter(FilterArgs {
            input: in_cram,
            output: out_cram,
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            ..FilterArgs::default_for_test()
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--cram-reference"), "got: {msg}");
        assert!(msg.contains("reference sequences"), "got: {msg}");
    }

    #[test]
    fn test_filter_without_metrics_file_succeeds_and_writes_no_file() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::RecordBuf;

        let dir = tempfile::TempDir::new().unwrap();

        let in_cram = dir.path().join("input.cram");
        let header = sam::Header::default();
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &make_ti_record("r1", 9606))
                .unwrap();
            w.write_alignment_record(&header, &make_ti_record("r2", 1234))
                .unwrap();
            w.try_finish(&header).unwrap();
        }

        let out_cram = dir.path().join("output.cram");

        run_filter(FilterArgs {
            input: in_cram,
            output: out_cram.clone(),
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            ..FilterArgs::default_for_test()
        })
        .unwrap();

        let mut reader = crate::open_cram_reader(&out_cram, None).unwrap();
        let out_header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .records(&out_header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name().unwrap().as_ref() as &[u8], b"r1");
    }

    #[test]
    fn test_filter_cram_basic() {
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::RecordBuf;

        let dir = tempfile::TempDir::new().unwrap();

        // Write input CRAM: r1->9606, r2->1234
        let in_cram = dir.path().join("input.cram");
        let header = sam::Header::default();
        {
            let mut w = crate::open_cram_writer(&in_cram, None).unwrap();
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &make_ti_record("r1", 9606))
                .unwrap();
            w.write_alignment_record(&header, &make_ti_record("r2", 1234))
                .unwrap();
            w.try_finish(&header).unwrap();
        }

        let out_cram = dir.path().join("output.cram");
        let metrics_path = dir.path().join("metrics.json");

        run_filter(FilterArgs {
            input: in_cram,
            output: out_cram.clone(),
            metrics: Some(metrics_path),
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            ..FilterArgs::default_for_test()
        })
        .unwrap();

        // Only r1 (taxon 9606) should appear in output
        let mut reader = crate::open_cram_reader(&out_cram, None).unwrap();
        let out_header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .records(&out_header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        let name = records[0].name().unwrap();
        assert_eq!(name.as_ref() as &[u8], b"r1");
    }

    #[test]
    fn test_run_filter_extensionless_path_with_gzipped_fastq() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::{Read as _, Write as _};

        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("piped"); // no extension
        {
            let f = std::fs::File::create(&in_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
                .unwrap();
            enc.finish().unwrap();
        }
        let kraken = dir.path().join("kraken.tsv");
        std::fs::write(&kraken, "C\tr1\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\n").unwrap();

        let out = dir.path().join("out.fq");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let mut got = String::new();
        std::fs::File::open(&out)
            .unwrap()
            .read_to_string(&mut got)
            .unwrap();
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_run_filter_extensionless_path_with_bam() {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };

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
            r.data_mut()
                .insert(crate::TI_TAG.into(), Value::UInt32(9606));
            w.write_alignment_record(&header, &r).unwrap();
        }
        let in_path = dir.path().join("piped"); // no extension
        std::fs::rename(&bam_path, &in_path).unwrap();

        let out = dir.path().join("out.bam");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out.clone(),
            taxon_ids: taxa,
            per_record: true,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();
        // Verify the record was kept.
        let mut reader = bam::io::reader::Builder.build_from_path(&out).unwrap();
        let hdr = reader.read_header().unwrap();
        let recs: Vec<RecordBuf> = reader
            .record_bufs(&hdr)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name().unwrap().as_ref() as &[u8], b"r1");
    }

    #[test]
    fn test_run_filter_bam_extension_unambiguous_path() {
        // Exercises the .bam-extension fast path in run_filter_dispatch
        // (skips the sniffer entirely, opens via open_bam_reader).
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };

        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("input.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_path).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("kept".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            r.data_mut()
                .insert(crate::TI_TAG.into(), Value::Int32(9606));
            w.write_alignment_record(&header, &r).unwrap();
            let mut r2 = RecordBuf::default();
            *r2.name_mut() = Some("rejected".as_bytes().into());
            *r2.sequence_mut() = Sequence::from(b"TTTT".to_vec());
            *r2.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            r2.data_mut()
                .insert(crate::TI_TAG.into(), Value::Int32(1234));
            w.write_alignment_record(&header, &r2).unwrap();
        }

        let out = dir.path().join("out.bam");
        let rejects = dir.path().join("rejects.bam");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out.clone(),
            taxon_ids: taxa,
            rejects: Some(rejects.clone()),
            per_record: true,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let mut keep = bam::io::reader::Builder.build_from_path(&out).unwrap();
        let h = keep.read_header().unwrap();
        let recs: Vec<RecordBuf> = keep
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name().unwrap().as_ref() as &[u8], b"kept");

        let mut rej = bam::io::reader::Builder.build_from_path(&rejects).unwrap();
        let h = rej.read_header().unwrap();
        let recs: Vec<RecordBuf> = rej
            .record_bufs(&h)
            .collect::<std::io::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name().unwrap().as_ref() as &[u8], b"rejected");
    }

    #[test]
    fn test_run_filter_sam_extension_unambiguous_path() {
        // .sam-extension fast-path through run_filter -> filter_sam.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("input.sam");
        {
            let mut f = std::fs::File::create(&in_path).unwrap();
            writeln!(f, "@HD\tVN:1.6").unwrap();
            // ti:i:9606 → on-taxa for 9606
            writeln!(f, "kept\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tti:i:9606").unwrap();
            writeln!(f, "rejected\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tti:i:1234").unwrap();
        }

        let out = dir.path().join("out.sam");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out.clone(),
            taxon_ids: taxa,
            per_record: true,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let body = std::fs::read_to_string(&out).unwrap();
        assert!(
            body.contains("kept"),
            "expected kept record in output: {body}"
        );
        assert!(
            !body.contains("rejected"),
            "rejected must not be in output: {body}"
        );
    }

    #[test]
    fn test_run_filter_classifications_with_alignment_extension_errors() {
        // --classifications is only valid for FASTX input. .bam input + -k must
        // bail with a clear error.
        use noodles::bam;
        use noodles::sam;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("input.bam");
        {
            let mut w = bam::io::writer::Builder.build_from_path(&in_path).unwrap();
            let header = sam::Header::default();
            w.write_header(&header).unwrap();
        }
        let kraken_path = dir.path().join("kraken.tsv");
        std::fs::write(&kraken_path, b"").unwrap();

        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out.bam"),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("--classifications"));
    }

    #[test]
    fn test_run_filter_fasta_without_classifications_errors() {
        // FASTA input requires --classifications (-c). Omitting it must error.
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("input.fa");
        std::fs::write(&in_path, b">r1\nACGT\n").unwrap();

        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out.fa"),
            taxon_ids: taxa,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("--classifications"));
    }

    #[test]
    fn test_run_filter_unknown_sniff_errors() {
        // A file whose first bytes match no known signature must produce a
        // clear error (not a panic) when filter is asked to consume it.
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("garbage");
        std::fs::write(&in_path, b"\x00\x01\x02not-a-known-format\n").unwrap();

        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out"),
            taxon_ids: taxa,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("could not infer format") || msg.contains("Unknown"),
            "expected unknown-format error, got: {msg}"
        );
    }

    #[test]
    fn test_filter_fastq_gz_roundtrip() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::{Read as _, Write as _};

        let dir = tempfile::TempDir::new().unwrap();

        let in_path = dir.path().join("in.fq.gz");
        {
            let f = std::fs::File::create(&in_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
                .unwrap();
            enc.finish().unwrap();
        }

        let kraken_path = dir.path().join("kraken.tsv");
        std::fs::write(
            &kraken_path,
            "C\tr1\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\n",
        )
        .unwrap();

        let out_path = dir.path().join("out.fq.gz");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let f = std::fs::File::open(&out_path).unwrap();
        let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
        let mut got = String::new();
        dec.read_to_string(&mut got).unwrap();
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n");
    }

    /// Write a `.fq.gz` input, kraken assignments file, and run `run_filter`
    /// keeping taxon 9606 with the given `threads` + `compression_level`.
    /// Returns the (compressed file size, decompressed content) for the output.
    fn run_filter_fastq_gz(
        tmpdir: &std::path::Path,
        out_name: &str,
        threads: usize,
        compression_level: u32,
        n_records: usize,
    ) -> (u64, String) {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::{Read as _, Write as _};

        // Synthesize n_records FASTQ records with deterministic but varied
        // payloads so compression has something to work on. Half assigned to
        // 9606 (kept), half to 1234 (dropped) — exercises both branches.
        let mut fastq = String::new();
        let mut kraken = String::new();
        for i in 0..n_records {
            let name = format!("r{i}");
            let bases: String = (0..96)
                .map(|j| match (i + j) % 4 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                })
                .collect();
            let quals: String = "I".repeat(bases.len());
            fastq.push_str(&format!("@{name}\n{bases}\n+\n{quals}\n"));
            let taxon = if i % 2 == 0 { 9606 } else { 1234 };
            kraken.push_str(&format!("C\t{name}\t{taxon}\t{}\t{taxon}:1\n", bases.len()));
        }

        let in_path = tmpdir.join("in.fq.gz");
        {
            let f = std::fs::File::create(&in_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(fastq.as_bytes()).unwrap();
            enc.finish().unwrap();
        }
        let kraken_path = tmpdir.join("kraken.tsv");
        std::fs::write(&kraken_path, &kraken).unwrap();

        let out_path = tmpdir.join(out_name);
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            threads,
            compression_level,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let size = std::fs::metadata(&out_path).unwrap().len();
        let f = std::fs::File::open(&out_path).unwrap();
        let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
        let mut got = String::new();
        dec.read_to_string(&mut got).unwrap();
        (size, got)
    }

    #[test]
    fn test_filter_fastq_gz_higher_compression_level_yields_smaller_file() {
        // At a fixed thread count, raising --compression-level must shrink
        // the file while leaving the decoded payload byte-for-byte identical.
        let dir = tempfile::TempDir::new().unwrap();
        let (size_low, content_low) = run_filter_fastq_gz(dir.path(), "low.fq.gz", 1, 1, 200);
        let (size_high, content_high) = run_filter_fastq_gz(dir.path(), "high.fq.gz", 1, 9, 200);

        assert_eq!(
            content_low, content_high,
            "decoded output must match across compression levels"
        );
        assert!(
            size_high < size_low,
            "expected level 9 ({size_high} bytes) < level 1 ({size_low} bytes)"
        );
    }

    #[test]
    fn test_filter_fastq_gz_threads_one_and_many_produce_identical_output() {
        // At a fixed compression level, --threads must not change the decoded
        // payload. (File sizes may differ slightly because the worker pool
        // schedules block boundaries differently — that's expected and fine.)
        let dir = tempfile::TempDir::new().unwrap();
        let (_, content_serial) = run_filter_fastq_gz(dir.path(), "t1.fq.gz", 1, 5, 200);
        let (_, content_parallel) = run_filter_fastq_gz(dir.path(), "t4.fq.gz", 4, 5, 200);

        assert_eq!(
            content_serial, content_parallel,
            "decoded output must match across thread counts"
        );
    }

    #[test]
    fn test_filter_fasta_gz_input_plain_output() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();

        let in_path = dir.path().join("in.fa.gz");
        {
            let f = std::fs::File::create(&in_path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b">s1\nACGT\n>s2\nTTTT\n").unwrap();
            enc.finish().unwrap();
        }

        let kraken_path = dir.path().join("kraken.tsv");
        std::fs::write(
            &kraken_path,
            "C\ts1\t9606\t4\t9606:1\nC\ts2\t1234\t4\t1234:1\n",
        )
        .unwrap();

        let out_path = dir.path().join("out.fa");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let got = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(got, ">s1\nACGT\n");
    }

    /// Helper: write a plain FASTQ and a Kraken assignments file, then run
    /// `run_filter` with the given `unordered` flag. Returns the output path
    /// contents.
    fn run_fastq_filter(unordered: bool, fastq: &str, kraken: &str) -> String {
        use std::io::Write as _;

        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, fastq).unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        f.write_all(kraken.as_bytes()).unwrap();

        let out_path = dir.path().join("out.fq");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            unordered,
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();
        std::fs::read_to_string(&out_path).unwrap()
    }

    #[test]
    fn test_filter_fastq_first_record_slash2_errors() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, "@r1/2\nACGT\n+\nIIII\n").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        std::fs::write(&kraken_path, "C\tr1\t9606\t4\t9606:1\n").unwrap();
        let out_path = dir.path().join("out.fq");
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path,
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: false,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("ends with /2"));
    }

    #[test]
    fn test_filter_fastq_missing_read_counts_against_metric() {
        // r2 is absent from the Kraken file → missing_classifications bumps and
        // the read is treated as taxon 0 (rejected unless a target matches 0).
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(
            &in_path,
            "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n",
        )
        .unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        // r2 absent; r1 and r3 are 9606 so they survive.
        f.write_all(b"C\tr1\t9606\t4\t9606:1\nC\tr3\t9606\t4\t9606:1\n")
            .unwrap();
        let out_path = dir.path().join("out.fq");
        let metrics_path = dir.path().join("metrics.tsv");

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            metrics: Some(metrics_path.clone()),
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();
        let got = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n@r3\nGGGG\n+\nKKKK\n");
        let metrics = std::fs::read_to_string(&metrics_path).unwrap();
        // missing_classifications column should be 1.
        assert!(
            metrics.lines().any(|l| l.split('\t').any(|c| c == "1")),
            "metrics:\n{metrics}"
        );
    }

    #[test]
    fn test_filter_fastq_streaming_disagreeing_buffered_taxa_errors() {
        // Out-of-order Kraken file: lookups for r1 stream past r2 and r2 (with
        // disagreeing taxa for the same name); the streaming lookup buffers
        // the first r2 entry, then sees a conflicting second r2 entry and
        // bails.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, "@r1\nACGT\n+\nIIII\n").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        // Streaming lookup must consume r2:9606 (buffered as out-of-order),
        // then r2:1234 with a different taxon → conflict.
        f.write_all(b"C\tr2\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\nC\tr1\t9606\t4\t9606:1\n")
            .unwrap();
        let out_path = dir.path().join("out.fq");

        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path,
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("disagreeing taxon"));
    }

    #[test]
    fn test_filter_fastq_empty_input_returns_empty_metric() {
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, b"").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        std::fs::write(&kraken_path, b"").unwrap();
        let out_path = dir.path().join("out.fq");
        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: [9606u32].into_iter().collect(),
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();
        assert_eq!(std::fs::read_to_string(&out_path).unwrap(), "");
    }

    #[test]
    fn test_filter_fastq_streaming_in_order() {
        // Default `unordered=false`: streaming lookup with assignments in the
        // same order as the FASTQ input.
        let fastq = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n";
        let kraken = "C\tr1\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\nC\tr3\t9606\t4\t9606:1\n";
        let got = run_fastq_filter(false, fastq, kraken);
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n@r3\nGGGG\n+\nKKKK\n");
    }

    #[test]
    fn test_filter_fastq_streaming_handles_modest_disorder() {
        // Default `unordered=false`: streaming lookup absorbs out-of-order
        // assignments via the lookahead buffer (here r3 is emitted before r2
        // in the kraken file, but both reads are still resolved correctly).
        let fastq = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n";
        let kraken = "C\tr1\t9606\t4\t9606:1\nC\tr3\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\n";
        let got = run_fastq_filter(false, fastq, kraken);
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n@r3\nGGGG\n+\nKKKK\n");
    }

    #[test]
    fn test_filter_fastq_unordered_matches_streaming() {
        // `unordered=true` loads the full map upfront. With completely
        // reversed assignments order, the unordered path still yields the
        // same selection as the streaming path.
        let fastq = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n";
        let kraken = "C\tr3\t9606\t4\t9606:1\nC\tr2\t1234\t4\t1234:1\nC\tr1\t9606\t4\t9606:1\n";
        let got = run_fastq_filter(true, fastq, kraken);
        assert_eq!(got, "@r1\nACGT\n+\nIIII\n@r3\nGGGG\n+\nKKKK\n");
    }

    #[test]
    fn test_run_filter_n_adjustment_zero_errors() {
        // run_filter must reject --n-adjustment 0 outright (line ~220).
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fa");
        std::fs::write(&in_path, b">r1\nACGT\n").unwrap();
        let kraken_path = dir.path().join("k.tsv");
        std::fs::write(&kraken_path, b"C\tr1\t9606\t4\t9606:1\n").unwrap();
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out.fa"),
            taxon_ids: taxa,
            rescue_n_adjustment: Some(0),
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("--rescue-n-adjustment must be >= 1"));
    }

    #[test]
    fn test_validate_tree_requirements_descendants_without_tree_errors() {
        let mut args = FilterArgs::default_for_test();
        args.include_descendants = true;
        let tree: Option<KrakenTaxonomyTree> = None;
        let err = validate_tree_requirements(&tree, &args).unwrap_err();
        assert!(format!("{err:#}").contains("--include-descendants"));
    }

    #[test]
    fn test_validate_tree_requirements_ancestors_without_tree_errors() {
        let mut args = FilterArgs::default_for_test();
        args.allow_ancestors = true;
        let tree: Option<KrakenTaxonomyTree> = None;
        let err = validate_tree_requirements(&tree, &args).unwrap_err();
        assert!(format!("{err:#}").contains("--allow-ancestors"));
    }

    #[test]
    fn test_lookup_state_stream_missing_file_errors() {
        // LookupState::from_source(FastxSource::Stream(...)) must surface a
        // clear error when the assignments file is missing.
        let result = LookupState::from_source(FastxSource::Stream(std::path::Path::new(
            "/nonexistent.tsv",
        )));
        assert!(result.is_err());
        let err = match result {
            Ok(_) => unreachable!(),
            Err(e) => e,
        };
        assert!(format!("{err:#}").contains("failed to open"));
    }

    impl FilterArgs {
        fn default_for_test() -> Self {
            FilterArgs {
                input: std::path::PathBuf::new(),
                output: std::path::PathBuf::new(),
                kraken_report: None,
                metrics: None,
                taxon_ids: AHashSet::new(),
                rejects: None,
                allow_ancestors: false,
                rescue_max_edit_distance: None,
                rescue_max_indels: None,
                rescue_max_indel_length: None,
                rescue_n_adjustment: None,
                per_record: false,
                classifications: None,
                include_descendants: false,
                include_unclassified: false,
                cram_reference: None,
                keep_unannotated: false,
                unordered: false,
                threads: 1,
                compression_level: 5,
            }
        }
    }

    #[test]
    fn test_run_filter_fastq_interleaved_with_rejects() {
        // Interleaved FASTQ filter with a rejects path: r1 is on-taxa,
        // r2 is off-taxa. Both must be processed as a single template.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(
            &in_path,
            "@kept/1\nACGT\n+\nIIII\n@kept/2\nTGCA\n+\nIIII\n@drop/1\nTTTT\n+\nJJJJ\n@drop/2\nAAAA\n+\nJJJJ\n",
        )
        .unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        writeln!(f, "C\tkept\t9606\t150\t9606:1").unwrap();
        writeln!(f, "C\tdrop\t1234\t150\t1234:1").unwrap();
        drop(f);

        let out = dir.path().join("out.fq");
        let rejects = dir.path().join("rejects.fq");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out.clone(),
            taxon_ids: taxa,
            rejects: Some(rejects.clone()),
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let kept = std::fs::read_to_string(&out).unwrap();
        let rej = std::fs::read_to_string(&rejects).unwrap();
        assert!(
            kept.contains("@kept/1"),
            "kept output should include kept/1: {kept}"
        );
        assert!(
            kept.contains("@kept/2"),
            "kept output should include kept/2: {kept}"
        );
        assert!(
            !kept.contains("@drop"),
            "kept output must not include drop: {kept}"
        );
        assert!(
            rej.contains("@drop/1"),
            "rejects should include drop/1: {rej}"
        );
        assert!(
            rej.contains("@drop/2"),
            "rejects should include drop/2: {rej}"
        );
    }

    #[test]
    fn test_run_filter_fastq_interleaved_name_mismatch_errors() {
        // Detected interleaved (first record has /1) but R2 is from a
        // different template; error path at line ~889.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, "@a/1\nACGT\n+\nIIII\n@b/2\nTGCA\n+\nIIII\n").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        writeln!(f, "C\ta\t9606\t150\t9606:1").unwrap();
        drop(f);

        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out.fq"),
            taxon_ids: taxa,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("not immediately followed by") || msg.contains("interleaved FASTQ"),
            "expected interleaved name-mismatch error, got: {msg}"
        );
    }

    #[test]
    fn test_run_filter_fastq_interleaved_truncated_after_r1_errors() {
        // Detected interleaved (first record has /1) but the file ends after
        // r1; error path at line ~881.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fq");
        std::fs::write(&in_path, "@a/1\nACGT\n+\nIIII\n").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        writeln!(f, "C\ta\t9606\t150\t9606:1").unwrap();
        drop(f);

        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);
        let err = super::run_filter(super::FilterArgs {
            input: in_path,
            output: dir.path().join("out.fq"),
            taxon_ids: taxa,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("has no following /2") || msg.contains("truncated"),
            "expected truncated-interleaved error, got: {msg}"
        );
    }

    #[test]
    fn test_fastx_sink_plain_flush_passes_through() {
        use std::io::Write;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("plain.fq");
        let mut sink = FastxSink::create(&p, 1, 5).unwrap();
        sink.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        sink.flush().unwrap();
        sink.finalize().unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_fastx_sink_gz_finalize_emits_footer() {
        // Without `finalize` the gzip stream might be missing its footer; this
        // test ensures finalize() emits a complete decompressible gzip output.
        use std::io::Read;
        use std::io::Write;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("out.fq.gz");
        let mut sink = FastxSink::create(&p, 1, 5).unwrap();
        sink.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        sink.finalize().unwrap();

        let f = std::fs::File::open(&p).unwrap();
        let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
        let mut out = String::new();
        dec.read_to_string(&mut out).unwrap();
        assert_eq!(out, "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_filter_fasta_streaming_handles_missing_read() {
        // A FASTA record absent from the kraken stream resolves to taxon 0
        // (rejected when not in the target set), matching the map-based path.
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let in_path = dir.path().join("in.fa");
        std::fs::write(&in_path, ">s1\nACGT\n>s2\nTTTT\n>s3\nGGGG\n").unwrap();
        let kraken_path = dir.path().join("kraken.tsv");
        let mut f = std::fs::File::create(&kraken_path).unwrap();
        // s2 omitted entirely from the stream.
        f.write_all(b"C\ts1\t9606\t4\t9606:1\nC\ts3\t9606\t4\t9606:1\n")
            .unwrap();

        let out_path = dir.path().join("out.fa");
        let mut taxa = ahash::AHashSet::new();
        taxa.insert(9606u32);

        super::run_filter(super::FilterArgs {
            input: in_path,
            output: out_path.clone(),
            taxon_ids: taxa,
            per_record: true,
            classifications: Some(kraken_path),
            ..super::FilterArgs::default_for_test()
        })
        .unwrap();

        let got = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(got, ">s1\nACGT\n>s3\nGGGG\n");
    }
}
