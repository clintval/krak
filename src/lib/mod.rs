//! `kraklib`; the Kraken toolkit you'll get addicted to.
#![warn(missing_docs)]

pub mod annotate;
pub mod filter;
pub mod kraken_report;
pub(crate) mod kraken_report_embed;
pub mod kraken_result;
pub mod kraken_taxonomy;
pub mod n2ref;
pub mod prep;
pub mod report2tsv;

pub use annotate::{run_annotate, AnnotateArgs};
pub use filter::{run_filter, FilterArgs};
pub use n2ref::{run_n2ref, N2RefArgs};
pub use prep::{run_prep, PrepArgs};
pub use report2tsv::{run_report2tsv, Report2TsvArgs};

/// SAM auxiliary tag for the Kraken taxonomic ID: `ti`.
pub(crate) const TI_TAG: [u8; 2] = [b't', b'i'];

/// Return `true` if the SAM/BAM/CRAM header indicates query-grouped ordering
/// (`SO:queryname` or `GO:query` in the `@HD` line).
pub(crate) fn is_query_grouped(header: &noodles::sam::Header) -> bool {
    use noodles::sam::header::record::value::map::header::{group_order, sort_order, tag};
    let Some(hd) = header.header() else {
        return false;
    };
    let fields = hd.other_fields();
    let is_query_sort = fields
        .get(&tag::SORT_ORDER)
        .map(|v| v.as_slice() == sort_order::QUERY_NAME)
        .unwrap_or(false);
    let is_query_group = fields
        .get(&tag::GROUP_ORDER)
        .map(|v| v.as_slice() == group_order::QUERY)
        .unwrap_or(false);
    is_query_sort || is_query_group
}

/// Open a BAM file backed by at most two BGZF decompression workers.
pub(crate) fn open_bam_reader(
    path: &std::path::Path,
) -> anyhow::Result<noodles::bam::io::Reader<noodles::bgzf::io::MultithreadedReader<std::fs::File>>>
{
    use anyhow::Context as _;
    let workers = std::num::NonZero::new(
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
            .min(2),
    )
    .unwrap();
    let file = std::fs::File::open(path)
        .with_context(|| format!("failed to open BAM: {}", path.display()))?;
    Ok(noodles::bam::io::Reader::from(
        noodles::bgzf::io::MultithreadedReader::with_worker_count(workers, file),
    ))
}

/// The format inferred by sniffing the first bytes of an input stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SniffedFormat {
    /// FASTA; first byte `>`.
    Fasta,
    /// FASTQ; first byte `@`, no tab in the first line.
    Fastq,
    /// SAM (text); first byte `@`, first line contains a tab.
    Sam,
    /// BAM; `BAM\x01` after BGZF decompression.
    Bam,
    /// CRAM; `CRAM` magic prefix.
    Cram,
    /// First few bytes don't match any known signature.
    Unknown,
}

/// Classify a slice from the start of an input stream (after any gzip layer
/// has been removed). Tolerates short slices and returns `Unknown` when it
/// cannot decide.
pub fn sniff_bytes(bytes: &[u8]) -> SniffedFormat {
    if bytes.starts_with(b"BAM\x01") {
        return SniffedFormat::Bam;
    }
    if bytes.starts_with(b"CRAM") {
        return SniffedFormat::Cram;
    }
    match bytes.first() {
        Some(b'>') => SniffedFormat::Fasta,
        Some(b'@') => sniff_at_prefixed(bytes),
        _ => SniffedFormat::Unknown,
    }
}

/// Helper: input begins with `@`. Decide between SAM (header line; contains a
/// tab in the first line) and FASTQ (read name; no tab in the first line).
fn sniff_at_prefixed(bytes: &[u8]) -> SniffedFormat {
    if let Some(eol) = memchr::memchr(b'\n', bytes) {
        let line = &bytes[..eol];
        if line.contains(&b'\t') {
            SniffedFormat::Sam
        } else {
            SniffedFormat::Fastq
        }
    } else {
        // No newline within the buffered chunk. Fall back to the
        // byte-position rule: SAM headers are always `@<2 uppercase>\t...`.
        if bytes.len() >= 4
            && bytes[1].is_ascii_uppercase()
            && bytes[2].is_ascii_uppercase()
            && bytes[3] == b'\t'
        {
            SniffedFormat::Sam
        } else {
            SniffedFormat::Fastq
        }
    }
}

/// Open `path` and identify the underlying format by peeking the first
/// bytes, transparently decoding a gzip layer for the *peek only*.
///
/// The returned `BufReader<File>` is left **at byte 0**; its internal
/// buffer holds the bytes that were peeked, so the caller can hand it
/// straight to `bgzf::Reader::new`, `MultiGzDecoder::new`, or a plain
/// `fastq::io::Reader::new` and start reading from the start of the file.
///
/// `gzipped` is `true` iff the outer stream starts with `1f 8b`. The caller
/// dispatches by `(format, gzipped)`:
///
/// - `Fasta`/`Fastq`/`Sam` + `gzipped=false` -> wrap in the noodles text reader directly.
/// - `Fasta`/`Fastq`/`Sam` + `gzipped=true`  -> wrap the reader in `MultiGzDecoder` first.
/// - `Bam` (gzipped flag ignored) -> `bam::Reader::from(bgzf::Reader::new(reader))`.
/// - `Cram` (gzipped must be `false`) -> `cram::Reader::new(reader)`.
///
/// Errors:
/// - `io::Error` from the underlying open / fill_buf.
/// - `InvalidData` when the user has gzip-wrapped a CRAM file (not a real format).
pub fn sniff_input(
    path: &std::path::Path,
) -> std::io::Result<(SniffedFormat, bool, std::io::BufReader<std::fs::File>)> {
    use std::io::{BufRead, BufReader, Read};

    // 64 KiB is large enough to hold any single BGZF block plus headroom.
    let file = std::fs::File::open(path)?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);

    let head: &[u8] = reader.fill_buf()?;
    let is_gzip = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;

    let format = if is_gzip {
        // Decompress a *copy* of the buffered bytes into a small scratch
        // buffer. The original BufReader is never `consume()`d, so its bytes
        // remain available to the caller. A decompression failure (truncated
        // member, bad CRC, etc.) is swallowed: we fall back to `Unknown`
        // rather than propagating, since the caller's downstream open will
        // surface a clearer "failed to read FASTQ/BAM" error with full path.
        let mut decompressed = [0u8; 256];
        let cursor = std::io::Cursor::new(head);
        let mut dec = flate2::bufread::MultiGzDecoder::new(cursor);
        let n = dec.read(&mut decompressed).unwrap_or(0);
        let inner = sniff_bytes(&decompressed[..n]);
        if inner == SniffedFormat::Cram {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "gzipped CRAM is not a real format; remove the .gz wrapping ({})",
                    path.display()
                ),
            ));
        }
        inner
    } else {
        sniff_bytes(head)
    };

    Ok((format, is_gzip, reader))
}

/// Strip a Kraken v1 `/1` or `/2` read-pair suffix from a read name, if present.
///
/// Kraken v1 appends `/1` and `/2` to distinguish mates when producing paired
/// output. Stripping these before comparison lets paired modes handle that
/// naming convention without erroring on mismatched names.
pub(crate) fn strip_pair_suffix(name: &str) -> &str {
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}

/// `true` when the path refers to a process pseudo-file (`/dev/stdin`,
/// `/dev/stdout`, `/dev/stderr`, `/dev/fd/N`, `/proc/self/fd/N`). Used to
/// skip behaviors that only make sense on real on-disk files (e.g. building
/// a sidecar index next to the output).
pub(crate) fn is_pseudo_path(path: &std::path::Path) -> bool {
    let s = path.to_str().unwrap_or("");
    s == "/dev/stdin"
        || s == "/dev/stdout"
        || s == "/dev/stderr"
        || s.starts_with("/dev/fd/")
        || s.starts_with("/proc/self/fd/")
}

/// `true` when the SAM header marks the file as `SO:coordinate`.
pub(crate) fn is_coordinate_sorted(header: &noodles::sam::Header) -> bool {
    use noodles::sam::header::record::value::map::header::{sort_order, tag};
    header
        .header()
        .and_then(|hd| hd.other_fields().get(&tag::SORT_ORDER))
        .map(|v| v.as_slice() == sort_order::COORDINATE)
        .unwrap_or(false)
}

/// Build a sibling index (`.bai` for BAM, `.crai` for CRAM) next to a
/// just-finalized alignment output, if and only if the output is a real
/// on-disk file in a format that supports indexing and the header is marked
/// `SO:coordinate`.
///
/// Skips silently for SAM (no native index format), for stdin/stdout/fd
/// pseudo-paths, and for non-coordinate-sorted outputs. Emits an `info!` log
/// before doing the work; the index is built by re-reading the file via
/// noodles' built-in `bam::fs::index` / `cram::fs::index` and written next to
/// the output as `<path>.bai` or `<path>.crai` (samtools convention).
pub(crate) fn maybe_index_alignment_output(
    path: &std::path::Path,
    header: &noodles::sam::Header,
    fmt: AlignmentFormat,
) -> anyhow::Result<()> {
    use anyhow::Context as _;

    if is_pseudo_path(path) {
        return Ok(());
    }
    if matches!(fmt, AlignmentFormat::Sam) {
        return Ok(());
    }
    if !is_coordinate_sorted(header) {
        return Ok(());
    }

    match fmt {
        AlignmentFormat::Bam => {
            let bai_path = sibling_index_path(path, "bai");
            log::info!("Indexing BAM: {} -> {}", path.display(), bai_path.display());
            let index = noodles::bam::fs::index(path)
                .with_context(|| format!("failed to build BAM index for {}", path.display()))?;
            noodles::bam::bai::fs::write(&bai_path, &index)
                .with_context(|| format!("failed to write BAM index: {}", bai_path.display()))?;
        }
        AlignmentFormat::Cram => {
            let crai_path = sibling_index_path(path, "crai");
            log::info!(
                "Indexing CRAM: {} -> {}",
                path.display(),
                crai_path.display()
            );
            let index = noodles::cram::fs::index(path)
                .with_context(|| format!("failed to build CRAM index for {}", path.display()))?;
            noodles::cram::crai::fs::write(&crai_path, &index)
                .with_context(|| format!("failed to write CRAM index: {}", crai_path.display()))?;
        }
        AlignmentFormat::Sam => unreachable!("guarded above"),
    }

    Ok(())
}

/// Build `<path>.<ext>` losslessly: appends to the OS string rather than
/// formatting through `Display`, so non-UTF-8 byte sequences in the path
/// survive intact.
fn sibling_index_path(path: &std::path::Path, ext: &str) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    std::path::PathBuf::from(s)
}

/// Build a `BufRead` for a sniffed text input. `gzipped=true` wraps the reader
/// in a `MultiGzDecoder`. The result is ready to hand to a noodles text reader
/// (`fasta::io::Reader::new`, `fastq::io::Reader::new`, `sam::io::Reader::new`).
///
/// The inner `BufReader<File>` already provides 64 KiB of buffering for the
/// gzip decoder's input side. The decoder's output is `Read`-only, so we still
/// need a `BufReader` wrapper to satisfy the noodles `BufRead` requirement;
/// but we use a small (8 KiB) buffer here to avoid duplicating the inner one.
pub fn into_text_bufread(
    reader: std::io::BufReader<std::fs::File>,
    gzipped: bool,
) -> Box<dyn std::io::BufRead> {
    if gzipped {
        Box::new(std::io::BufReader::with_capacity(
            8 * 1024,
            flate2::bufread::MultiGzDecoder::new(reader),
        ))
    } else {
        Box::new(reader)
    }
}

/// Open a FASTX file, transparently decompressing gzip and BGZF input.
///
/// Detection is byte-based: peeks the first two bytes for the gzip magic
/// number `1f 8b`, so this works for `/dev/stdin` and other paths without a
/// meaningful extension. `flate2::bufread::MultiGzDecoder` is used because it
/// reads concatenated gzip members, which transparently handles both plain
/// gzip and BGZF.
///
/// File-open errors are propagated raw so callers can wrap them with
/// format-specific context (e.g. `failed to open FASTQ: <path>`).
pub fn open_fastx_reader(path: &std::path::Path) -> std::io::Result<Box<dyn std::io::BufRead>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)?;
    let mut reader = BufReader::new(file);
    let buf = reader.fill_buf()?;
    let is_gzip = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;
    if is_gzip {
        let dec = flate2::bufread::MultiGzDecoder::new(reader);
        Ok(Box::new(BufReader::new(dec)))
    } else {
        Ok(Box::new(reader))
    }
}

/// Open a FASTX file for writing.
///
/// If the path's extension is `gz` (case-insensitive), the writer is wrapped
/// in a gzip encoder at default compression (level 6). Otherwise output is
/// plain. `/dev/stdout` has no extension and is therefore never auto-gzipped;
/// users pipe through `gzip` if they want that.
///
/// File-create errors are propagated raw so callers can wrap them with
/// format-specific context.
pub fn open_fastx_writer(path: &std::path::Path) -> std::io::Result<Box<dyn std::io::Write>> {
    use std::io::BufWriter;

    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    let file = std::fs::File::create(path)?;
    if is_gz {
        let enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        Ok(Box::new(BufWriter::new(enc)))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

/// Alignment file format inferred from file extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    /// Binary Alignment Map.
    Bam,
    /// CRAM compressed alignment.
    Cram,
    /// Sequence Alignment Map (text).
    Sam,
}

impl AlignmentFormat {
    /// Detect from a file path's extension (case-insensitive).
    ///
    /// Paths that represent file descriptors or standard streams
    /// (`/dev/stdin`, `/dev/stdout`, `/dev/fd/*`, `/proc/self/fd/*`) default
    /// to `Bam` because piped BAM is the conventional format for streaming
    /// alignment data between tools. All other unrecognised extensions fall
    /// back to `Sam`. This asymmetry matters for SAM-via-stdin: callers that
    /// accept stdin must complement this helper with a content-sniff (see
    /// `sniff_input`) and reroute to the SAM handler when the head bytes
    /// disagree with the path-based default.
    pub fn from_path(path: &std::path::Path) -> Self {
        match infer_format(path) {
            InferredFormat::Alignment(a) => a,
            InferredFormat::Fastx(_) => Self::Sam,
        }
    }
}

/// FASTX subkind: FASTA or FASTQ. Distinguishes `.fa`/`.fasta` from
/// `.fq`/`.fastq` (with optional `.gz`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastxKind {
    /// FASTA (`.fa` / `.fasta`, optionally `.gz`).
    Fasta,
    /// FASTQ (`.fq` / `.fastq`, optionally `.gz`).
    Fastq,
}

/// Broad input format inferred from a file extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InferredFormat {
    /// FASTA / FASTQ (`.fa` / `.fasta` / `.fq` / `.fastq`, optionally `.gz`).
    Fastx(FastxKind),
    /// SAM / BAM / CRAM (defaults to `Sam` for unrecognised extensions, with
    /// pipe/fd paths defaulting to `Bam` for the same reasons documented on
    /// `AlignmentFormat::from_path`).
    Alignment(AlignmentFormat),
}

/// Single source of truth for inferring an input format from a file path.
///
/// Strips a trailing `.gz` (case-insensitive) so paths like `reads.fq.gz`
/// resolve to `Fastx(Fastq)`. Falls back to `Alignment(Sam)` for
/// unrecognised extensions, with pipe/fd paths defaulting to `Alignment(Bam)`.
pub fn infer_format(path: &std::path::Path) -> InferredFormat {
    let stem = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
    let lower = stem.to_ascii_lowercase();
    let inner = lower.strip_suffix(".gz").unwrap_or(&lower);
    let ext = std::path::Path::new(inner)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");
    match ext {
        "fa" | "fasta" => InferredFormat::Fastx(FastxKind::Fasta),
        "fq" | "fastq" => InferredFormat::Fastx(FastxKind::Fastq),
        "bam" => InferredFormat::Alignment(AlignmentFormat::Bam),
        "cram" => InferredFormat::Alignment(AlignmentFormat::Cram),
        "sam" => InferredFormat::Alignment(AlignmentFormat::Sam),
        _ if is_pseudo_path(path) => InferredFormat::Alignment(AlignmentFormat::Bam),
        _ => InferredFormat::Alignment(AlignmentFormat::Sam),
    }
}

/// Build a FASTA repository from an optionally indexed reference path.
///
/// Pass `Some(ref_path)` to open an indexed FASTA (`.fai` must exist alongside it).
/// Pass `None` to obtain an empty (default) repository, suitable for CRAM files
/// with embedded references or when decoding tags only.
pub(crate) fn build_fasta_repo(
    reference: Option<&std::path::Path>,
) -> anyhow::Result<noodles::fasta::Repository> {
    use anyhow::Context as _;
    use noodles::fasta;
    use noodles::fasta::repository::adapters::IndexedReader as FastaIndexedAdapter;

    match reference {
        Some(ref_path) => {
            let adapter = fasta::io::indexed_reader::Builder::default()
                .build_from_path(ref_path)
                .with_context(|| {
                    format!(
                        "failed to open indexed (.fai present) FASTA: {}",
                        ref_path.display()
                    )
                })
                .map(FastaIndexedAdapter::new)?;
            Ok(fasta::Repository::new(adapter))
        }
        None => Ok(fasta::Repository::default()),
    }
}

/// Open a CRAM reader, optionally wiring in an external reference FASTA.
///
/// Pass `Some(ref_path)` for CRAM files whose sequences were encoded against
/// an external reference (the common case from `samtools view -T ref.fa -C`).
/// The reference must have a `.fai` index alongside it.
///
/// Pass `None` for CRAM files with embedded references or when reading tags
/// only (sequences will decode as N if the file requires an external reference
/// and none is supplied).
pub(crate) fn open_cram_reader(
    path: &std::path::Path,
    reference: Option<&std::path::Path>,
) -> anyhow::Result<noodles::cram::io::Reader<std::fs::File>> {
    use anyhow::Context as _;

    let repo = build_fasta_repo(reference)?;
    noodles::cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_path(path)
        .with_context(|| format!("failed to open CRAM: {}", path.display()))
}

/// Bail with a clear error when a CRAM header declares reference sequences
/// but the caller did not provide `--cram-reference`. Without this, noodles
/// panics deep inside its slice decoder with "invalid slice reference
/// sequence name".
pub(crate) fn require_cram_reference_if_mapped(
    header: &noodles::sam::Header,
    reference: Option<&std::path::Path>,
) -> anyhow::Result<()> {
    if reference.is_none() && !header.reference_sequences().is_empty() {
        anyhow::bail!(
            "CRAM input has reference sequences in its header; \
             pass --cram-reference <path/to/ref.fa> (with a .fai index) \
             to decode reference-compressed records"
        );
    }
    Ok(())
}

/// Create a CRAM writer, optionally wiring in an external reference FASTA.
///
/// The noodles CRAM writer requires a reference repository to encode mapped
/// (aligned) records using reference-based compression. Pass `Some(ref_path)`
/// whenever the output may contain mapped reads; the reference must have a
/// `.fai` index alongside it.
///
/// Pass `None` only when writing fully unmapped records (e.g. test fixtures
/// with no alignment data).
pub(crate) fn open_cram_writer(
    path: &std::path::Path,
    reference: Option<&std::path::Path>,
) -> anyhow::Result<noodles::cram::io::Writer<std::io::BufWriter<std::fs::File>>> {
    use anyhow::Context as _;

    let repo = build_fasta_repo(reference)?;
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create CRAM: {}", path.display()))?;
    Ok(noodles::cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_writer(std::io::BufWriter::new(file)))
}

#[cfg(test)]
pub(crate) fn write_minimal_taxonomy_dmp(db_path: &std::path::Path) {
    let tax_dir = db_path.join("taxonomy");
    std::fs::create_dir_all(&tax_dir).unwrap();

    let nodes = concat!(
        "1\t|\t1\t|\tno rank\t|\n",
        "9989\t|\t1\t|\torder\t|\n",
        "10116\t|\t9989\t|\tspecies\t|\n",
    );
    std::fs::write(tax_dir.join("nodes.dmp"), nodes).unwrap();

    let names = concat!(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n",
        "9989\t|\tRodentia\t|\t\t|\tscientific name\t|\n",
        "10116\t|\tRattus norvegicus\t|\t\t|\tscientific name\t|\n",
    );
    std::fs::write(tax_dir.join("names.dmp"), names).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// Write a minimal, zero-record CRAM to a temp file and return the path.
    fn write_empty_cram(path: &std::path::Path) {
        use noodles::sam;
        let mut w = open_cram_writer(path, None).unwrap();
        w.write_header(&sam::Header::default()).unwrap();
        w.try_finish(&sam::Header::default()).unwrap();
    }

    /// Write a minimal FASTA + .fai to a temp dir. Returns fa_path.
    ///
    /// FASTA content:  >chr1\nACGT\n
    /// FAI:            chr1\t4\t6\t4\t5\n
    fn write_indexed_fasta(dir: &std::path::Path) -> std::path::PathBuf {
        let fa_path = dir.join("ref.fa");
        let fai_path = dir.join("ref.fa.fai");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();
        // name  len  offset  bases_per_line  bytes_per_line
        std::fs::write(&fai_path, b"chr1\t4\t6\t4\t5\n").unwrap();
        fa_path
    }

    #[test]
    fn test_open_cram_reader_no_reference() {
        let dir = tempfile::TempDir::new().unwrap();
        let cram_path = dir.path().join("test.cram");
        write_empty_cram(&cram_path);

        let mut reader = open_cram_reader(&cram_path, None).unwrap();
        reader.read_header().unwrap(); // must not panic or error
    }

    #[test]
    fn test_open_cram_reader_with_reference() {
        let dir = tempfile::TempDir::new().unwrap();
        let cram_path = dir.path().join("test.cram");
        write_empty_cram(&cram_path);
        let fa_path = write_indexed_fasta(dir.path());

        let mut reader = open_cram_reader(&cram_path, Some(&fa_path)).unwrap();
        reader.read_header().unwrap();
    }

    #[test]
    fn test_open_cram_reader_missing_file() {
        let result = open_cram_reader(std::path::Path::new("/nonexistent/path/to.cram"), None);
        assert!(result.is_err());
        let err = result.err().unwrap();
        let msg = format!("{:#}", err);
        assert!(msg.contains("failed to open CRAM"));
        assert!(msg.contains("/nonexistent/path/to.cram"));
    }

    #[test]
    fn test_require_cram_reference_if_mapped() {
        use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let mut header = noodles::sam::Header::default();
        // Empty header → no reference required.
        require_cram_reference_if_mapped(&header, None).unwrap();

        let len = NonZeroUsize::new(10).unwrap();
        header
            .reference_sequences_mut()
            .insert("chr1".as_bytes().into(), Map::<ReferenceSequence>::new(len));
        // With @SQ entries → must error when no reference supplied.
        let err = require_cram_reference_if_mapped(&header, None).unwrap_err();
        assert!(format!("{err:#}").contains("--cram-reference"));
        // With reference supplied → ok.
        require_cram_reference_if_mapped(&header, Some(std::path::Path::new("ref.fa"))).unwrap();
    }

    #[test]
    fn test_open_cram_reader_missing_fai() {
        let dir = tempfile::TempDir::new().unwrap();
        let cram_path = dir.path().join("test.cram");
        write_empty_cram(&cram_path);
        // Write FASTA but deliberately omit the .fai
        let fa_path = dir.path().join("ref.fa");
        std::fs::write(&fa_path, b">chr1\nACGT\n").unwrap();

        let result = open_cram_reader(&cram_path, Some(&fa_path));
        assert!(result.is_err());
        let msg = format!("{:#}", result.err().unwrap());
        assert!(msg.contains("failed to open indexed (.fai present) FASTA"));
    }

    #[test]
    fn test_open_fastx_reader_plain_passthrough() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("plain.fq");
        std::fs::write(&p, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let mut r = super::open_fastx_reader(&p).unwrap();
        let mut got = Vec::new();
        r.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_open_fastx_reader_decompresses_gzip() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("data.fq.gz");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            enc.finish().unwrap();
        }
        let mut r = super::open_fastx_reader(&p).unwrap();
        let mut got = Vec::new();
        r.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_open_fastx_reader_decompresses_bgzf() {
        use noodles::bgzf;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("data.fq.gz");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut w = bgzf::io::Writer::new(f);
            w.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            w.finish().unwrap();
        }
        let mut r = super::open_fastx_reader(&p).unwrap();
        let mut got = Vec::new();
        r.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_open_fastx_reader_short_file_is_plain() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("tiny");
        std::fs::write(&p, b"X").unwrap();
        let mut r = super::open_fastx_reader(&p).unwrap();
        let mut got = Vec::new();
        r.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"X");
    }

    #[test]
    fn test_open_fastx_writer_plain() {
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("out.fq");
        {
            let mut w = super::open_fastx_writer(&p).unwrap();
            w.write_all(b"hello\n").unwrap();
            w.flush().unwrap();
        }
        let got = std::fs::read(&p).unwrap();
        assert_eq!(got, b"hello\n");
    }

    #[test]
    fn test_open_fastx_writer_gzip_roundtrip() {
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("out.fq.gz");
        {
            let mut w = super::open_fastx_writer(&p).unwrap();
            w.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            w.flush().unwrap();
        }
        let f = std::fs::File::open(&p).unwrap();
        let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
        let mut got = Vec::new();
        dec.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn test_sniff_bytes_fasta() {
        assert_eq!(
            super::sniff_bytes(b">read1\nACGT\n"),
            super::SniffedFormat::Fasta
        );
    }

    #[test]
    fn test_sniff_bytes_fastq() {
        assert_eq!(
            super::sniff_bytes(b"@read1\nACGT\n+\n!!!!\n"),
            super::SniffedFormat::Fastq
        );
    }

    #[test]
    fn test_sniff_bytes_fastq_with_long_name_and_description() {
        assert_eq!(
            super::sniff_bytes(
                b"@HWI-D00179:5:H21RAADXX:1:1101:1110:2173 1:N:0:1\nACGT\n+\n!!!!\n"
            ),
            super::SniffedFormat::Fastq
        );
    }

    #[test]
    fn test_sniff_bytes_sam_header_lines() {
        assert_eq!(
            super::sniff_bytes(b"@HD\tVN:1.6\n"),
            super::SniffedFormat::Sam
        );
        assert_eq!(
            super::sniff_bytes(b"@SQ\tSN:chr1\tLN:100\n"),
            super::SniffedFormat::Sam
        );
        assert_eq!(
            super::sniff_bytes(b"@PG\tID:foo\n"),
            super::SniffedFormat::Sam
        );
        assert_eq!(
            super::sniff_bytes(b"@RG\tID:rg1\n"),
            super::SniffedFormat::Sam
        );
        assert_eq!(
            super::sniff_bytes(b"@CO\tcomment text\n"),
            super::SniffedFormat::Sam
        );
    }

    #[test]
    fn test_sniff_bytes_sam_with_unconventional_tag() {
        // First-line-tab rule classifies any `@…\t…\n` as SAM, robust to tag
        // case or non-standard tags.
        assert_eq!(
            super::sniff_bytes(b"@xy\tfoo:bar\n"),
            super::SniffedFormat::Sam
        );
    }

    #[test]
    fn test_sniff_bytes_bam() {
        assert_eq!(
            super::sniff_bytes(b"BAM\x01anything"),
            super::SniffedFormat::Bam
        );
    }

    #[test]
    fn test_sniff_bytes_cram() {
        assert_eq!(
            super::sniff_bytes(b"CRAM\x03\x00rest"),
            super::SniffedFormat::Cram
        );
    }

    #[test]
    fn test_sniff_bytes_at_with_no_newline_uses_position_fallback() {
        assert_eq!(
            super::sniff_bytes(b"@HD\tVN:1.6"),
            super::SniffedFormat::Sam
        );
        assert_eq!(
            super::sniff_bytes(b"@somelongnamewithoutanewline"),
            super::SniffedFormat::Fastq
        );
    }

    #[test]
    fn test_sniff_bytes_unknown() {
        assert_eq!(super::sniff_bytes(b""), super::SniffedFormat::Unknown);
        assert_eq!(
            super::sniff_bytes(b"random bytes"),
            super::SniffedFormat::Unknown
        );
    }

    #[test]
    fn test_sniff_input_plain_fasta_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x.fa");
        std::fs::write(&p, b">r1\nACGT\n>r2\nTTTT\n").unwrap();
        let (fmt, gz, mut r) = super::sniff_input(&p).unwrap();
        assert_eq!(fmt, super::SniffedFormat::Fasta);
        assert!(!gz);
        let mut got = Vec::new();
        r.read_to_end(&mut got).unwrap();
        assert_eq!(got, b">r1\nACGT\n>r2\nTTTT\n");
    }

    #[test]
    fn test_sniff_input_gzipped_fastq_no_extension() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\n!!!!\n").unwrap();
            enc.finish().unwrap();
        }
        let (fmt, gz, mut r) = super::sniff_input(&p).unwrap();
        assert_eq!(fmt, super::SniffedFormat::Fastq);
        assert!(gz);
        // Reader is at byte 0 of the raw gzip stream; re-decompressing must
        // reproduce the original bytes.
        let mut decoded = Vec::new();
        flate2::bufread::MultiGzDecoder::new(&mut r)
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, b"@r1\nACGT\n+\n!!!!\n");
    }

    #[test]
    fn test_sniff_input_bam() {
        use noodles::bgzf;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x.bam");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut w = bgzf::io::Writer::new(f);
            w.write_all(b"BAM\x01<rest>").unwrap();
            w.finish().unwrap();
        }
        let (fmt, gz, mut r) = super::sniff_input(&p).unwrap();
        assert_eq!(fmt, super::SniffedFormat::Bam);
        assert!(gz);
        let mut head = [0u8; 2];
        r.read_exact(&mut head).unwrap();
        assert_eq!(head, [0x1f, 0x8b]);
    }

    #[test]
    fn test_sniff_input_plain_cram_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x.cram");
        std::fs::write(&p, b"CRAM\x03\x00rest").unwrap();
        let (fmt, gz, _) = super::sniff_input(&p).unwrap();
        assert_eq!(fmt, super::SniffedFormat::Cram);
        assert!(!gz);
    }

    #[test]
    fn test_sniff_input_gzipped_cram_is_rejected() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"CRAM\x03\x00rest").unwrap();
            enc.finish().unwrap();
        }
        let err = super::sniff_input(&p).unwrap_err();
        assert!(err.to_string().contains("gzipped CRAM"));
    }

    #[test]
    fn test_into_text_bufread_plain_passthrough() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x.fa");
        std::fs::write(&p, b">r1\nACGT\n").unwrap();
        let (_, gz, reader) = super::sniff_input(&p).unwrap();
        let mut text = super::into_text_bufread(reader, gz);
        let mut got = Vec::new();
        text.read_to_end(&mut got).unwrap();
        assert_eq!(got, b">r1\nACGT\n");
    }

    #[test]
    fn test_into_text_bufread_gzipped_decodes() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("x");
        {
            let f = std::fs::File::create(&p).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\n!!!!\n").unwrap();
            enc.finish().unwrap();
        }
        let (_, gz, reader) = super::sniff_input(&p).unwrap();
        let mut text = super::into_text_bufread(reader, gz);
        let mut got = Vec::new();
        text.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"@r1\nACGT\n+\n!!!!\n");
    }

    #[test]
    fn test_sniff_input_empty_file_is_unknown() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("empty");
        std::fs::write(&p, b"").unwrap();
        let (fmt, _, _) = super::sniff_input(&p).unwrap();
        assert_eq!(fmt, super::SniffedFormat::Unknown);
    }

    #[test]
    fn test_open_fastx_writer_extension_case_insensitive() {
        use std::io::Write as _;
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("out.fq.GZ");
        {
            let mut w = super::open_fastx_writer(&p).unwrap();
            w.write_all(b"x").unwrap();
            w.flush().unwrap();
        }
        let f = std::fs::File::open(&p).unwrap();
        let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
        let mut got = Vec::new();
        dec.read_to_end(&mut got).unwrap();
        assert_eq!(got, b"x");
    }

    #[test]
    fn test_is_pseudo_path() {
        assert!(super::is_pseudo_path(std::path::Path::new("/dev/stdin")));
        assert!(super::is_pseudo_path(std::path::Path::new("/dev/stdout")));
        assert!(super::is_pseudo_path(std::path::Path::new("/dev/stderr")));
        assert!(super::is_pseudo_path(std::path::Path::new("/dev/fd/1")));
        assert!(super::is_pseudo_path(std::path::Path::new(
            "/proc/self/fd/0"
        )));
        assert!(!super::is_pseudo_path(std::path::Path::new(
            "/tmp/output.bam"
        )));
        assert!(!super::is_pseudo_path(std::path::Path::new("output.bam")));
    }

    #[test]
    fn test_is_coordinate_sorted() {
        use noodles::sam;
        let plain = sam::Header::default();
        assert!(!super::is_coordinate_sorted(&plain));

        let header_co: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().unwrap();
        assert!(super::is_coordinate_sorted(&header_co));

        let header_qn: sam::Header = "@HD\tVN:1.6\tSO:queryname\n".parse().unwrap();
        assert!(!super::is_coordinate_sorted(&header_qn));
    }

    #[test]
    fn test_maybe_index_alignment_output_skips_pseudo_path() {
        // /dev/stdout doesn't exist as a regular file; helper must skip
        // without ever reading or attempting to write a sidecar.
        let header: noodles::sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().unwrap();
        let p = std::path::Path::new("/dev/stdout");
        super::maybe_index_alignment_output(p, &header, super::AlignmentFormat::Bam).unwrap();
    }

    #[test]
    fn test_maybe_index_alignment_output_skips_sam() {
        // SAM has no native index format; helper must skip even if the
        // header is coordinate-sorted.
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("out.sam");
        std::fs::write(&p, b"@HD\tVN:1.6\tSO:coordinate\n").unwrap();
        let header: noodles::sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().unwrap();
        super::maybe_index_alignment_output(&p, &header, super::AlignmentFormat::Sam).unwrap();
        assert!(!dir.path().join("out.sam.bai").exists());
        assert!(!dir.path().join("out.sam.crai").exists());
    }

    #[test]
    fn test_maybe_index_alignment_output_skips_non_coordinate_sorted() {
        // Queryname-sorted BAM must NOT produce an index.
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("queryname.bam");
        let header: sam::Header = "@HD\tVN:1.6\tSO:queryname\n".parse().unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&p).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("r1".as_bytes().into());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        super::maybe_index_alignment_output(&p, &header, super::AlignmentFormat::Bam).unwrap();
        assert!(!dir.path().join("queryname.bam.bai").exists());
    }

    #[test]
    fn test_maybe_index_alignment_output_writes_bai_for_coord_sorted_bam() {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("coord.bam");
        // Header parsed from text: SO:coordinate + one @SQ so we have a
        // reference for the index.
        let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:8\n"
            .parse()
            .unwrap();
        {
            let mut w = bam::io::writer::Builder.build_from_path(&p).unwrap();
            w.write_header(&header).unwrap();
            let mut r = RecordBuf::default();
            *r.name_mut() = Some("r1".as_bytes().into());
            *r.flags_mut() = noodles::sam::alignment::record::Flags::default();
            *r.reference_sequence_id_mut() = Some(0);
            *r.alignment_start_mut() = Some(noodles::core::Position::MIN);
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            w.write_alignment_record(&header, &r).unwrap();
        }
        super::maybe_index_alignment_output(&p, &header, super::AlignmentFormat::Bam).unwrap();
        let bai = dir.path().join("coord.bam.bai");
        assert!(bai.exists(), "expected sibling .bai at {}", bai.display());
        // Sanity: index must be readable.
        noodles::bam::bai::fs::read(&bai).expect("written .bai must be readable");
    }

    #[test]
    fn test_maybe_index_alignment_output_writes_crai_for_coord_sorted_cram() {
        // Counterpart to the BAM/.bai test for the CRAM/.crai branch.
        // CRAM with no @SQ avoids the reference-resolution requirement;
        // SO:coordinate is what gates the indexer call.
        use noodles::sam;

        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path().join("coord.cram");
        let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().unwrap();
        {
            let mut w = super::open_cram_writer(&p, None).unwrap();
            w.write_header(&header).unwrap();
            w.try_finish(&header).unwrap();
        }
        super::maybe_index_alignment_output(&p, &header, super::AlignmentFormat::Cram).unwrap();
        let crai = dir.path().join("coord.cram.crai");
        assert!(
            crai.exists(),
            "expected sibling .crai at {}",
            crai.display()
        );
    }
}
