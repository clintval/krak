//! An addicting set of Kraken-enhancing tools.
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};
use env_logger::Env;
use log::*;

use clap::builder::styling::{AnsiColor, Effects, Style, Styles};

use kraklib::annotate::AnnotateArgs;
use kraklib::filter::FilterArgs;
use kraklib::n2ref::N2RefArgs;
use kraklib::prep::PrepArgs;
use kraklib::report2tsv::Report2TsvArgs;

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// Cargo's color style.
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// An addicting set of Kraken-enhancing tools.
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    color = clap::ColorChoice::Always,
    term_width = 80
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Convert FASTX/SAM/BAM/CRAM for Kraken classification. Single-end and
    /// interleaved FASTA/FASTQ are both accepted; interleaved is auto-detected
    /// from read names (`/1`/`/2` suffixes or matching mate names, including
    /// Casava 1.8+). Query-grouped SAM/BAM/CRAM is auto-detected from the
    /// `@HD` header. Pass `--per-record` to disable auto-detection.
    Prep(PrepCmd),
    /// Annotate SAM/BAM/CRAM records with Kraken classifications.
    Annotate(AnnotateCmd),
    /// Filter FASTX/SAM/BAM/CRAM records by Kraken classifications.
    Filter(FilterCmd),
    /// Revert aligned N-calls in SAM/BAM/CRAM to reference bases.
    #[command(name = "n2ref")]
    N2Ref(N2RefCmd),
    /// Convert a Kraken report to a flat TSV.
    #[command(name = "report2tsv")]
    Report2Tsv(Report2TsvCmd),
}

/// Arguments for the `prep` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct PrepCmd {
    /// Input file(s). For FASTA/FASTQ, one file (single-end or interleaved) or
    /// two files (R1, R2 paired-end). For SAM/BAM/CRAM, exactly one file.
    /// Accepted positionally or as -i/--input. Use `-` or omit for stdin
    /// (single-input only).
    #[arg(index = 1, value_name = "FILE", num_args = 0..=2)]
    input_positional: Vec<PathBuf>,

    /// Input file(s) (flag form; equivalent to the positional arguments).
    /// Accepts 1 file (single-end / interleaved FASTA/FASTQ or SAM/BAM/CRAM)
    /// or 2 files (paired-end FASTA/FASTQ R1, R2). Use `-` for stdin
    /// (single-input only).
    #[arg(
        short = 'i',
        long = "input",
        value_name = "FILE",
        num_args = 1..=2,
        conflicts_with = "input_positional"
    )]
    input_flag: Vec<PathBuf>,

    /// Disable auto pair-detection. Each FASTQ/FASTA record (or SAM/BAM/CRAM
    /// primary record) is emitted as its own single-end template, even when
    /// the input looks interleaved (`/1`/`/2` suffixes or matching mate names)
    /// or query-grouped. Mutually exclusive with two-file input. Secondary
    /// (0x100) and supplementary (0x800) alignments are always dropped,
    /// regardless of this flag.
    #[arg(long)]
    per_record: bool,

    /// Output FASTA file. Use `-` or omit for stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: PathBuf,

    /// Reference FASTA for CRAM decompression (requires a `.fai` index alongside
    /// it). Not needed for CRAM files with embedded references or for
    /// FASTX/SAM/BAM input.
    #[arg(long)]
    cram_reference: Option<PathBuf>,
}

/// Arguments for the `annotate` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct AnnotateCmd {
    /// Input SAM/BAM/CRAM file. Use `-` or omit for stdin.
    #[arg(short = 'i', long, default_value = "-")]
    input: PathBuf,

    /// Kraken classification output file (tab-delimited, 5 columns).
    /// Use `-` for stdin.
    #[arg(short = 'a', long)]
    assignments: PathBuf,

    /// Output SAM/BAM/CRAM file with `ti` tags added. Use `-` or omit for stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: PathBuf,

    /// Kraken report file. When provided, the taxonomy tree is embedded in the
    /// output header as a `@CO krak:report:` line, making `filter --kraken-report`
    /// unnecessary for downstream filtering. Mutually exclusive with --kraken-db.
    #[arg(short = 'R', long, conflicts_with = "kraken_db")]
    kraken_report: Option<PathBuf>,

    /// Kraken database directory. Reads the Kraken DB files in the directory
    /// and embeds the taxonomy tree in the output header, replacing the need
    /// for both `kraken2 --report` and `--kraken-report`. Mutually exclusive
    /// with --kraken-report.
    #[arg(short = 'd', long, conflicts_with = "kraken_report")]
    kraken_db: Option<PathBuf>,

    /// Load all assignments into memory before reading the input file. Use when
    /// the assignments file is substantially out of QNAME order relative to the
    /// input. By default, assignments are streamed record-by-record with a
    /// lookahead buffer: modest disorder (such as Kraken v1 multi-threaded
    /// output which flushes work-unit buffers in completion order) is handled
    /// automatically and the buffer grows only as deep as the actual disorder.
    /// Use `--unordered` when disorder is large or unpredictable (e.g. a
    /// completely unsorted file) which will load all assignments into memory
    /// upfront.
    #[arg(long)]
    unordered: bool,

    /// Reference FASTA for CRAM decompression (requires a `.fai` index alongside
    /// it). Not needed for CRAM files with embedded references or for SAM/BAM
    /// input.
    #[arg(long)]
    cram_reference: Option<PathBuf>,

    /// Number of bgzf compression worker threads for BAM output. At 1
    /// (default), one compressor + one writer thread pipeline with the
    /// annotation loop. Ignored for SAM (no compression) and CRAM (per-block
    /// codecs).
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    threads: u16,

    /// bgzf compression level (0-9) for BAM output. Ignored for SAM and CRAM.
    #[arg(long, default_value_t = 5, value_parser = clap::value_parser!(u8).range(0..=9))]
    compression_level: u8,
}

/// Arguments for the `filter` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct FilterCmd {
    /// Input file(s). For SAM/BAM/CRAM, one file annotated with `ti` tags
    /// (via `krak annotate`). For FASTA/FASTQ, one file (single-end /
    /// interleaved) or two files (R1, R2 paired-end); supply taxon IDs with
    /// `--classifications` (-c). Accepted positionally or as -i/--input.
    /// Use `-` or omit for stdin (single-input only).
    #[arg(index = 1, value_name = "FILE", num_args = 0..=2)]
    input_positional: Vec<PathBuf>,

    /// Input file(s) (flag form; equivalent to the positional arguments).
    /// Accepts 1 file (SAM/BAM/CRAM or single-end / interleaved FASTA/FASTQ)
    /// or 2 files (paired-end FASTA/FASTQ R1, R2). Use `-` for stdin
    /// (single-input only).
    #[arg(
        short = 'i',
        long = "input",
        value_name = "FILE",
        num_args = 1..=2,
        conflicts_with = "input_positional"
    )]
    input_flag: Vec<PathBuf>,

    /// Output file(s). Format matches input: SAM/BAM/CRAM for alignment input,
    /// FASTA/FASTQ for FASTX input. One file for single-input mode, two files
    /// for paired-end FASTA/FASTQ. Use `-` or omit for stdout (single-output
    /// only).
    #[arg(short = 'o', long = "output", value_name = "FILE", num_args = 1..=2)]
    output: Vec<PathBuf>,

    /// Kraken report file. For SAM/BAM/CRAM, serves as fallback when no
    /// taxonomy tree is embedded in the header (embed one via `krak annotate
    /// --kraken-report` or `--kraken-db`). For FASTA/FASTQ, required when
    /// using `--allow-ancestors` or `--include-descendants`.
    #[arg(short = 'R', long)]
    kraken_report: Option<PathBuf>,

    /// TSV metrics output file. If omitted, metrics are only logged.
    #[arg(short = 'm', long)]
    metrics: Option<PathBuf>,

    /// Taxon IDs to retain (repeat for multiple).
    #[arg(short = 't', long = "taxon-id", required = true)]
    taxon_ids: Vec<u32>,

    /// Output file(s) for rejected records. Format matches input: SAM/BAM/CRAM
    /// for alignment input, FASTA/FASTQ for FASTX input. One file for
    /// single-input mode, two files for paired-end FASTA/FASTQ.
    #[arg(short = 'r', long, value_name = "FILE", num_args = 1..=2)]
    rejects: Vec<PathBuf>,

    /// Also keep reads assigned to ancestors of target taxon IDs.
    #[arg(short = 'a', long)]
    allow_ancestors: bool,

    /// Maximum edit distance for rescuing off-taxa reads. Not applicable to
    /// FASTA/FASTQ input (requires MD tag and CIGAR from alignment).
    #[arg(long)]
    rescue_max_edit_distance: Option<u32>,

    /// Maximum number of indel events allowed in off-taxa rescue. Not
    /// applicable to FASTA/FASTQ input.
    #[arg(long)]
    rescue_max_indels: Option<u32>,

    /// Maximum length of any single indel in off-taxa rescue. Not applicable
    /// to FASTA/FASTQ input.
    #[arg(long)]
    rescue_max_indel_length: Option<u32>,

    /// For every COUNT 'N' bases in a read, reduce the rescue-max-edit-distance
    /// threshold by 1 (integer division: 4 Ns with `--rescue-n-adjustment 5`
    /// reduces by 0; 5 Ns reduces by 1; 10 Ns reduces by 2). Must be >= 1 if
    /// set. Not applicable to FASTA/FASTQ input.
    #[arg(long)]
    rescue_n_adjustment: Option<u32>,

    /// Disable auto template grouping. Each FASTQ/FASTA record (or SAM/BAM/CRAM
    /// primary record) is filtered as its own single-record template, even when
    /// the input looks interleaved (`/1`/`/2` suffixes or matching mate names)
    /// or query-grouped. Required when the input BAM/SAM/CRAM is
    /// coordinate-sorted (use `samtools sort -n` to sort by queryname instead).
    /// Secondary (0x100) and supplementary (0x800) alignments; which normally
    /// share their primary's keep/reject decision and whose `ti` tags are
    /// ignored for classification; are instead classified independently by
    /// their own `ti` tag under this mode.
    #[arg(long)]
    per_record: bool,

    /// Kraken2 per-read classification output file. Required with FASTA/FASTQ input;
    /// mutually exclusive with SAM/BAM/CRAM input.
    #[arg(short = 'c', long)]
    classifications: Option<PathBuf>,

    /// Also keep reads classified at any taxon in the clade of each target
    /// (expands the target set to all descendants). Requires a Kraken report
    /// via `--kraken-report`; for SAM/BAM/CRAM a report embedded in the header
    /// (via `krak annotate`) also suffices.
    #[arg(short = 'd', long)]
    include_descendants: bool,

    /// Also keep reads classified as unclassified (taxon ID 0).
    #[arg(short = 'u', long)]
    include_unclassified: bool,

    /// Reference FASTA for CRAM decompression (requires a `.fai` index alongside
    /// it). Not needed for CRAM files with embedded references, SAM/BAM input,
    /// or FASTA/FASTQ input.
    #[arg(long)]
    cram_reference: Option<PathBuf>,

    /// Keep records that lack a `ti` tag (unannotated reads). By default,
    /// records with no `ti` tag are rejected. Not applicable to FASTA/FASTQ
    /// input (taxon IDs come from `--classifications`, not from SAM tags).
    #[arg(long)]
    keep_unannotated: bool,

    /// Load all assignments into memory before reading the input file. Use
    /// when the `--classifications` file is substantially out of QNAME order
    /// relative to the input. By default, assignments are streamed
    /// record-by-record with a lookahead buffer: modest disorder (such as
    /// Kraken v1 multi-threaded output which flushes work-unit buffers in
    /// completion order) is handled automatically and the buffer grows only
    /// as deep as the actual disorder. Use `--unordered` when disorder is
    /// large or unpredictable (e.g. a completely unsorted file) which will
    /// load all assignments into memory upfront. Only applies to FASTA/FASTQ
    /// input; SAM/BAM/CRAM input reads taxon IDs from `ti` tags.
    #[arg(long)]
    unordered: bool,

    /// Number of bgzf compression worker threads for `.gz` outputs. At 1
    /// (default), one compressor + one writer thread pipeline with the main
    /// filter loop; higher values fan compression out across more workers.
    /// SAM and CRAM outputs ignore this value.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    threads: u16,

    /// bgzf compression level (0-9) for `.gz` outputs.
    #[arg(long, default_value_t = 5, value_parser = clap::value_parser!(u8).range(0..=9))]
    compression_level: u8,
}

/// Arguments for the `n2ref` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct N2RefCmd {
    /// Input SAM/BAM/CRAM file. Accepted positionally (first argument) or as
    /// -i. Use `-` or omit for stdin.
    #[arg(index = 1, value_name = "INPUT")]
    input_positional: Option<PathBuf>,

    /// Output SAM/BAM/CRAM file. Accepted positionally (second argument) or
    /// as -o. Use `-` or omit for stdout.
    #[arg(index = 2, value_name = "OUTPUT")]
    output_positional: Option<PathBuf>,

    /// Input SAM/BAM/CRAM file (flag form; equivalent to the first positional
    /// argument). Use `-` for stdin.
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    input_flag: Option<PathBuf>,

    /// Output SAM/BAM/CRAM file (flag form; equivalent to the second
    /// positional argument). Use `-` for stdout.
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output_flag: Option<PathBuf>,

    /// Reference FASTA file (must match SAM/BAM/CRAM reference dictionary).
    #[arg(short = 'r', long)]
    reference: PathBuf,

    /// Replacement base quality score for converted N-calls (0–93). Defaults to original quality.
    #[arg(short = 'q', long)]
    qual: Option<u8>,

    /// Number of bgzf compression worker threads for BAM output. At 1
    /// (default), one compressor + one writer thread pipeline with the
    /// n2ref loop. Ignored for SAM (no compression) and CRAM (per-block
    /// codecs).
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..))]
    threads: u16,

    /// bgzf compression level (0-9) for BAM output. Ignored for SAM and CRAM.
    #[arg(long, default_value_t = 5, value_parser = clap::value_parser!(u8).range(0..=9))]
    compression_level: u8,
}

/// Arguments for the `report2tsv` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct Report2TsvCmd {
    /// Input Kraken report file. Accepted positionally (first argument) or as
    /// -i. Use `-` or omit for stdin.
    #[arg(index = 1, value_name = "INPUT")]
    input_positional: Option<PathBuf>,

    /// Output TSV file. Accepted positionally (second argument) or as -o.
    /// Use `-` or omit for stdout.
    #[arg(index = 2, value_name = "OUTPUT")]
    output_positional: Option<PathBuf>,

    /// Input Kraken report file (flag form; equivalent to the first positional
    /// argument). Use `-` for stdin.
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    input_flag: Option<PathBuf>,

    /// Output TSV file (flag form; equivalent to the second positional
    /// argument). Use `-` for stdout.
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output_flag: Option<PathBuf>,
}

/// Replace `-` with the given pseudo-path (`/dev/stdin` or `/dev/stdout`).
fn resolve_dash(path: PathBuf, default: &'static str) -> PathBuf {
    if path.as_os_str() == "-" {
        PathBuf::from(default)
    } else {
        path
    }
}

fn resolve_input(path: PathBuf) -> PathBuf {
    resolve_dash(path, "/dev/stdin")
}

fn resolve_output(path: PathBuf) -> PathBuf {
    resolve_dash(path, "/dev/stdout")
}

/// Pick exactly one of a positional arg or its `--flag` equivalent. Defaults
/// to the given pseudo-path (`/dev/stdin` / `/dev/stdout`) if both are absent;
/// errors out the process if both are provided.
fn pick_one(
    positional: Option<PathBuf>,
    flag: Option<PathBuf>,
    flag_label: &str,
    default: &'static str,
) -> PathBuf {
    match (positional, flag) {
        (Some(p), None) | (None, Some(p)) => resolve_dash(p, default),
        (None, None) => PathBuf::from(default),
        (Some(_), Some(_)) => {
            error!("cannot specify {flag_label} both positionally and as a flag");
            process::exit(1);
        }
    }
}

/// Resolve a 1..=2 output list to a `(output, output2)` pair. Errors out if
/// `paired` and the list has 1 element (or vice versa). Defaults to a single
/// stdout pseudo-path when the list is empty.
fn pick_output_pair(output: Vec<PathBuf>, paired: bool) -> (PathBuf, Option<PathBuf>) {
    match (output.len(), paired) {
        (0, false) => (PathBuf::from("/dev/stdout"), None),
        (0, true) => {
            error!("paired input requires two output paths via -o/--output");
            process::exit(1);
        }
        (1, false) => (resolve_output(output.into_iter().next().unwrap()), None),
        (1, true) => {
            error!("paired input requires two output paths via -o/--output");
            process::exit(1);
        }
        (2, false) => {
            error!("two output paths supplied to -o/--output but input is single-end");
            process::exit(1);
        }
        (2, true) => {
            let mut it = output.into_iter();
            let o1 = resolve_output(it.next().unwrap());
            let o2 = it.next().unwrap();
            if o2.as_os_str() == "-" {
                error!("R2 output cannot be stdout (`-`)");
                process::exit(1);
            }
            (o1, Some(o2))
        }
        _ => unreachable!("clap enforces num_args = 1..=2"),
    }
}

/// Resolve a 0..=2 rejects list to `(rejects, rejects2)`. Same matching rules
/// as outputs, except an empty list means "no rejects file".
fn pick_rejects_pair(rejects: Vec<PathBuf>, paired: bool) -> (Option<PathBuf>, Option<PathBuf>) {
    match (rejects.len(), paired) {
        (0, _) => (None, None),
        (1, false) => (
            Some(resolve_output(rejects.into_iter().next().unwrap())),
            None,
        ),
        (1, true) => {
            error!("paired input requires two reject paths via -r/--rejects");
            process::exit(1);
        }
        (2, false) => {
            error!("two reject paths supplied to -r/--rejects but input is single-end");
            process::exit(1);
        }
        (2, true) => {
            let mut it = rejects.into_iter();
            let r1 = resolve_output(it.next().unwrap());
            let r2 = it.next().unwrap();
            if r2.as_os_str() == "-" {
                error!("R2 rejects cannot be stdout (`-`)");
                process::exit(1);
            }
            (Some(r1), Some(r2))
        }
        _ => unreachable!("clap enforces num_args = 1..=2"),
    }
}

/// Resolve a 1..=2 input list (positional or `--input`/`-i` flag form) to a
/// `(input, input2)` pair. `clap`'s `conflicts_with` guarantees the two
/// sources are mutually exclusive, so we pick whichever is non-empty.
/// Defaults to a single stdin pseudo-path when both are empty. Errors out
/// if `per_record` is set with two inputs.
fn pick_input_pair(
    positional: Vec<PathBuf>,
    flag: Vec<PathBuf>,
    flag_label: &str,
    per_record: bool,
) -> (PathBuf, Option<PathBuf>) {
    let files = if !positional.is_empty() {
        positional
    } else {
        flag
    };
    match files.len() {
        0 => (PathBuf::from("/dev/stdin"), None),
        1 => (resolve_input(files.into_iter().next().unwrap()), None),
        2 => {
            if per_record {
                error!(
                    "cannot use --per-record with two-file paired input \
                     (positional or {flag_label})"
                );
                process::exit(1);
            }
            let mut it = files.into_iter();
            let r1 = resolve_input(it.next().unwrap());
            let r2 = it.next().unwrap();
            if r2.as_os_str() == "-" {
                error!("R2 input cannot be stdin (`-`)");
                process::exit(1);
            }
            (r1, Some(r2))
        }
        _ => unreachable!("clap enforces num_args = 0..=2 / 1..=2"),
    }
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or("info");
    env_logger::Builder::from_env(env).init();

    let matches = Cli::command().term_width(80).get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    let result = match cli.command {
        Commands::Prep(cmd) => {
            let (input, input2) =
                pick_input_pair(cmd.input_positional, cmd.input_flag, "-i", cmd.per_record);
            kraklib::run_prep(PrepArgs {
                input,
                input2,
                per_record: cmd.per_record,
                output: resolve_output(cmd.output),
                cram_reference: cmd.cram_reference,
            })
        }
        Commands::Annotate(cmd) => kraklib::run_annotate(AnnotateArgs {
            input: resolve_input(cmd.input),
            assignments: resolve_input(cmd.assignments),
            output: resolve_output(cmd.output),
            kraken_report: cmd.kraken_report,
            kraken_db: cmd.kraken_db,
            unordered: cmd.unordered,
            cram_reference: cmd.cram_reference,
            threads: cmd.threads as usize,
            compression_level: cmd.compression_level as u32,
        }),
        Commands::Filter(cmd) => {
            let (input, input2) =
                pick_input_pair(cmd.input_positional, cmd.input_flag, "-i", cmd.per_record);
            let (output, output2) = pick_output_pair(cmd.output, input2.is_some());
            let (rejects, rejects2) = pick_rejects_pair(cmd.rejects, input2.is_some());
            kraklib::run_filter(FilterArgs {
                input,
                input2,
                output,
                output2,
                kraken_report: cmd.kraken_report,
                metrics: cmd.metrics,
                taxon_ids: cmd.taxon_ids.into_iter().collect(),
                rejects,
                rejects2,
                allow_ancestors: cmd.allow_ancestors,
                rescue_max_edit_distance: cmd.rescue_max_edit_distance,
                rescue_max_indels: cmd.rescue_max_indels,
                rescue_max_indel_length: cmd.rescue_max_indel_length,
                rescue_n_adjustment: cmd.rescue_n_adjustment,
                per_record: cmd.per_record,
                classifications: cmd.classifications,
                include_descendants: cmd.include_descendants,
                include_unclassified: cmd.include_unclassified,
                cram_reference: cmd.cram_reference,
                keep_unannotated: cmd.keep_unannotated,
                unordered: cmd.unordered,
                threads: cmd.threads as usize,
                compression_level: cmd.compression_level as u32,
            })
        }
        Commands::N2Ref(cmd) => {
            let input = pick_one(cmd.input_positional, cmd.input_flag, "-i", "/dev/stdin");
            let output = pick_one(cmd.output_positional, cmd.output_flag, "-o", "/dev/stdout");
            kraklib::run_n2ref(N2RefArgs {
                input,
                output,
                reference: cmd.reference,
                qual: cmd.qual,
                threads: cmd.threads as usize,
                compression_level: cmd.compression_level as u32,
            })
        }
        Commands::Report2Tsv(cmd) => {
            let input = pick_one(cmd.input_positional, cmd.input_flag, "-i", "/dev/stdin");
            let output = pick_one(cmd.output_positional, cmd.output_flag, "-o", "/dev/stdout");
            kraklib::run_report2tsv(Report2TsvArgs { input, output })
        }
    };

    match result {
        Ok(()) => process::exit(0),
        Err(e) => {
            error!("{e:#}");
            process::exit(1);
        }
    }
}
