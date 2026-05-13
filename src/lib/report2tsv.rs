//! Convert a Kraken report to a flat TSV.

use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::kraken_report::KrakenReportEntry;

/// Arguments for the `report2tsv` subcommand.
///
/// Paths are passed through verbatim; `-` is resolved to `/dev/stdin` /
/// `/dev/stdout` by the CLI front-end before constructing this struct.
pub struct Report2TsvArgs {
    /// Input Kraken report file. Use `/dev/stdin` for stdin.
    pub input: PathBuf,
    /// Output TSV file. Use `/dev/stdout` for stdout.
    pub output: PathBuf,
}

/// TSV column header names, in output order.
const TSV_HEADER: &[&str] = &[
    "pct_fragments",
    "num_fragments_clade",
    "num_fragments_direct",
    "rank_code",
    "taxon_id",
    "name",
];

/// Run the `report2tsv` subcommand.
pub fn run_report2tsv(args: Report2TsvArgs) -> Result<()> {
    let entries = KrakenReportEntry::read_file(&args.input)?;

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args.output)
        .with_context(|| format!("failed to create output: {}", args.output.display()))?;

    writer
        .write_record(TSV_HEADER)
        .context("failed to write TSV header")?;

    for entry in &entries {
        let pct = format!("{:.2}", entry.pct_fragments);
        let clade = entry.num_fragments_clade.to_string();
        let direct = entry.num_fragments_direct.to_string();
        let tid = entry.taxon_id.to_string();
        writer
            .write_record([&pct, &clade, &direct, &entry.rank_code, &tid, &entry.name])
            .context("failed to write TSV row")?;
    }

    writer.flush().context("failed to flush TSV output")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const REPORT_LINES: &str = "100.00\t2000\t0\tR\t1\troot\n 99.95\t1999\t50\tD\t2\t  Bacteria\n  0.05\t1\t1\tS\t9606\t    Homo sapiens\n";

    fn run_stock_report() -> Vec<String> {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        write!(input, "{REPORT_LINES}").unwrap();
        let output = tempfile::NamedTempFile::new().unwrap();
        run_report2tsv(Report2TsvArgs {
            input: input.path().to_path_buf(),
            output: output.path().to_path_buf(),
        })
        .unwrap();
        std::fs::read_to_string(output.path())
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect()
    }

    #[test]
    fn test_report2tsv_columns_and_header() {
        let lines = run_stock_report();
        assert_eq!(
            lines[0],
            "pct_fragments\tnum_fragments_clade\tnum_fragments_direct\trank_code\ttaxon_id\tname"
        );
        assert_eq!(lines.len(), 4, "header + 3 data rows");
    }

    #[test]
    fn test_report2tsv_name_is_trimmed() {
        let lines = run_stock_report();
        let fields: Vec<&str> = lines.last().unwrap().split('\t').collect();
        assert_eq!(fields[5], "Homo sapiens");
    }

    #[test]
    fn test_report2tsv_pct_is_decimal() {
        let lines = run_stock_report();
        // Use the second data row (Bacteria, 99.95%) which has a fractional part.
        let pct = lines[2].split('\t').next().unwrap();
        assert!(pct.contains('.'), "pct should be a decimal: {pct}");
        assert!(pct.parse::<f64>().is_ok(), "pct should parse as f64: {pct}");
    }

    #[test]
    fn test_report2tsv_row_values() {
        let lines = run_stock_report();
        let row: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(row[0], "99.95");
        assert_eq!(row[1], "1999");
        assert_eq!(row[2], "50");
        assert_eq!(row[3], "D");
        assert_eq!(row[4], "2");
        assert_eq!(row[5], "Bacteria");
    }
}
