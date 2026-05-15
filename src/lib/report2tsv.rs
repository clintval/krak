//! Convert a Kraken report to a flat TSV.

use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::kraken_report::{KrakenReportEntry, SPACES_PER_DEPTH};

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
    "tax_id",
    "name",
    "rank",
    "level",
    "parent_tax_id",
    "parent_rank",
    "clade_count",
    "direct_count",
    "descendant_count",
    "frac_clade",
    "frac_direct",
    "frac_descendant",
    "minimizer_count",
    "distinct_minimizer_count",
];

/// Sum of clade counts across all depth-0 rows (i.e. unclassified + root).
fn total_sequences(entries: &[KrakenReportEntry]) -> u64 {
    entries
        .iter()
        .filter(|e| e.indent == 0)
        .map(|e| e.num_fragments_clade)
        .sum()
}

/// For each entry, derive `(parent_tax_id, parent_rank)` strings using an
/// indent-decreasing stack. Rows at indent 0 (root, unclassified) yield empty
/// strings.
fn parent_fields(entries: &[KrakenReportEntry]) -> Vec<(String, String)> {
    let mut out = Vec::with_capacity(entries.len());
    let mut stack: Vec<(usize, u32, String)> = Vec::new();
    for entry in entries {
        while let Some(top) = stack.last() {
            if top.0 < entry.indent {
                break;
            }
            stack.pop();
        }
        let parent = match stack.last() {
            Some((_, tid, rank)) => (tid.to_string(), rank.clone()),
            None => (String::new(), String::new()),
        };
        out.push(parent);
        stack.push((entry.indent, entry.taxon_id, entry.rank_code.clone()));
    }
    out
}

/// Run the `report2tsv` subcommand.
pub fn run_report2tsv(args: Report2TsvArgs) -> Result<()> {
    let entries = KrakenReportEntry::read_file(&args.input)?;
    let parents = parent_fields(&entries);
    let total = total_sequences(&entries);

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args.output)
        .with_context(|| format!("failed to create output: {}", args.output.display()))?;

    writer
        .write_record(TSV_HEADER)
        .context("failed to write TSV header")?;

    for (entry, (parent_tax_id, parent_rank)) in entries.iter().zip(parents.iter()) {
        let descendant_count = entry.num_fragments_clade - entry.num_fragments_direct;
        let (frac_clade, frac_direct, frac_descendant) = if total > 0 {
            let t = total as f64;
            (
                entry.num_fragments_clade as f64 / t,
                entry.num_fragments_direct as f64 / t,
                descendant_count as f64 / t,
            )
        } else {
            (0.0, 0.0, 0.0)
        };
        let level = entry.indent / SPACES_PER_DEPTH;
        let minimizer = entry
            .minimizer_count
            .map_or(String::new(), |v| v.to_string());
        let distinct_minimizer = entry
            .distinct_minimizer_count
            .map_or(String::new(), |v| v.to_string());

        writer
            .write_record([
                &entry.taxon_id.to_string(),
                &entry.name,
                &entry.rank_code,
                &level.to_string(),
                parent_tax_id,
                parent_rank,
                &entry.num_fragments_clade.to_string(),
                &entry.num_fragments_direct.to_string(),
                &descendant_count.to_string(),
                &format!("{frac_clade}"),
                &format!("{frac_direct}"),
                &format!("{frac_descendant}"),
                &minimizer,
                &distinct_minimizer,
            ])
            .context("failed to write TSV row")?;
    }

    writer.flush().context("failed to flush TSV output")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Stock 6-column report: unclassified(50) + root(2000), Bacteria(1999) with
    /// direct=50, Homo sapiens(1) under Bacteria.
    const REPORT_6COL: &str = "  2.44\t50\t50\tU\t0\tunclassified\n\
                                97.56\t2000\t0\tR\t1\troot\n\
                                97.51\t1999\t50\tD\t2\t  Bacteria\n\
                                 0.05\t1\t1\tS\t9606\t    Homo sapiens\n";

    /// Stock 8-column report (extended minimizer-data layout): same structure
    /// as `REPORT_6COL` but with `mc`/`dmc` columns between direct_count and
    /// rank_code.
    const REPORT_8COL: &str = "  2.44\t50\t50\t0\t0\tU\t0\tunclassified\n\
                                97.56\t2000\t0\t8000\t4000\tR\t1\troot\n\
                                97.51\t1999\t50\t7800\t3900\tD\t2\t  Bacteria\n\
                                 0.05\t1\t1\t10\t5\tS\t9606\t    Homo sapiens\n";

    fn run_report(report_text: &str) -> Vec<Vec<String>> {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        write!(input, "{report_text}").unwrap();
        let output = tempfile::NamedTempFile::new().unwrap();
        run_report2tsv(Report2TsvArgs {
            input: input.path().to_path_buf(),
            output: output.path().to_path_buf(),
        })
        .unwrap();
        std::fs::read_to_string(output.path())
            .unwrap()
            .lines()
            .map(|line| line.split('\t').map(str::to_owned).collect())
            .collect()
    }

    #[test]
    fn test_6col_report_has_empty_minimizer_cols() {
        let rows = run_report(REPORT_6COL);
        for row in rows.iter().skip(1) {
            assert_eq!(
                row[12], "",
                "minimizer_count must be empty for 6-col reports"
            );
            assert_eq!(
                row[13], "",
                "distinct_minimizer_count must be empty for 6-col reports"
            );
        }
    }

    #[test]
    fn test_8col_report_carries_minimizer_values() {
        let rows = run_report(REPORT_8COL);
        // Header + 4 data rows.
        assert_eq!(rows.len(), 5);
        // unclassified row
        assert_eq!(rows[1][12], "0");
        assert_eq!(rows[1][13], "0");
        // root
        assert_eq!(rows[2][12], "8000");
        assert_eq!(rows[2][13], "4000");
        // leaf
        assert_eq!(rows[4][12], "10");
        assert_eq!(rows[4][13], "5");
    }

    #[test]
    fn test_fractions_are_unit_scale() {
        let rows = run_report(REPORT_6COL);
        for row in rows.iter().skip(1) {
            let frac_clade: f64 = row[9].parse().unwrap();
            assert!(
                (0.0..=1.0).contains(&frac_clade),
                "frac_clade out of [0,1]: {frac_clade}"
            );
        }
    }

    #[test]
    fn test_frac_clade_matches_count_over_total() {
        // total = unclassified.clade(50) + root.clade(2000) = 2050.
        let rows = run_report(REPORT_6COL);
        // Bacteria row has clade=1999, direct=50, descendant=1949.
        let bacteria = &rows[3];
        let total = 50.0 + 2000.0;
        let expected_clade: f64 = 1999.0 / total;
        let expected_direct: f64 = 50.0 / total;
        let expected_desc: f64 = 1949.0 / total;
        assert!((bacteria[9].parse::<f64>().unwrap() - expected_clade).abs() < 1e-12);
        assert!((bacteria[10].parse::<f64>().unwrap() - expected_direct).abs() < 1e-12);
        assert!((bacteria[11].parse::<f64>().unwrap() - expected_desc).abs() < 1e-12);
    }

    #[test]
    fn test_descendant_count_is_clade_minus_direct() {
        let rows = run_report(REPORT_6COL);
        for row in rows.iter().skip(1) {
            let clade: u64 = row[6].parse().unwrap();
            let direct: u64 = row[7].parse().unwrap();
            let descendant: u64 = row[8].parse().unwrap();
            assert_eq!(descendant, clade - direct, "row: {row:?}");
        }
    }

    #[test]
    fn test_level_is_indent_divided_by_spaces_per_depth() {
        let rows = run_report(REPORT_6COL);
        // unclassified (indent 0), root (indent 0), Bacteria (indent 2 -> level 1),
        // Homo sapiens (indent 4 -> level 2).
        assert_eq!(rows[1][3], "0"); // unclassified
        assert_eq!(rows[2][3], "0"); // root
        assert_eq!(rows[3][3], "1"); // Bacteria
        assert_eq!(rows[4][3], "2"); // Homo sapiens
    }

    #[test]
    fn test_root_has_empty_parent_fields() {
        let rows = run_report(REPORT_6COL);
        assert_eq!(rows[2][4], "");
        assert_eq!(rows[2][5], "");
    }

    #[test]
    fn test_unclassified_has_empty_parent_fields() {
        let rows = run_report(REPORT_6COL);
        // Unclassified is at indent 0; we must not treat it as another
        // depth-0 row's parent — its parent fields are empty.
        assert_eq!(rows[1][0], "0", "unclassified tax_id is 0");
        assert_eq!(rows[1][4], "");
        assert_eq!(rows[1][5], "");
    }

    #[test]
    fn test_descendants_have_parent_fields_filled() {
        let rows = run_report(REPORT_6COL);
        // Bacteria (taxid 2, depth 1) -> parent is root (taxid 1, rank R).
        assert_eq!(rows[3][4], "1");
        assert_eq!(rows[3][5], "R");
        // Homo sapiens (taxid 9606, depth 2) -> parent is Bacteria.
        assert_eq!(rows[4][4], "2");
        assert_eq!(rows[4][5], "D");
    }

    #[test]
    fn test_name_is_trimmed() {
        let rows = run_report(REPORT_6COL);
        assert_eq!(rows[4][1], "Homo sapiens");
    }

    #[test]
    fn test_zero_total_produces_zero_fractions() {
        // Single root row with zero counts -> total_sequences == 0 ->
        // fractions must be exactly 0.0 (not NaN, not divide-by-zero).
        let report = "0.00\t0\t0\tR\t1\troot\n";
        let rows = run_report(report);
        assert_eq!(rows[1][9], "0");
        assert_eq!(rows[1][10], "0");
        assert_eq!(rows[1][11], "0");
    }

    #[test]
    fn test_total_sequences_sums_depth0_clade() {
        // total = unclassified(50) + root(2000) = 2050.
        // frac_clade for root = 2000 / 2050.
        let rows = run_report(REPORT_6COL);
        let root_frac_clade: f64 = rows[2][9].parse().unwrap();
        assert!((root_frac_clade - (2000.0 / 2050.0)).abs() < 1e-12);
    }

    #[test]
    fn test_tax_id_and_rank_columns() {
        let rows = run_report(REPORT_6COL);
        assert_eq!(rows[2][0], "1"); // root tax_id
        assert_eq!(rows[2][2], "R"); // root rank
        assert_eq!(rows[3][0], "2"); // Bacteria tax_id
        assert_eq!(rows[3][2], "D"); // Bacteria rank
    }

    #[test]
    fn test_sibling_subtrees_get_correct_parents() {
        // Tree: root(1) -> A(2) -> A1(3); root(1) -> B(4) -> B1(5).
        // After visiting A1 (indent 4), the indent-stack must pop A1 *and* A
        // when we reach B (indent 2) so B's parent is root, not A.
        let report = "100.00\t10\t0\tR\t1\troot\n\
                       50.00\t5\t1\tD\t2\t  A\n\
                       40.00\t4\t4\tS\t3\t    A1\n\
                       50.00\t5\t1\tD\t4\t  B\n\
                       40.00\t4\t4\tS\t5\t    B1\n";
        let rows = run_report(report);
        // B (row index 4): parent is root(1, R).
        assert_eq!(rows[4][0], "4");
        assert_eq!(rows[4][4], "1");
        assert_eq!(rows[4][5], "R");
        // B1 (row index 5): parent is B(4, D).
        assert_eq!(rows[5][0], "5");
        assert_eq!(rows[5][4], "4");
        assert_eq!(rows[5][5], "D");
    }

    #[test]
    fn test_emits_only_header_for_empty_report() {
        let rows = run_report("");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0][0], "tax_id");
    }
}
