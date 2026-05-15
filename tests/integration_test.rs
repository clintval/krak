use std::io::Write;
use std::path::Path;

use assert_cmd::Command;
use tempfile::NamedTempFile;

fn krak() -> Command {
    Command::cargo_bin("krak").unwrap()
}

fn write_tmp(content: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    write!(f, "{content}").unwrap();
    f
}

fn write_tmp_fastq(records: &[(&str, &str, &str)]) -> NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    for (name, seq, qual) in records {
        writeln!(f, "@{name}").unwrap();
        writeln!(f, "{seq}").unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{qual}").unwrap();
    }
    f
}

fn write_tmp_fasta(records: &[(&str, &str)]) -> NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".fasta")
        .tempfile()
        .unwrap();
    for (name, seq) in records {
        writeln!(f, ">{name}").unwrap();
        writeln!(f, "{seq}").unwrap();
    }
    f
}

fn write_tmp_sam(header: &[&str], records: &[&str]) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    for line in header {
        writeln!(f, "{line}").unwrap();
    }
    for line in records {
        writeln!(f, "{line}").unwrap();
    }
    f
}

fn write_tmp_classifications(entries: &[(&str, u32)]) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    for (name, taxon_id) in entries {
        let (code, kmer_col) = if *taxon_id == 0 {
            ('U', "0:100".to_owned())
        } else {
            ('C', format!("{taxon_id}:100"))
        };
        writeln!(f, "{code}\t{name}\t{taxon_id}\t100\t{kmer_col}").unwrap();
    }
    f
}

fn write_tmp_kraken_report(entries: &[(u32, usize, &str)]) -> NamedTempFile {
    // entries: (taxon_id, indent_spaces, name)
    let mut f = NamedTempFile::new().unwrap();
    for (taxon_id, level, name) in entries {
        let spaces = " ".repeat(*level);
        writeln!(f, "0.00\t0\t0\tS\t{taxon_id}\t{spaces}{name}").unwrap();
    }
    f
}

// Keep in sync with `kraklib::write_minimal_taxonomy_dmp` in src/lib/mod.rs.
pub fn write_minimal_taxonomy_dmp(db_path: &Path) {
    let tax_dir = db_path.join("taxonomy");
    std::fs::create_dir_all(&tax_dir).unwrap();

    // nodes.dmp: taxon_id | parent_id | rank | (remaining fields omitted)
    let nodes = concat!(
        "1\t|\t1\t|\tno rank\t|\n",
        "9989\t|\t1\t|\torder\t|\n",
        "10116\t|\t9989\t|\tspecies\t|\n",
    );
    std::fs::write(tax_dir.join("nodes.dmp"), nodes).unwrap();

    // names.dmp: taxon_id | name_txt | unique_name | name_class |
    let names = concat!(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n",
        "9989\t|\tRodentia\t|\t\t|\tscientific name\t|\n",
        "10116\t|\tRattus norvegicus\t|\t\t|\tscientific name\t|\n",
    );
    std::fs::write(tax_dir.join("names.dmp"), names).unwrap();
}

pub fn minimal_taxo_k2d_bytes() -> Vec<u8> {
    const MAGIC: &[u8; 8] = b"K2TAXDAT";

    // Name data:
    //   offset  0: "root\0"              (5 bytes)
    //   offset  5: "Rodentia\0"          (9 bytes)
    //   offset 14: "Rattus norvegicus\0" (18 bytes)
    let name_data: &[u8] = b"root\0Rodentia\0Rattus norvegicus\0";
    // Rank data:
    //   offset  0: "no rank\0"  (8 bytes)
    //   offset  8: "order\0"    (6 bytes)
    //   offset 14: "species\0"  (8 bytes)
    let rank_data: &[u8] = b"no rank\0order\0species\0";

    // (parent_id, first_child, child_count, name_offset, rank_offset, external_id, godparent_id)
    let node_specs: &[(u64, u64, u64, u64, u64, u64, u64)] = &[
        (0, 0, 0, 0, 0, 0, 0),       // Node 0: placeholder
        (0, 2, 1, 0, 0, 1, 0),       // Node 1: root (ext_id=1)
        (1, 3, 1, 5, 8, 9989, 0),    // Node 2: Rodentia (ext_id=9989)
        (2, 0, 0, 14, 14, 10116, 0), // Node 3: Rattus norvegicus (ext_id=10116)
    ];

    let mut buf: Vec<u8> = Vec::new();
    buf.extend_from_slice(MAGIC);
    buf.extend_from_slice(&(node_specs.len() as u64).to_le_bytes());
    buf.extend_from_slice(&(name_data.len() as u64).to_le_bytes());
    buf.extend_from_slice(&(rank_data.len() as u64).to_le_bytes());
    for &(p, fc, cc, no, ro, ei, gp) in node_specs {
        for v in [p, fc, cc, no, ro, ei, gp] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
    }
    buf.extend_from_slice(name_data);
    buf.extend_from_slice(rank_data);
    buf
}

#[test]
fn test_help_exits_zero() {
    krak().arg("--help").assert().success();
}

#[test]
fn test_version_exits_zero() {
    krak().arg("--version").assert().success();
}

#[test]
fn test_no_subcommand_fails() {
    krak().assert().failure();
}

#[test]
fn test_prep_help() {
    krak().args(["prep", "--help"]).assert().success();
}

#[test]
fn test_prep_single_end_fastq_positional() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">r1\nACGT\n"), "got: {out}");
    assert!(out.contains(">r2\nTTTT\n"), "got: {out}");
}

#[test]
fn test_prep_single_end_fastq_flag() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">r1\nACGT\n"), "got: {out}");
}

#[test]
fn test_prep_single_end_fasta() {
    let input = write_tmp_fasta(&[("s1", "ACGT"), ("s2", "GGCC")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">s1\nACGT\n"), "got: {out}");
    assert!(out.contains(">s2\nGGCC\n"), "got: {out}");
}

#[test]
fn test_prep_interleaved_fastq() {
    let input = write_tmp_fastq(&[("pair1", "AAAA", "IIII"), ("pair1", "TTTT", "IIII")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">pair1\nAAAANTTTT\n"), "got: {out}");
}

#[test]
fn test_prep_paired_fastq() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "IIII")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "JJJJ"), ("p2", "GGGG", "JJJJ")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">p1\nAAAANTTTT\n"), "got: {out}");
    assert!(out.contains(">p2\nCCCCNGGGG\n"), "got: {out}");
}

#[test]
fn test_prep_single_end_sam() {
    // No @HD SO/GO; treated as single-end.
    let header = ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:8"];
    let records = [
        "read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII",
        "read2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&header, &records);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">read1\nACGT\n"), "got: {out}");
}

#[test]
fn test_prep_single_end_sam_paired_flag_is_error() {
    // Flag 1 = SEGMENTED (paired in sequencing).
    let header = ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:8"];
    let records = ["read1\t1\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("read1"));
}

#[test]
fn test_prep_query_grouped_sam_auto() {
    // SO:queryname in @HD -> auto query-grouped mode.
    let header = ["@HD\tVN:1.6\tSO:queryname", "@SQ\tSN:chr1\tLN:8"];
    // Flag 65 = 0x41 = SEGMENTED | READ_1; flag 129 = 0x81 = SEGMENTED | READ_2.
    let records = [
        "pair1\t65\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII",
        "pair1\t129\tchr1\t1\t60\t4M\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&header, &records);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">pair1\nAAAANTTTT\n"), "got: {out}");
}

#[test]
fn test_prep_per_record_and_input2_together_is_error() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "JJJJ")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "--per-record",
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure();
}

#[test]
fn test_prep_input2_with_sam_is_error() {
    let header = ["@HD\tVN:1.6"];
    let input = write_tmp_sam(&header, &[]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "JJJJ")]);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            "-i",
            input.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .failure();
}

#[test]
fn test_annotate_help() {
    krak().args(["annotate", "--help"]).assert().success();
}

#[test]
fn test_annotate_with_report_embeds_co_line() {
    let header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606)]);
    let report = write_tmp_kraken_report(&[(1, 0, "root"), (9606, 2, "Homo sapiens")]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "--kraken-report",
            report.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.lines().any(|l| l.starts_with("@CO\tkrak:report:")),
        "expected @CO krak:report: line in header; got:\n{out}"
    );
}

#[test]
fn test_annotate_without_report_has_no_co_line() {
    let header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606)]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        !out.lines().any(|l| l.starts_with("@CO\tkrak:report:")),
        "expected no @CO krak:report: line without --kraken-report; got:\n{out}"
    );
}

#[test]
fn test_annotate_kraken_db_embeds_taxonomy_from_taxo_k2d() {
    use std::io::Write as _;

    // Write a synthetic taxo.k2d into a temp directory (simulating a Kraken DB).
    let db_dir = tempfile::tempdir().unwrap();
    let taxo_path = db_dir.path().join("taxo.k2d");
    let mut f = std::fs::File::create(&taxo_path).unwrap();
    f.write_all(&minimal_taxo_k2d_bytes()).unwrap();

    let header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read1", 10116)]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "--kraken-db",
            db_dir.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.lines().any(|l| l.starts_with("@CO\tkrak:report:")),
        "expected @CO krak:report: line from taxo.k2d; got:\n{out}"
    );
}

#[test]
fn test_annotate_kraken_db_v1_embeds_taxonomy_from_nodes_dmp() {
    // Simulate a Kraken v1 DB: taxonomy/nodes.dmp + names.dmp, but NO taxo.k2d.
    let db_dir = tempfile::tempdir().unwrap();
    write_minimal_taxonomy_dmp(db_dir.path());

    let header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read1", 10116)]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "--kraken-db",
            db_dir.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.lines().any(|l| l.starts_with("@CO\tkrak:report:")),
        "expected @CO krak:report: line from taxonomy/nodes.dmp; got:\n{out}"
    );
}

#[test]
fn test_annotate_kraken_db_and_kraken_report_together_is_error() {
    let db_dir = tempfile::tempdir().unwrap();
    let header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606)]);
    let report = write_tmp_kraken_report(&[(1, 0, "root")]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "--kraken-report",
            report.path().to_str().unwrap(),
            "--kraken-db",
            db_dir.path().to_str().unwrap(),
        ])
        .assert()
        .failure();
}

#[test]
fn test_annotate_default_mode_handles_out_of_order_assignments() {
    let header = ["@HD\tVN:1.6"];
    let records = [
        "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
        "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read2", 9606), ("read1", 10116)]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.lines()
            .any(|l| l.starts_with("read1\t") && l.contains("ti:i:10116")),
        "read1 should have ti:i:10116 on same line: {out}"
    );
    assert!(
        out.lines()
            .any(|l| l.starts_with("read2\t") && l.contains("ti:i:9606")),
        "read2 should have ti:i:9606 on same line: {out}"
    );
}

#[test]
fn test_annotate_unordered_flag_succeeds_with_out_of_order_assignments() {
    let header = ["@HD\tVN:1.6"];
    let records = [
        "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
        "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&header, &records);
    let assignments = write_tmp_classifications(&[("read2", 9606), ("read1", 10116)]);
    let output = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "--unordered",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.lines()
            .any(|l| l.starts_with("read1\t") && l.contains("ti:i:10116")),
        "read1 should have ti:i:10116 on same line: {out}"
    );
    assert!(
        out.lines()
            .any(|l| l.starts_with("read2\t") && l.contains("ti:i:9606")),
        "read2 should have ti:i:9606 on same line: {out}"
    );
}

#[test]
fn test_filter_help() {
    krak().args(["filter", "--help"]).assert().success();
}

#[test]
fn test_filter_requires_taxon_id() {
    let input = NamedTempFile::new().unwrap();
    let output = NamedTempFile::new().unwrap();
    let metrics = NamedTempFile::new().unwrap();
    krak()
        .args([
            "filter",
            "--input",
            input.path().to_str().unwrap(),
            "--output",
            output.path().to_str().unwrap(),
            "--metrics",
            metrics.path().to_str().unwrap(),
            // missing --taxon-id
        ])
        .assert()
        .failure();
}

#[test]
fn test_filter_allow_ancestors_requires_report() {
    let input = write_tmp("");
    let output = NamedTempFile::new().unwrap();
    let metrics = NamedTempFile::new().unwrap();
    krak()
        .args([
            "filter",
            "--input",
            input.path().to_str().unwrap(),
            "--output",
            output.path().to_str().unwrap(),
            "--metrics",
            metrics.path().to_str().unwrap(),
            "--taxon-id",
            "9606",
            "--allow-ancestors",
            // missing --kraken-report
        ])
        .assert()
        .failure();
}

#[test]
fn test_filter_single_end_fastq() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("r1", 9606), ("r2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.contains("@r1\n"),
        "kept read should be in output: {out}"
    );
    assert!(
        !out.contains("@r2\n"),
        "rejected read should not be in output: {out}"
    );
}

#[test]
fn test_filter_fastq_unordered_flag_succeeds_with_out_of_order_assignments() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("r2", 1234), ("r1", 9606)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "--unordered",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains("@r1\n"), "r1 should be kept: {out}");
    assert!(!out.contains("@r2\n"), "r2 should be rejected: {out}");
}

#[test]
fn test_filter_fastq_streaming_handles_modest_disorder() {
    let input = write_tmp_fastq(&[
        ("r1", "ACGT", "IIII"),
        ("r2", "TTTT", "JJJJ"),
        ("r3", "GGGG", "KKKK"),
    ]);
    let kraken_out = write_tmp_classifications(&[("r1", 9606), ("r3", 9606), ("r2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains("@r1\n"), "r1 should be kept: {out}");
    assert!(out.contains("@r3\n"), "r3 should be kept: {out}");
    assert!(!out.contains("@r2\n"), "r2 should be rejected: {out}");
}

#[test]
fn test_filter_single_end_fasta() {
    let input = write_tmp_fasta(&[("s1", "ACGT"), ("s2", "GGCC")]);
    let kraken_out = write_tmp_classifications(&[("s1", 9606), ("s2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fasta")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.contains(">s1\n"),
        "kept record should be in output: {out}"
    );
    assert!(
        !out.contains(">s2\n"),
        "rejected record should not be in output: {out}"
    );
}

#[test]
fn test_filter_fastq_rejects_file() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("r1", 9606), ("r2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let rejects = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-r",
            rejects.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let rejects_out = std::fs::read_to_string(rejects.path()).unwrap();
    assert!(
        rejects_out.contains("@r2\n"),
        "rejected read should be in rejects file: {rejects_out}"
    );
    assert!(
        !rejects_out.contains("@r1\n"),
        "kept read should not be in rejects: {rejects_out}"
    );
}

#[test]
fn test_filter_paired_fastq() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "JJJJ")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII"), ("p2", "GGGG", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606), ("p2", 1234)]);
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out1.path().to_str().unwrap(),
            out2.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out1_body = std::fs::read_to_string(out1.path()).unwrap();
    let out2_body = std::fs::read_to_string(out2.path()).unwrap();
    assert!(
        out1_body.contains("@p1\n"),
        "R1 kept missing p1: {out1_body}"
    );
    assert!(
        out2_body.contains("@p1\n"),
        "R2 kept missing p1: {out2_body}"
    );
    assert!(
        !out1_body.contains("@p2\n"),
        "R1 rejected p2 leaked: {out1_body}"
    );
    assert!(
        !out2_body.contains("@p2\n"),
        "R2 rejected p2 leaked: {out2_body}"
    );
}

#[test]
fn test_filter_paired_fastq_positional() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "JJJJ")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII"), ("p2", "GGGG", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606), ("p2", 1234)]);
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out1.path().to_str().unwrap(),
            out2.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let out1_body = std::fs::read_to_string(out1.path()).unwrap();
    let out2_body = std::fs::read_to_string(out2.path()).unwrap();
    assert!(out1_body.contains("@p1\n"));
    assert!(out2_body.contains("@p1\n"));
    assert!(!out1_body.contains("@p2\n"));
    assert!(!out2_body.contains("@p2\n"));
}

#[test]
fn test_filter_paired_fastq_with_rejects() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "JJJJ")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII"), ("p2", "GGGG", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606), ("p2", 1234)]);
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let rej1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let rej2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out1.path().to_str().unwrap(),
            out2.path().to_str().unwrap(),
            "-r",
            rej1.path().to_str().unwrap(),
            rej2.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .success();

    let rej1_body = std::fs::read_to_string(rej1.path()).unwrap();
    let rej2_body = std::fs::read_to_string(rej2.path()).unwrap();
    assert!(
        rej1_body.contains("@p2\n"),
        "R1 rejects should have p2: {rej1_body}"
    );
    assert!(
        rej2_body.contains("@p2\n"),
        "R2 rejects should have p2: {rej2_body}"
    );
    assert!(!rej1_body.contains("@p1\n"));
    assert!(!rej2_body.contains("@p1\n"));
}

#[test]
fn test_filter_paired_requires_two_outputs() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606)]);
    let out = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .failure();
}

#[test]
fn test_filter_paired_per_record_is_error() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606)]);
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out1.path().to_str().unwrap(),
            out2.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "--per-record",
            "-t",
            "9606",
        ])
        .assert()
        .failure();
}

#[test]
fn test_filter_paired_unequal_lengths_is_error() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII"), ("p2", "CCCC", "JJJJ")]);
    let r2 = write_tmp_fastq(&[("p1", "TTTT", "IIII")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606), ("p2", 1234)]);
    let out1 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let out2 = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            "-i",
            r1.path().to_str().unwrap(),
            r2.path().to_str().unwrap(),
            "-o",
            out1.path().to_str().unwrap(),
            out2.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("unequal record counts"));
}

#[test]
fn test_filter_positional_and_flag_conflict_fails() {
    let r1 = write_tmp_fastq(&[("p1", "AAAA", "IIII")]);
    let kraken_out = write_tmp_classifications(&[("p1", 9606)]);
    let out = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();

    krak()
        .args([
            "filter",
            r1.path().to_str().unwrap(),
            "-i",
            r1.path().to_str().unwrap(),
            "-o",
            out.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .failure();
}

#[test]
fn test_filter_fastq_classifications_required() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII")]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("--classifications"));
}

#[test]
fn test_filter_sam_classifications_is_error() {
    let header = ["@HD\tVN:1.6"];
    let input = write_tmp_sam(&header, &[]);
    let kraken_out = write_tmp_classifications(&[]);
    let output = NamedTempFile::new().unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("not valid for SAM/BAM/CRAM"));
}

#[test]
fn test_filter_include_unclassified() {
    // r1 is unclassified (taxon 0); with -u it should be kept.
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("r1", 0), ("r2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
            "--include-unclassified",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.contains("@r1\n"),
        "unclassified read should be kept: {out}"
    );
    assert!(
        !out.contains("@r2\n"),
        "off-taxa read should be rejected: {out}"
    );
}

#[test]
fn test_filter_include_descendants() {
    // Taxonomy: root(1) -> human(9606). Filter for taxon 1 + descendants -> keeps human reads.
    let report = write_tmp_kraken_report(&[(1, 0, "root"), (9606, 2, "human")]);
    let input = write_tmp_fastq(&[("read1", "ACGT", "IIII"), ("read2", "TTTT", "JJJJ")]);
    let kraken_out = write_tmp_classifications(&[("read1", 9606), ("read2", 1234)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-R",
            report.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "1",
            "--include-descendants",
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(
        out.contains("@read1\n"),
        "descendant read should be kept: {out}"
    );
    assert!(
        !out.contains("@read2\n"),
        "unrelated read should be rejected: {out}"
    );
}

#[test]
fn test_filter_include_descendants_requires_report() {
    let input = write_tmp_fastq(&[("r1", "ACGT", "IIII")]);
    let kraken_out = write_tmp_classifications(&[("r1", 9606)]);
    let output = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .unwrap();
    let metrics = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            input.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-c",
            kraken_out.path().to_str().unwrap(),
            "-m",
            metrics.path().to_str().unwrap(),
            "-t",
            "9606",
            "--include-descendants",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("taxonomy tree"));
}

#[test]
fn test_filter_uses_embedded_report_for_descendants() {
    // Annotate a SAM and embed the report in the header.
    let sam_header = ["@HD\tVN:1.6"];
    let records = [
        "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", // will get ti:9606
        "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII", // will get ti:1
        "read3\t4\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII", // will get ti:1234 (off-tree)
    ];
    let input = write_tmp_sam(&sam_header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606), ("read2", 1), ("read3", 1234)]);
    // Tree: root(1) -> Homo sapiens(9606). taxon 1234 is not in the tree.
    let report = write_tmp_kraken_report(&[(1, 0, "root"), (9606, 2, "Homo sapiens")]);
    let annotated = NamedTempFile::new().unwrap();

    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            annotated.path().to_str().unwrap(),
            "--kraken-report",
            report.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // Filter keeping taxon 1 and all descendants; no --kraken-report flag needed.
    let metrics = NamedTempFile::new().unwrap();
    let filtered = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            annotated.path().to_str().unwrap(),
            "-o",
            filtered.path().to_str().unwrap(),
            "-t",
            "1",
            "--include-descendants",
            "--per-record",
            "-m",
            metrics.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(filtered.path()).unwrap();
    assert!(
        out.contains("read1\t"),
        "read1 (descendant taxon 9606) should be kept: {out}"
    );
    assert!(
        out.contains("read2\t"),
        "read2 (target taxon 1) should be kept: {out}"
    );
    assert!(
        !out.contains("read3\t"),
        "read3 (taxon 1234, off-tree) should be filtered: {out}"
    );
}

#[test]
fn test_filter_hidden_kraken_report_flag_still_works() {
    // An annotated SAM without an embedded report; hidden --kraken-report fallback.
    let sam_header = ["@HD\tVN:1.6"];
    let records = [
        "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
        "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&sam_header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606), ("read2", 1234)]);
    let report = write_tmp_kraken_report(&[(1, 0, "root"), (9606, 2, "Homo sapiens")]);

    // Annotate WITHOUT embedding the report.
    let annotated = NamedTempFile::new().unwrap();
    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            annotated.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // Filter using the hidden --kraken-report fallback.
    let metrics = NamedTempFile::new().unwrap();
    let filtered = NamedTempFile::new().unwrap();

    krak()
        .args([
            "filter",
            "-i",
            annotated.path().to_str().unwrap(),
            "-o",
            filtered.path().to_str().unwrap(),
            "-t",
            "1",
            "--include-descendants",
            "--kraken-report",
            report.path().to_str().unwrap(),
            "--per-record",
            "-m",
            metrics.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let out = std::fs::read_to_string(filtered.path()).unwrap();
    assert!(
        out.contains("read1\t"),
        "read1 (descendant of 1) should be kept: {out}"
    );
    assert!(
        !out.contains("read2\t"),
        "read2 (taxon 1234, off-tree) should be filtered: {out}"
    );
}

#[test]
fn test_filter_include_descendants_without_tree_fails() {
    let sam_header = ["@HD\tVN:1.6"];
    let records = ["read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"];
    let input = write_tmp_sam(&sam_header, &records);
    let assignments = write_tmp_classifications(&[("read1", 9606)]);
    let annotated = NamedTempFile::new().unwrap();

    // Annotate without embedding the report.
    krak()
        .args([
            "annotate",
            "-i",
            input.path().to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            annotated.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    let metrics = NamedTempFile::new().unwrap();
    let output = NamedTempFile::new().unwrap();

    // Filter with --include-descendants but no tree source; must fail.
    krak()
        .args([
            "filter",
            "-i",
            annotated.path().to_str().unwrap(),
            "-o",
            output.path().to_str().unwrap(),
            "-t",
            "1",
            "--include-descendants",
            "-m",
            metrics.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("taxonomy tree"));
}

/// End-to-end-test, a full kraken2 pipeline.
///
/// Build a minimal Kraken2 database in a temp directory.
///
/// Taxonomy: root (1) -> Rodentia (9989) -> Rattus norvegicus (10116)
/// Library: one 65 bp sequence tagged as taxon 10116
///
/// Built with `kraken2-build --build --no-masking` so no dustmasker is needed.
/// Callers should keep the returned `TempDir` alive for the duration of the test.
fn build_minidb() -> tempfile::TempDir {
    let db = tempfile::tempdir().unwrap();
    let db_path = db.path();

    // Taxonomy files; minimal 3-field NCBI format
    std::fs::create_dir_all(db_path.join("taxonomy")).unwrap();
    std::fs::write(
        db_path.join("taxonomy/nodes.dmp"),
        "1\t|\t1\t|\tno rank\t|\n\
         9989\t|\t1\t|\torder\t|\n\
         10116\t|\t9989\t|\tspecies\t|\n",
    )
    .unwrap();
    std::fs::write(
        db_path.join("taxonomy/names.dmp"),
        "1\t|\troot\t|\t\t|\tscientific name\t|\n\
         9989\t|\tRodentia\t|\t\t|\tscientific name\t|\n\
         10116\t|\tRattus norvegicus\t|\t\t|\tscientific name\t|\n",
    )
    .unwrap();
    // Optional files that kraken2-build may open (empty is fine)
    std::fs::write(db_path.join("taxonomy/merged.dmp"), "").unwrap();
    std::fs::write(db_path.join("taxonomy/delnodes.dmp"), "").unwrap();

    // Library: one sequence > k=35 bp so at least one k-mer indexes into the DB
    std::fs::create_dir_all(db_path.join("library/added")).unwrap();
    std::fs::write(
        db_path.join("library/added/sequences.fna"),
        ">kraken:taxid|10116|rn6_mock\n\
         ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n",
    )
    .unwrap();
    // Required by kraken2-build --build: scan_fasta_file.pl output format is
    // "TAXID\t<seqid>\t<taxid>". build_kraken2_db.sh first globs
    // "prelim_map_*.txt" in library/added/ and cats them into prelim_map.txt,
    // so the file must match that pattern (not just "prelim_map.txt").
    std::fs::write(
        db_path.join("library/added/prelim_map_0001.txt"),
        "TAXID\tkraken:taxid|10116|rn6_mock\t10116\n",
    )
    .unwrap();

    let status = std::process::Command::new("kraken2-build")
        .args([
            "--build",
            "--db",
            db_path.to_str().unwrap(),
            "--no-masking",
            "--threads",
            "1",
        ])
        .status()
        .expect("failed to spawn kraken2-build");
    assert!(status.success(), "kraken2-build --build exited non-zero");

    db
}

/// End-to-end pipeline: prep -> kraken2 -> annotate --kraken-db -> filter.
///
/// Skipped when `kraken2` or `kraken2-build` are not on PATH (e.g. in a plain
/// `cargo test` without pixi). In CI the Testing job runs via `pixi run` so
/// both binaries are available.
///
/// The test uses 16 bp reads. Joined (R1 N R2) FASTA records are 33 bp, which
/// is below Kraken2's default k = 35, guaranteeing every read is unclassified.
#[test]
fn test_e2e_kraken2_pipeline() {
    if which::which("kraken2").is_err() {
        eprintln!("SKIP: kraken2 not on PATH");
        return;
    }
    if which::which("kraken2-build").is_err() {
        eprintln!("SKIP: kraken2-build not on PATH");
        return;
    }

    let _db = build_minidb();
    let db_path = _db.path();

    let dir = tempfile::tempdir().unwrap();

    let qname_sam = dir.path().join("qname.sam");
    {
        use std::io::Write as _;
        let mut f = std::fs::File::create(&qname_sam).unwrap();
        writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
        // flag 77  = paired | unmapped | mate-unmapped | read1
        // flag 141 = paired | unmapped | mate-unmapped | read2
        writeln!(
            f,
            "template1\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
        writeln!(
            f,
            "template1\t141\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTTTTTTTT\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
        writeln!(
            f,
            "template2\t77\t*\t0\t0\t*\t*\t0\t0\tCCCCCCCCCCCCCCCC\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
        writeln!(
            f,
            "template2\t141\t*\t0\t0\t*\t*\t0\t0\tGGGGGGGGGGGGGGGG\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
    }

    let fasta_out = dir.path().join("prep.fa");
    krak()
        .args([
            "prep",
            qname_sam.to_str().unwrap(),
            "-o",
            fasta_out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let fasta_content = std::fs::read_to_string(&fasta_out).unwrap();
    let fasta_records = fasta_content.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(
        fasta_records, 2,
        "expected 2 FASTA records (one per template): {fasta_content}"
    );

    let kraken_output = dir.path().join("kraken2.out");
    let k2_status = std::process::Command::new("kraken2")
        .args([
            "--db",
            db_path.to_str().unwrap(),
            "--output",
            kraken_output.to_str().unwrap(),
            fasta_out.to_str().unwrap(),
        ])
        .status()
        .expect("failed to spawn kraken2");
    assert!(k2_status.success(), "kraken2 exited non-zero: {k2_status}");

    let kraken_lines = std::fs::read_to_string(&kraken_output).unwrap();
    assert_eq!(
        kraken_lines.lines().count(),
        2,
        "expected 2 kraken2 output lines: {kraken_lines}"
    );
    let unclassified = kraken_lines.lines().filter(|l| l.starts_with('U')).count();
    assert_eq!(
        unclassified, 2,
        "expected both reads unclassified (16 bp < k=35): {kraken_lines}"
    );

    let annotated_sam = dir.path().join("annotated.sam");
    krak()
        .args([
            "annotate",
            "-i",
            qname_sam.to_str().unwrap(),
            "-a",
            kraken_output.to_str().unwrap(),
            "--kraken-db",
            db_path.to_str().unwrap(),
            "-o",
            annotated_sam.to_str().unwrap(),
        ])
        .assert()
        .success();

    let annotated_content = std::fs::read_to_string(&annotated_sam).unwrap();
    // Every alignment record must carry ti:i:0 (unclassified).
    let ti_zero = annotated_content
        .lines()
        .filter(|l| !l.starts_with('@') && l.contains("ti:i:0"))
        .count();
    assert_eq!(
        ti_zero, 4,
        "expected all 4 records to have ti:i:0: {annotated_content}"
    );
    // The full taxonomy tree must be embedded as a @CO comment.
    assert!(
        annotated_content
            .lines()
            .any(|l| l.starts_with("@CO\tkrak:report:")),
        "@CO krak:report: not found in annotated header: {annotated_content}"
    );

    let filtered_u = dir.path().join("filtered_unclassified.sam");
    let metrics1 = dir.path().join("metrics1.tsv");
    krak()
        .args([
            "filter",
            "-i",
            annotated_sam.to_str().unwrap(),
            "-o",
            filtered_u.to_str().unwrap(),
            "-m",
            metrics1.to_str().unwrap(),
            "-t",
            "10116",
            "--include-unclassified",
        ])
        .assert()
        .success();

    let out_u = std::fs::read_to_string(&filtered_u).unwrap();
    let kept_u = out_u.lines().filter(|l| !l.starts_with('@')).count();
    assert_eq!(
        kept_u, 4,
        "expected all 4 records kept with --include-unclassified: {out_u}"
    );

    let filtered_none = dir.path().join("filtered_none.sam");
    let metrics2 = dir.path().join("metrics2.tsv");
    krak()
        .args([
            "filter",
            "-i",
            annotated_sam.to_str().unwrap(),
            "-o",
            filtered_none.to_str().unwrap(),
            "-m",
            metrics2.to_str().unwrap(),
            "-t",
            "10116",
        ])
        .assert()
        .success();

    let out_none = std::fs::read_to_string(&filtered_none).unwrap();
    let kept_none = out_none.lines().filter(|l| !l.starts_with('@')).count();
    assert_eq!(
        kept_none, 0,
        "expected 0 records when reads are unclassified and -u not set: {out_none}"
    );
}

#[test]
fn test_n2ref_help() {
    krak().args(["n2ref", "--help"]).assert().success();
}

#[test]
fn test_n2ref_missing_reference_fails() {
    let input = NamedTempFile::new().unwrap();
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "n2ref",
            "--input",
            input.path().to_str().unwrap(),
            "--output",
            output.path().to_str().unwrap(),
            "--reference",
            "/nonexistent/ref.fa",
        ])
        .assert()
        .failure();
}

#[test]
fn n2ref_sniffs_sam_from_stdin() {
    let dir = tempfile::TempDir::new().unwrap();
    let fa = dir.path().join("ref.fa");
    let fai = dir.path().join("ref.fa.fai");
    std::fs::write(&fa, b">chr1\nACGT\n").unwrap();
    std::fs::write(&fai, b"chr1\t4\t6\t4\t5\n").unwrap();

    let sam_in = "@HD\tVN:1.6\tSO:unsorted\n\
                  @SQ\tSN:chr1\tLN:4\n\
                  r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNCGT\tIIII\n";
    let out = dir.path().join("out.sam");

    krak()
        .arg("n2ref")
        .arg("-")
        .arg(&out)
        .arg("-r")
        .arg(&fa)
        .write_stdin(sam_in)
        .assert()
        .success();

    let got = std::fs::read_to_string(&out).unwrap();
    let seq = got
        .lines()
        .find(|l| !l.starts_with('@'))
        .and_then(|l| l.split('\t').nth(9))
        .unwrap();
    assert_eq!(seq, "ACGT");
}

#[test]
fn test_report2tsv_positional_input_and_output() {
    let report = "100.00\t2000\t0\tR\t1\troot\n\
                  100.00\t2000\t0\tD\t2\t  Bacteria\n";
    let input = write_tmp(report);
    let output = NamedTempFile::new().unwrap();

    krak()
        .arg("report2tsv")
        .arg(input.path())
        .arg(output.path())
        .assert()
        .success();

    let got = std::fs::read_to_string(output.path()).unwrap();
    let mut lines = got.lines();
    assert_eq!(
        lines.next().unwrap(),
        "tax_id\tname\trank\tlevel\tparent_tax_id\tparent_rank\t\
         clade_count\tdirect_count\tdescendant_count\t\
         frac_clade\tfrac_direct\tfrac_descendant\t\
         minimizer_count\tdistinct_minimizer_count"
    );
    // root: clade=2000, direct=0, descendant=2000, fractions all over 2000.
    assert_eq!(
        lines.next().unwrap(),
        "1\troot\tR\t0\t\t\t2000\t0\t2000\t1\t0\t1\t\t"
    );
    // Bacteria: depth 1, parent root.
    assert_eq!(
        lines.next().unwrap(),
        "2\tBacteria\tD\t1\t1\tR\t2000\t0\t2000\t1\t0\t1\t\t"
    );
    assert!(lines.next().is_none());
}

#[test]
fn test_report2tsv_positional_and_flag_input_conflict_fails() {
    let input = write_tmp("100.00\t1\t0\tR\t1\troot\n");
    let output = NamedTempFile::new().unwrap();

    krak()
        .arg("report2tsv")
        .arg(input.path())
        .arg("-i")
        .arg(input.path())
        .arg(output.path())
        .assert()
        .failure();
}

fn write_tmp_fastq_gz(records: &[(&str, &str, &str)]) -> NamedTempFile {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    let f = tempfile::Builder::new()
        .suffix(".fq.gz")
        .tempfile()
        .unwrap();
    {
        let file = f.reopen().unwrap();
        let mut enc = GzEncoder::new(file, Compression::default());
        for (name, seq, qual) in records {
            writeln!(enc, "@{name}").unwrap();
            writeln!(enc, "{seq}").unwrap();
            writeln!(enc, "+").unwrap();
            writeln!(enc, "{qual}").unwrap();
        }
        enc.finish().unwrap();
    }
    f
}

#[test]
fn prep_reads_gzipped_fastq_file() {
    let in_fq_gz = write_tmp_fastq_gz(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let out = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();

    krak()
        .arg("prep")
        .arg("-i")
        .arg(in_fq_gz.path())
        .arg("-o")
        .arg(out.path())
        .assert()
        .success();

    let got = std::fs::read_to_string(out.path()).unwrap();
    assert_eq!(got, ">r1\nACGT\n>r2\nTTTT\n");
}

#[test]
fn filter_roundtrips_gzipped_fastq_files() {
    use flate2::bufread::MultiGzDecoder;
    use std::io::Read as _;

    let in_fq_gz = write_tmp_fastq_gz(&[("r1", "ACGT", "IIII"), ("r2", "TTTT", "JJJJ")]);
    let kraken = write_tmp_classifications(&[("r1", 9606), ("r2", 1234)]);
    let out = tempfile::Builder::new()
        .suffix(".fq.gz")
        .tempfile()
        .unwrap();

    krak()
        .arg("filter")
        .arg("-i")
        .arg(in_fq_gz.path())
        .arg("-o")
        .arg(out.path())
        .arg("-c")
        .arg(kraken.path())
        .arg("-t")
        .arg("9606")
        .arg("--per-record")
        .assert()
        .success();

    let f = std::fs::File::open(out.path()).unwrap();
    let mut dec = MultiGzDecoder::new(std::io::BufReader::new(f));
    let mut got = String::new();
    dec.read_to_string(&mut got).unwrap();
    assert_eq!(got, "@r1\nACGT\n+\nIIII\n");
}

#[test]
fn filter_reads_gzipped_stdin_writes_plain_stdout() {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let mut gz_bytes = Vec::new();
    {
        let mut enc = GzEncoder::new(&mut gz_bytes, Compression::default());
        enc.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
            .unwrap();
        enc.finish().unwrap();
    }

    let kraken = write_tmp_classifications(&[("r1", 9606), ("r2", 1234)]);

    let assert = krak()
        .arg("filter")
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg("-")
        .arg("-c")
        .arg(kraken.path().to_str().unwrap())
        .arg("-t")
        .arg("9606")
        .arg("--per-record")
        .write_stdin(gz_bytes)
        .assert()
        .success();

    let stdout = String::from_utf8(assert.get_output().stdout.clone()).unwrap();
    assert_eq!(stdout, "@r1\nACGT\n+\nIIII\n");
}

#[test]
fn prep_reads_gzipped_stdin() {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let mut gz_bytes = Vec::new();
    {
        let mut enc = GzEncoder::new(&mut gz_bytes, Compression::default());
        enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        enc.finish().unwrap();
    }

    let out = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();

    krak()
        .arg("prep")
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg(out.path())
        .write_stdin(gz_bytes)
        .assert()
        .success();

    let got = std::fs::read_to_string(out.path()).unwrap();
    assert_eq!(got, ">r1\nACGT\n");
}

#[test]
fn prep_reads_plain_fastq_stdin() {
    let out = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();

    krak()
        .arg("prep")
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg(out.path())
        .write_stdin(&b"@r1\nACGT\n+\nIIII\n"[..])
        .assert()
        .success();

    let got = std::fs::read_to_string(out.path()).unwrap();
    assert_eq!(got, ">r1\nACGT\n");
}

fn write_coord_sorted_bam(path: &Path, ti_tag: Option<u32>) {
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record_buf::{
        data::field::Value, QualityScores, RecordBuf, Sequence,
    };

    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:8\n"
        .parse()
        .unwrap();
    let mut w = bam::io::writer::Builder.build_from_path(path).unwrap();
    w.write_header(&header).unwrap();
    let mut r = RecordBuf::default();
    *r.name_mut() = Some("r1".as_bytes().into());
    *r.flags_mut() = noodles::sam::alignment::record::Flags::default();
    *r.reference_sequence_id_mut() = Some(0);
    *r.alignment_start_mut() = Some(noodles::core::Position::MIN);
    *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
    *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
    if let Some(ti) = ti_tag {
        r.data_mut().insert([b't', b'i'].into(), Value::UInt32(ti));
    }
    w.write_alignment_record(&header, &r).unwrap();
}

fn write_queryname_sorted_bam(path: &Path) {
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

    let header: sam::Header = "@HD\tVN:1.6\tSO:queryname\n".parse().unwrap();
    let mut w = bam::io::writer::Builder.build_from_path(path).unwrap();
    w.write_header(&header).unwrap();
    let mut r = RecordBuf::default();
    *r.name_mut() = Some("r1".as_bytes().into());
    *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
    *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
    w.write_alignment_record(&header, &r).unwrap();
}

#[test]
fn annotate_coord_sorted_bam_writes_sibling_bai() {
    let dir = tempfile::TempDir::new().unwrap();
    let in_bam = dir.path().join("in.bam");
    write_coord_sorted_bam(&in_bam, None);
    let assignments = write_tmp_classifications(&[("r1", 9606)]);
    let out_bam = dir.path().join("out.bam");

    krak()
        .args([
            "annotate",
            "-i",
            in_bam.to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            out_bam.to_str().unwrap(),
        ])
        .assert()
        .success();

    let bai = dir.path().join("out.bam.bai");
    assert!(bai.exists(), "expected sibling .bai at {}", bai.display());
    noodles::bam::bai::fs::read(&bai).expect(".bai must be readable");
}

#[test]
fn annotate_queryname_sorted_bam_does_not_write_bai() {
    let dir = tempfile::TempDir::new().unwrap();
    let in_bam = dir.path().join("in.bam");
    write_queryname_sorted_bam(&in_bam);
    let assignments = write_tmp_classifications(&[("r1", 9606)]);
    let out_bam = dir.path().join("out.bam");

    krak()
        .args([
            "annotate",
            "-i",
            in_bam.to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            out_bam.to_str().unwrap(),
        ])
        .assert()
        .success();

    let bai = dir.path().join("out.bam.bai");
    assert!(
        !bai.exists(),
        "queryname-sorted output must NOT produce a .bai sidecar; found {}",
        bai.display()
    );
}

#[test]
fn filter_coord_sorted_bam_writes_bai_for_output_and_rejects() {
    let dir = tempfile::TempDir::new().unwrap();
    let in_bam = dir.path().join("in.bam");
    // Write two coord-sorted records, one with ti=9606 (kept), one with ti=1234 (rejected).
    {
        use noodles::bam;
        use noodles::sam;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::{
            data::field::Value, QualityScores, RecordBuf, Sequence,
        };
        let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:64\n"
            .parse()
            .unwrap();
        let mut w = bam::io::writer::Builder.build_from_path(&in_bam).unwrap();
        w.write_header(&header).unwrap();
        for (name, pos, ti) in [("r1", 1u32, 9606u32), ("r2", 5u32, 1234u32)] {
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(name.as_bytes().into());
            *r.flags_mut() = noodles::sam::alignment::record::Flags::default();
            *r.reference_sequence_id_mut() = Some(0);
            *r.alignment_start_mut() =
                Some(noodles::core::Position::try_from(pos as usize).unwrap());
            *r.sequence_mut() = Sequence::from(b"ACGT".to_vec());
            *r.quality_scores_mut() = QualityScores::from(vec![30u8; 4]);
            r.data_mut().insert([b't', b'i'].into(), Value::UInt32(ti));
            w.write_alignment_record(&header, &r).unwrap();
        }
    }
    let out_bam = dir.path().join("out.bam");
    let rej_bam = dir.path().join("rej.bam");

    krak()
        .args([
            "filter",
            "-i",
            in_bam.to_str().unwrap(),
            "-o",
            out_bam.to_str().unwrap(),
            "-r",
            rej_bam.to_str().unwrap(),
            "-t",
            "9606",
            "--per-record",
        ])
        .assert()
        .success();

    let out_bai = dir.path().join("out.bam.bai");
    let rej_bai = dir.path().join("rej.bam.bai");
    assert!(
        out_bai.exists(),
        "expected sibling .bai for output: {}",
        out_bai.display()
    );
    assert!(
        rej_bai.exists(),
        "expected sibling .bai for rejects: {}",
        rej_bai.display()
    );
    noodles::bam::bai::fs::read(&out_bai).expect("output .bai must be readable");
    noodles::bam::bai::fs::read(&rej_bai).expect("rejects .bai must be readable");
}

#[test]
fn annotate_to_stdout_does_not_write_bai() {
    // -o - resolves to /dev/stdout; the indexer must skip silently.
    let dir = tempfile::TempDir::new().unwrap();
    let in_bam = dir.path().join("in.bam");
    write_coord_sorted_bam(&in_bam, None);
    let assignments = write_tmp_classifications(&[("r1", 9606)]);

    krak()
        .args([
            "annotate",
            "-i",
            in_bam.to_str().unwrap(),
            "-a",
            assignments.path().to_str().unwrap(),
            "-o",
            "-",
        ])
        .assert()
        .success();

    // No sidecar should be created next to the binary either.
    let cwd_bai = std::path::Path::new("/dev/stdout.bai");
    assert!(!cwd_bai.exists(), "/dev/stdout.bai must not exist");
}

#[test]
fn prep_reads_bam_stdin() {
    use noodles::bam;
    use noodles::sam;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

    // Build a real BAM in a temp file, then read its bytes for stdin.
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
    }
    let bam_bytes = std::fs::read(&bam_path).unwrap();

    let out = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();

    krak()
        .arg("prep")
        .arg("-i")
        .arg("-")
        .arg("-o")
        .arg(out.path())
        .write_stdin(bam_bytes)
        .assert()
        .success();

    let got = std::fs::read_to_string(out.path()).unwrap();
    assert_eq!(got, ">r1\nACGT\n");
}

#[test]
fn test_prep_per_record_on_paired_sam_emits_each_record() {
    // Header has no SO:queryname. Records have the paired flag (0x1) set.
    // Without --per-record this would bail; with --per-record each record is
    // emitted as a single-end template.
    let header = ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:8"];
    let records = [
        "pair1\t77\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII",
        "pair1\t141\tchr1\t1\t60\t4M\t*\t0\t0\tTTTT\tIIII",
    ];
    let input = write_tmp_sam(&header, &records);
    let output = NamedTempFile::new().unwrap();
    krak()
        .args([
            "prep",
            input.path().to_str().unwrap(),
            "--per-record",
            "-o",
            output.path().to_str().unwrap(),
        ])
        .assert()
        .success();
    let out = std::fs::read_to_string(output.path()).unwrap();
    assert!(out.contains(">pair1/1\nAAAA\n"), "got: {out}");
    assert!(out.contains(">pair1/2\nTTTT\n"), "got: {out}");
    // Two distinct templates, not one paired template
    assert_eq!(out.matches(">pair1/").count(), 2);
}
