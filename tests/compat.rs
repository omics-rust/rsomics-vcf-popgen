use std::path::PathBuf;
use std::process::Command;

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bin() -> String {
    env!("CARGO_BIN_EXE_rsomics-vcf-popgen").to_string()
}

fn run(args: &[&str]) -> String {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("run rsomics-vcf-popgen");
    assert!(
        out.status.success(),
        "non-zero exit:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

// ─── freq ────────────────────────────────────────────────────────────────────

#[test]
fn freq_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["freq", vcf.to_str().unwrap()]);
    assert!(out.starts_with("CHROM\tPOS\tN_ALLELES\tN_CHR"));
}

#[test]
fn freq_counts_correct() {
    let vcf = golden("small.vcf");
    let out = run(&["freq", vcf.to_str().unwrap()]);
    // chr1:100 has 5 samples, all called. alt (G) count = 4 (0/1 + 1/1 + 0/1 = 1+2+1 = 4)
    // n_chr = 10, n_alt = 4, freq_alt = 0.4000
    let lines: Vec<&str> = out.lines().collect();
    let site100 = lines
        .iter()
        .find(|l| l.contains("chr1\t100\t"))
        .expect("site chr1:100");
    assert!(site100.contains("G:0.4000"), "alt freq at 100: {site100}");
}

#[test]
fn freq_skips_monomorphic() {
    let vcf = golden("small.vcf");
    let out = run(&["freq", vcf.to_str().unwrap()]);
    // chr1:400 is all ref (0/0 for all 5 samples); alt freq = 0.0, still emitted
    // (vcftools --freq does emit monomorphic sites)
    let lines: Vec<&str> = out.lines().collect();
    let site400 = lines.iter().find(|l| l.contains("chr1\t400\t"));
    assert!(
        site400.is_some(),
        "monomorphic site 400 should still be emitted"
    );
}

// ─── het ─────────────────────────────────────────────────────────────────────

#[test]
fn het_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["het", vcf.to_str().unwrap()]);
    assert!(out.starts_with("INDV\tO(HOM)\tE(HOM)\tN_SITES\tF"));
}

#[test]
fn het_f_in_range() {
    let vcf = golden("small.vcf");
    let out = run(&["het", vcf.to_str().unwrap()]);
    let lines: Vec<&str> = out.lines().skip(1).collect();
    // 5 samples expected.
    assert_eq!(lines.len(), 5, "5 samples in het output");
    for line in &lines {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 5, "5 columns in het output: {line}");
        let f: f64 = cols[4].parse().unwrap_or(f64::NAN);
        if !f.is_nan() {
            assert!((-2.0..=2.0).contains(&f), "F out of plausible range: {f}");
        }
    }
}

// ─── hardy ───────────────────────────────────────────────────────────────────

#[test]
fn hardy_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["hardy", vcf.to_str().unwrap()]);
    assert!(out.starts_with("CHROM\tPOS\tOBS(HOM1/HET/HOM2)"));
}

#[test]
fn hardy_p_in_range() {
    let vcf = golden("small.vcf");
    let out = run(&["hardy", vcf.to_str().unwrap()]);
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(cols.len() >= 8, "8 cols in hardy output: {line}");
        let p_hwe: f64 = cols[5].parse().unwrap_or(f64::NAN);
        if !p_hwe.is_nan() {
            assert!((0.0..=1.0).contains(&p_hwe), "P_HWE out of [0,1]: {p_hwe}");
        }
    }
}

// ─── missing-site ────────────────────────────────────────────────────────────

#[test]
fn missing_site_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["missing-site", vcf.to_str().unwrap()]);
    assert!(out.starts_with("CHROM\tPOS\tN_DATA"));
}

#[test]
fn missing_site_counts_missing() {
    let vcf = golden("small.vcf");
    let out = run(&["missing-site", vcf.to_str().unwrap()]);
    // chr1:300 has 1 missing (./.); chr1:500 has 1 missing (./.); others 0.
    let lines: Vec<&str> = out.lines().skip(1).collect();
    // 8 sites total in small.vcf.
    assert_eq!(lines.len(), 8);
    let site300 = lines.iter().find(|l| l.contains("chr1\t300\t")).unwrap();
    let cols: Vec<&str> = site300.split('\t').collect();
    assert_eq!(cols[4], "1", "1 missing at chr1:300");
    assert_eq!(cols[5], "0.2000", "20% missing at chr1:300");
}

// ─── missing-indv ────────────────────────────────────────────────────────────

#[test]
fn missing_indv_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["missing-indv", vcf.to_str().unwrap()]);
    assert!(out.starts_with("INDV\tN_DATA"));
}

#[test]
fn missing_indv_counts() {
    let vcf = golden("small.vcf");
    let out = run(&["missing-indv", vcf.to_str().unwrap()]);
    let lines: Vec<&str> = out.lines().skip(1).collect();
    assert_eq!(lines.len(), 5, "5 samples");
    // SAMPLE2 has ./. at pos 500 → 1 missing.
    let s2 = lines.iter().find(|l| l.starts_with("SAMPLE2")).unwrap();
    let cols: Vec<&str> = s2.split('\t').collect();
    assert_eq!(cols[3], "1", "SAMPLE2 has 1 missing");
}

// ─── pi ──────────────────────────────────────────────────────────────────────

#[test]
fn pi_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["pi", vcf.to_str().unwrap(), "--window", "1000000"]);
    assert!(out.starts_with("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI"));
}

#[test]
fn pi_nonnegative() {
    let vcf = golden("small.vcf");
    let out = run(&["pi", vcf.to_str().unwrap(), "--window", "1000000"]);
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let pi: f64 = cols[4].parse().unwrap_or(f64::NAN);
        assert!(pi >= 0.0, "pi should be non-negative: {pi}");
    }
}

// ─── singleton ───────────────────────────────────────────────────────────────

#[test]
fn singleton_has_header() {
    let vcf = golden("small.vcf");
    let out = run(&["singleton", vcf.to_str().unwrap()]);
    assert!(out.starts_with("CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV"));
}

#[test]
fn singleton_type_is_s_or_d() {
    let vcf = golden("small.vcf");
    let out = run(&["singleton", vcf.to_str().unwrap()]);
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols[2] == "S" || cols[2] == "D",
            "type must be S or D: {}",
            cols[2]
        );
    }
}

// ─── fst (differential vs vcftools) ──────────────────────────────────────────

// Weir & Cockerham FST has no second self-checkable reference, so this runs the
// real vcftools when present and diffs per-site values. Pinned to 0.1.17 — other
// versions can drift, so loud-skip rather than fail on a mismatched binary.
#[test]
fn fst_matches_vcftools_0_1_17() {
    let Some(vcftools) = vcftools_0_1_17() else {
        eprintln!("SKIP fst differential: vcftools 0.1.17 not on PATH");
        return;
    };

    let dir = tempfile::tempdir_in(scratch()).unwrap();
    let d = dir.path();
    write_fst_fixture(d);

    let vt = Command::new(&vcftools)
        .args(["--vcf"])
        .arg(d.join("in.vcf"))
        .arg("--weir-fst-pop")
        .arg(d.join("pop1.txt"))
        .arg("--weir-fst-pop")
        .arg(d.join("pop2.txt"))
        .arg("--out")
        .arg(d.join("vt"))
        .output()
        .expect("run vcftools");
    assert!(vt.status.success(), "vcftools failed");

    let ours = run(&[
        "fst",
        d.join("in.vcf").to_str().unwrap(),
        "--pop",
        d.join("pop1.txt").to_str().unwrap(),
        "--pop",
        d.join("pop2.txt").to_str().unwrap(),
    ]);

    let want = parse_fst(&std::fs::read_to_string(d.join("vt.weir.fst")).unwrap());
    let got = parse_fst(&ours);
    assert_eq!(want.len(), got.len(), "site count differs");
    for ((wk, wv), (gk, gv)) in want.iter().zip(&got) {
        assert_eq!(wk, gk, "site key differs");
        match (wv, gv) {
            (Some(a), Some(b)) => assert!(
                (a - b).abs() < 1e-6 || (a.abs() < 1e-6 && b.abs() < 1e-6),
                "site {wk:?}: vcftools {a} vs ours {b}"
            ),
            (None, Some(b)) | (Some(b), None) => {
                assert!(b.abs() < 1e-6, "site {wk:?}: one nan, other {b}")
            }
            (None, None) => {}
        }
    }
}

fn vcftools_0_1_17() -> Option<String> {
    let out = Command::new("vcftools").arg("--version").output().ok()?;
    let text = String::from_utf8_lossy(&out.stderr) + String::from_utf8_lossy(&out.stdout);
    text.contains("0.1.17").then(|| "vcftools".to_string())
}

fn scratch() -> String {
    std::env::var("TMPDIR").unwrap_or_else(|_| "/Volumes/KIOXIA/tmp".to_string())
}

// Per-site (chrom, pos) -> Fst, with -nan/nan mapped to None.
fn parse_fst(s: &str) -> Vec<((String, String), Option<f64>)> {
    s.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let v = match c[2] {
                "-nan" | "nan" => None,
                x => x.parse::<f64>().ok(),
            };
            ((c[0].to_string(), c[1].to_string()), v)
        })
        .collect()
}

fn write_fst_fixture(d: &std::path::Path) {
    // Two populations of four, allele frequencies drawn per population so sites
    // span the full differentiation range (Fst near 0, negative, and near 1).
    let pop1 = ["S1", "S2", "S3", "S4"];
    let pop2 = ["S5", "S6", "S7", "S8"];
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\n",
    );
    let freqs1 = [0.1, 0.5, 0.9, 0.3, 0.7, 0.5, 0.2, 0.8];
    let freqs2 = [0.9, 0.5, 0.1, 0.7, 0.3, 0.5, 0.8, 0.2];
    let mut seed = 0x1234u64;
    let mut rng = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((seed >> 33) as f64) / (1u64 << 31) as f64
    };
    for s in 0..40usize {
        let f1 = freqs1[s % freqs1.len()];
        let f2 = freqs2[s % freqs2.len()];
        let mut gts = Vec::new();
        for _ in 0..4 {
            let a1 = u8::from(rng() < f1);
            let a2 = u8::from(rng() < f1);
            gts.push(format!("{a1}/{a2}"));
        }
        for _ in 0..4 {
            let a1 = u8::from(rng() < f2);
            let a2 = u8::from(rng() < f2);
            gts.push(format!("{a1}/{a2}"));
        }
        vcf.push_str(&format!(
            "1\t{}\t.\tA\tG\t.\t.\t.\tGT\t{}\n",
            100 + s * 100,
            gts.join("\t")
        ));
    }
    std::fs::write(d.join("in.vcf"), vcf).unwrap();
    std::fs::write(d.join("pop1.txt"), pop1.join("\n")).unwrap();
    std::fs::write(d.join("pop2.txt"), pop2.join("\n")).unwrap();
}

// ─── CLI ─────────────────────────────────────────────────────────────────────

#[test]
fn cli_help_exits_zero() {
    let out = Command::new(bin())
        .arg("--help")
        .output()
        .expect("run rsomics-vcf-popgen --help");
    assert!(out.status.success(), "help should exit 0");
}
