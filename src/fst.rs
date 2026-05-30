//! Weir & Cockerham (1984) FST per site — equivalent to `vcftools --weir-fst-pop`.
//!
//! Per-site theta = a/(a+b+c) using the variance components of W&C 1984 eqs 2-4.
//! Every term is symmetric under p -> 1-p, so a biallelic site needs no convention
//! for which allele is focal. The mean estimate averages per-site ratios; the
//! weighted estimate is sum(a)/sum(a+b+c) across sites — both reported like vcftools.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::vcf::{VcfReader, VcfRecord};

pub struct FstRecord {
    pub chrom: String,
    pub pos: u64,
    pub fst: f64, // NaN where undefined (a+b+c == 0, or fewer than two populations present)
}

pub struct FstSummary {
    pub records: Vec<FstRecord>,
    pub mean: f64,     // average of per-site ratios over defined sites
    pub weighted: f64, // sum(a) / sum(a+b+c)
}

fn read_pop(path: &Path) -> Result<Vec<String>> {
    let f = File::open(path)
        .with_context(|| format!("cannot open population file {}", path.display()))?;
    let mut ids = Vec::new();
    for line in BufReader::new(f).lines() {
        let id = line?.trim().to_string();
        if !id.is_empty() {
            ids.push(id);
        }
    }
    Ok(ids)
}

pub fn fst(input: &Path, pop_files: &[PathBuf]) -> Result<FstSummary> {
    let reader = VcfReader::open(input)?;
    let samples = reader.samples.clone();
    let pops: Vec<Vec<usize>> = pop_files
        .iter()
        .map(|pf| -> Result<Vec<usize>> {
            let ids: HashSet<String> = read_pop(pf)?.into_iter().collect();
            Ok(samples
                .iter()
                .enumerate()
                .filter(|(_, s)| ids.contains(*s))
                .map(|(i, _)| i)
                .collect())
        })
        .collect::<Result<_>>()?;

    let mut records = Vec::new();
    let mut sum_a = 0.0f64;
    let mut sum_abc = 0.0f64;
    let mut mean_acc = 0.0f64;
    let mut mean_n = 0u64;

    for rec in reader {
        let rec = rec?;
        if rec.alt_alleles.len() != 1 {
            continue; // vcftools computes Fst on biallelic sites only
        }
        let fst = match site_components(&rec, &pops) {
            Some((a, b, c)) => {
                let denom = a + b + c;
                sum_a += a;
                sum_abc += denom;
                if denom != 0.0 {
                    mean_acc += a / denom;
                    mean_n += 1;
                    a / denom
                } else {
                    f64::NAN
                }
            }
            None => f64::NAN,
        };
        records.push(FstRecord {
            chrom: rec.chrom,
            pos: rec.pos,
            fst,
        });
    }

    Ok(FstSummary {
        records,
        mean: if mean_n > 0 {
            mean_acc / mean_n as f64
        } else {
            f64::NAN
        },
        weighted: if sum_abc != 0.0 {
            sum_a / sum_abc
        } else {
            f64::NAN
        },
    })
}

/// W&C variance components (a, b, c) for one biallelic site. `None` when fewer than
/// two populations carry a non-missing genotype, or the sample-size terms degenerate.
fn site_components(rec: &VcfRecord, pops: &[Vec<usize>]) -> Option<(f64, f64, f64)> {
    let mut n = Vec::new(); // diploid sample size per population
    let mut p = Vec::new(); // alt-allele frequency per population
    let mut h = Vec::new(); // observed heterozygote proportion per population
    for pop in pops {
        let (mut ni, mut alt, mut het) = (0u64, 0u64, 0u64);
        for &si in pop {
            let gt = rec.gts[si];
            if gt.is_missing() {
                continue;
            }
            ni += 1;
            alt += u64::from(!gt.a1.is_ref()) + u64::from(!gt.a2.is_ref());
            het += u64::from(gt.is_het());
        }
        if ni > 0 {
            n.push(ni as f64);
            p.push(alt as f64 / (2.0 * ni as f64));
            h.push(het as f64 / ni as f64);
        }
    }

    let r = n.len() as f64;
    if n.len() < 2 {
        return None;
    }
    let n_tot: f64 = n.iter().sum();
    let n_bar = n_tot / r;
    if n_bar <= 1.0 {
        return None;
    }
    let sum_n2: f64 = n.iter().map(|x| x * x).sum();
    let n_c = (n_tot - sum_n2 / n_tot) / (r - 1.0);
    if n_c == 0.0 {
        return None;
    }
    let p_bar = n.iter().zip(&p).map(|(ni, pi)| ni * pi).sum::<f64>() / n_tot;
    let s2 = n
        .iter()
        .zip(&p)
        .map(|(ni, pi)| ni * (pi - p_bar).powi(2))
        .sum::<f64>()
        / ((r - 1.0) * n_bar);
    let h_bar = n.iter().zip(&h).map(|(ni, hi)| ni * hi).sum::<f64>() / n_tot;

    let pq = p_bar * (1.0 - p_bar);
    let a = (n_bar / n_c) * (s2 - (pq - ((r - 1.0) / r) * s2 - h_bar / 4.0) / (n_bar - 1.0));
    let b = (n_bar / (n_bar - 1.0))
        * (pq - ((r - 1.0) / r) * s2 - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar);
    let c = h_bar / 2.0;
    Some((a, b, c))
}

pub fn print_fst(summary: &FstSummary) {
    println!("CHROM\tPOS\tWEIR_AND_COCKERHAM_FST");
    for r in &summary.records {
        println!("{}\t{}\t{}", r.chrom, r.pos, fmt_g(r.fst));
    }
    eprintln!(
        "Weir and Cockerham mean Fst estimate: {}",
        fmt_g(summary.mean)
    );
    eprintln!(
        "Weir and Cockerham weighted Fst estimate: {}",
        fmt_g(summary.weighted)
    );
}

/// Six significant figures, %g-style — matches the C++ ostream default vcftools writes.
fn fmt_g(x: f64) -> String {
    if x.is_nan() {
        return "-nan".to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let exp = x.abs().log10().floor() as i32;
    if (-4..6).contains(&exp) {
        let decimals = (5 - exp).max(0) as usize;
        let mut s = format!("{x:.decimals$}");
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.pop();
            }
        }
        s
    } else {
        let mut s = format!("{x:.5e}");
        let e = s.find('e').unwrap();
        let real_exp: i32 = s[e + 1..].parse().unwrap();
        let mut mant = s[..e].to_string();
        if mant.contains('.') {
            while mant.ends_with('0') {
                mant.pop();
            }
            if mant.ends_with('.') {
                mant.pop();
            }
        }
        s = format!(
            "{mant}e{}{:02}",
            if real_exp < 0 { '-' } else { '+' },
            real_exp.abs()
        );
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Two populations, one perfectly-differentiated site: pop1 all 0/0, pop2 all 1/1.
    // Fixation gives Fst = 1 by every estimator.
    #[test]
    fn fixed_difference_is_one() {
        let dir = tempfile::tempdir_in(scratch()).unwrap();
        let vcf = dir.path().join("f.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n\
             1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t1/1\n",
        )
        .unwrap();
        let pop1 = dir.path().join("p1.txt");
        let pop2 = dir.path().join("p2.txt");
        std::fs::write(&pop1, "A\nB\n").unwrap();
        std::fs::write(&pop2, "C\nD\n").unwrap();

        let s = fst(&vcf, &[pop1, pop2]).unwrap();
        assert_eq!(s.records.len(), 1);
        assert!(
            (s.records[0].fst - 1.0).abs() < 1e-9,
            "got {}",
            s.records[0].fst
        );
        assert!((s.weighted - 1.0).abs() < 1e-9);
    }

    // No between-population structure (both pops identical 0/1) → numerator ~0.
    #[test]
    fn no_structure_is_near_zero_or_negative() {
        let dir = tempfile::tempdir_in(scratch()).unwrap();
        let vcf = dir.path().join("f.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n\
             1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\t0/1\n",
        )
        .unwrap();
        let pop1 = dir.path().join("p1.txt");
        let pop2 = dir.path().join("p2.txt");
        std::fs::write(&pop1, "A\nB\n").unwrap();
        std::fs::write(&pop2, "C\nD\n").unwrap();

        let s = fst(&vcf, &[pop1, pop2]).unwrap();
        assert!(s.records[0].fst <= 1e-9, "got {}", s.records[0].fst);
    }

    fn scratch() -> String {
        std::env::var("TMPDIR").unwrap_or_else(|_| "/Volumes/KIOXIA/tmp".to_string())
    }
}
