# dmrseqx

Detection of differentially methylated regions with the procedure of
[dmrseq](https://bioconductor.org/packages/dmrseq/) (Korthauer et al., Biostatistics
2019), written in C and run in threads. A whole-genome comparison takes minutes instead
of days, with the same regions and nearly the same calls. It compares two or more groups
or follows a continuous covariate, as dmrseq does.

dmrseqx is an independent reimplementation and is not affiliated with the authors of
dmrseq. Please also cite their paper when you use it.

## Installation

```r
# install.packages("BiocManager")
BiocManager::install(c("bsseq", "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics"))
install.packages("remotes")
remotes::install_github("JiantaoShi/dmrseqx")
```

A C compiler is needed. Threads use POSIX threads.

## Usage

The arguments are those of `dmrseq::dmrseq()`, with `threads` in place of `BPPARAM`:

```r
library(dmrseqx)
regions <- dmrseqx(bs, testCovariate = "group", cutoff = 0.1, threads = 8)
regions[which(regions$qval < 0.05)]
```

`bs` is a sorted `BSseq` object with raw counts; as for dmrseq, loci without coverage in
all samples of a group (for a continuous covariate, in all samples but one) must be removed
first. The result is a `GRanges` object with the
same columns as from dmrseq: `L`, `area`, `beta`, `stat`, `pval`, `qval`, `index`.

| Argument | Default | Meaning |
|---|---|---|
| `testCovariate` | | column of `pData(bs)`: two groups, several groups, or continuous |
| `adjustCovariate` | `NULL` | columns of covariates to adjust for |
| `matchCovariate` | `NULL` | column on which permutations are matched, for two groups |
| `cutoff` | `0.1` | cutoff on the smoothed difference for candidates |
| `minNumRegion` | `5` | smallest number of CpGs in a region |
| `smooth`, `bpSpan`, `minInSpan`, `maxGapSmooth` | `TRUE`, `1000`, `30`, `2500` | smoothing |
| `maxGap` | `1000` | largest distance between CpGs of a region |
| `maxPerms` | `10` | largest number of permutations |
| `stat` | `"stat"` | statistic for the p-values |
| `threads` | `1` | threads; results don't depend on them |
| `fit` | `"dmrseq"` | `"global"`: true optimum in short regions, see below |
| `perms` | `NULL` | matrix of permutations instead of drawing them |

Use `set.seed()` for reproducible results. `block = TRUE` is not supported; dmrseq handles it.

### Several groups and continuous covariates

```r
# a global test among three or more groups: a factor or character column
regions <- dmrseqx(bs, testCovariate = "tissue", cutoff = 0.2, threads = 8)

# a continuous covariate: a numeric column with more than two values
bsseq::pData(bs)$age10 <- (bsseq::pData(bs)$age - 50) / 10
regions <- dmrseqx(bs, testCovariate = "age10", cutoff = 0.05, threads = 8)
```

These follow the matching paths of dmrseq, whose code marks them as experimental. They behave
differently from two groups in ways that matter for the choice of `cutoff`:

- **Several groups.** The effect of a CpG is the largest minus the smallest fitted proportion
  among the groups. It is never negative, so regions have no direction, and `stat`, the square
  root of the sequential F statistic of the groups, is never negative either. Noise alone can
  reach the default cutoff when groups are small: with 2 samples in each of 4 groups, 86% of
  6,000 simulated CpGs fell into 10 regions. Raise `cutoff` in such a case. The result has one
  column `beta_<level>` per level except the first, in place of `beta`.
- **Continuous.** The effect of a CpG is the change of the fitted proportion between the values
  0 and 1 of the covariate. `cutoff` therefore depends on both the units and the origin of the
  covariate, and a covariate far from 0, such as age in years, gives effects at extrapolated
  values. Center and scale it so that this change means something. `stat` is the t statistic
  of the slope and `beta` the slope.
- **Both.** Permutations reorder the samples without restriction, and adjustment covariates are
  reordered with them for the region model. Extreme coefficients aren't refitted.
- **Missing coverage.** A permutation can leave a group without coverage at a CpG, even if the
  observed labels don't. Such a CpG has no effect of its own in that pass, and the smoother passes
  over it. dmrseq usually stops with an error here (`singular matrix in 'backsolve'`), and also
  when a chromosome has only one pattern of samples without coverage (`seq_len(nrow(NA.patterns))`).
  To run both on the same data, keep CpGs that lack coverage in fewer samples than the smallest
  group has.

## Comparison with dmrseq

ESCC tumors and normal tissue, 10 vs 9 samples, 25,321,525 CpGs covered in all samples,
default settings, and the permutation labels of the dmrseq run. As the reference we used
dmrseq with a sparse fit for regions of at least 40 CpGs, since the unmodified code was
projected to need about 4 days at 48 workers; on chromosomes 1 and 10 this reference agrees
with unmodified dmrseq, with identical statistics below 40 CpGs and differences up to 1.5e-4
above.

| | dmrseq with sparse fit | dmrseqx |
|---|---|---|
| Time, observed labels and 10 permutations | 32,112 s at 47 workers | 258 s at 1 thread, 69 s at 8 threads |
| Peak memory of the main process | 16.7 GB | 10.7 GB |
| Regions, observed labels | 288,254 | 288,254, all with the same coordinates |
| Regions, 10 permutations | 820,604 | 820,604, all with the same coordinates |
| Statistic, regions of at least 40 CpGs | | largest difference 2.5e-5 |
| Statistic, regions of less than 40 CpGs | | 97.9% to 98.2% agree to 1e-3 per labelling |
| Calls at q <= 0.05 | 161,242 | 161,253, of which 160,722 in both (0.65% differ) |

Reproduce with `inst/validation/check_genome.R`. The time of dmrseqx is the call of `dmrseqx()`,
without reading the data.

Chromosome 21 with 3 vs 3 samples, compared directly with unmodified dmrseq 1.32.0
(`inst/validation/check_original.R`), all 9 relabellings: all 3,434 regions of the observed
labels are the same; for at least 40 CpGs the largest difference of the statistic is 4.6e-5;
below that 98% agree to 1e-3; calls 1,525 and 1,522, of which 1,520 are in both. It took 2.4 s
at one thread against 493 s at 16 workers.

On smaller data (`inst/validation/check_dmrseq.R`), including missing coverage, adjustment
and matching covariates and unbalanced groups, compared directly with dmrseq 1.32.0 at one
thread each: all region coordinates are the same in all 11 settings, all calls are the same,
and dmrseqx is 50 to 300 times faster.

Several groups and a continuous covariate, compared pass by pass with `dmrseq:::bumphunt` at one
thread each (`inst/validation/check_general.R`): the observed labels and 2 permutations on simulated
data (up to 6,000 CpGs, 8 or 9 samples, cutoff 0.1 unless noted), the observed labels and 1
permutation on ESCC chromosome 21 (364,865 CpGs, 9 samples, cutoff 0.2). All regions have the same
coordinates in all settings, and the area differs by at most 1.3e-7.

| Setting | Regions | Statistic, at least 40 CpGs: largest difference | Statistic, less than 40 CpGs: agree to 1e-3 | dmrseq | dmrseqx |
|---|---|---|---|---|---|
| 3 groups | 191 | 1.6e-5 | 87 of 88 | 49 s | 0.14 s |
| 4 groups of 2 samples | 36 | 9.6e-6 | 8 of 8 | 16,905 s | 0.10 s |
| 3 groups, 10% missing (a) | 154 | 1.7e-5 | 61 of 62 | 131 s | 0.28 s |
| 3 groups, 30% missing (a) | 93 | 3.3e-5 | 24 of 24 | 551 s | 0.17 s |
| 3 groups, 10% missing, filtered (b) | 162 | 2.0e-5 | 73 of 73 | 55 s | 0.25 s |
| 3 groups, 30% missing, filtered (b) | 111 | 1.6e-4 | 62 of 62 | 40 s | 0.17 s |
| 3 groups, adjusted | 152 | 8.9e-6 | 56 of 58 | 377 s | 0.09 s |
| Continuous | 66 | 2.3e-4 | 33 of 34 | 7.0 s | 0.09 s |
| Continuous, values far from 0, cutoff 0.005 | 80 | 2.8e-4 | 53 of 53 | 10.6 s | 0.08 s |
| Continuous, 10% missing | 60 | 1.5e-5 | 37 of 37 | 16.8 s | 0.12 s |
| Continuous, 30% missing | 64 | 2.5e-4 | 45 of 48 | 25.6 s | 0.10 s |
| Continuous, adjusted | 64 | 2.1e-4 | 21 of 22 | 9.5 s | 0.07 s |
| ESCC chromosome 21, 3 groups | 8,931 | 3.6e-5 | 7,816 of 7,959 (98.2%) | 1,052 s | 3.1 s |
| ESCC chromosome 21, continuous | 4,405 | 3.4e-5 | 3,979 of 4,059 (98.0%) | 290 s | 1.2 s |

(a) dmrseq stops in both permutations, where a group has no coverage at 15 to 320 CpGs. The
reference for these passes is dmrseq with NA at these CpGs (`getEstimateNA` in the script).
(b) Only CpGs that lack coverage in fewer samples than the smallest group, so that dmrseq runs
unchanged.

Differences above 1e-4 are at statistics of 30 to 40, at most 7e-6 in relative terms. The time of
dmrseq with 4 groups comes from its dense fit of regions with about 2,000 CpGs, which noise produces
with so few samples at this cutoff. It isn't a typical gain in speed.

This shows that dmrseqx gives the results of dmrseq on these paths. We haven't evaluated the false
discovery rate on them, for either program.

### Regions of less than 40 CpGs

This is the only source of larger differences, and it comes from how dmrseq fits these
regions. Their errors have a correlation that decays with the distance in base pairs,
and `nlme::gls()` starts its search at a correlation parameter of 0.2. When the CpGs of
a region are about 8 bp or more apart, the correlation at this value is practically zero
between all CpGs, the likelihood is flat, and the optimizer stops at once. About a third
of these regions therefore keep the start value, which means no correlation, and have larger
statistics than with the optimum. Whether the optimizer stops also depends on rounding
errors for about 4% of them: with objective values that differ by only 1e-12, one region
went from t = -11.09 to -7.14. In dmrseq itself, 2 of 798 such regions changed when only
the order of the samples was changed.

`fit = "dmrseq"` reproduces this with a rule on the likelihood at the start value (stop if
g^2 / (2 |f|) < 3e-11 for the value f and the derivative g of the objective). This gives the
same statistic for about 98% of these regions, about what the optimizer of `gls()` itself gives
when run on an objective that is the same up to rounding.

`fit = "global"` uses the optimum of the restricted likelihood everywhere. This is what the
model asks for, but statistics and calls then differ from dmrseq for these regions.

## Method

The steps of a two-group comparison follow dmrseq:

1. **Effect of each CpG.** The difference of the pooled methylation of the groups, and a
   weight from coverage and a robust standard deviation.
2. **Smoothing.** Local quadratic regression with tricube kernel, evaluated at the vertices
   of an adaptive tree and interpolated between them with cubic Hermite polynomials (the
   evaluation structure of Loader 1999, used by locfit). This interpolation changes region
   boundaries, so that fitting at every CpG doesn't reproduce the regions of dmrseq.
3. **Candidates.** Runs of CpGs with a smoothed difference beyond `cutoff`, refined at the
   edges by the raw difference and by the shape of the effect.
4. **Test.** `asin(2 (M + 0.05) / (N + 0.1) - 1) ~ group + locus + covariates`, by restricted
   maximum likelihood, with variance proportional to the inverse median coverage of a locus
   and correlated errors within samples: AR(1) over rows for at least 40 CpGs, and a correlation
   decaying with the distance for less. Regions with extreme coefficients are refitted
   without variance weights, and coverage outliers in short regions remove them too.
5. **Significance.** The same steps on permutations of the labels give a pooled null,
   empirical p-values and Benjamini-Hochberg q-values.

For several groups and continuous covariates, step 1 is a linear regression of the logit of
the methylation proportion on the design, at each CpG, with a variance from its residuals. The
statistic of step 4 is the t statistic of the slope, or the square root of the sequential F
statistic of the groups.

The likelihood needs little work because the inverse of the correlation matrix of a sample
only connects neighbouring rows. If every sample covers every locus, the locus effects can
be profiled out in closed form, and an evaluation costs O(L), or O(1) for AR(1). Otherwise
the normal equations of the locus effects are sparse and are factorized on their envelope.

The C code doesn't use code of locfit, nlme or outliers: it is written from the published
descriptions, and checked against these packages (`tests/`, `inst/validation/`).

## References

Korthauer K, Chakraborty S, Benjamini Y, Irizarry RA (2019). Detection and accurate false
discovery rate control of differentially methylated regions from whole genome bisulfite
sequencing. *Biostatistics* 20(3), 367-383.

Loader C (1999). *Local Regression and Likelihood*. Springer.

## License

MIT. Parts of the R code are adapted from dmrseq (MIT License).
