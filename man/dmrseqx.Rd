\name{dmrseqx}
\alias{dmrseqx}
\title{Differentially methylated regions with the procedure of dmrseq, in compiled code}
\description{
Finds candidate regions from smoothed differences between two or more groups of
samples, or along a continuous covariate, tests each with a region level generalized
least squares model with correlated errors, and assesses significance with a pooled
permutation null. The steps and the results follow \code{dmrseq::dmrseq()}; the
computations are done in C and in threads.
}
\usage{
dmrseqx(bs, testCovariate, adjustCovariate = NULL, cutoff = 0.1, minNumRegion = 5,
        smooth = TRUE, bpSpan = 1000, minInSpan = 30, maxGapSmooth = 2500,
        maxGap = 1000, verbose = TRUE, maxPerms = 10, matchCovariate = NULL,
        stat = "stat", threads = 1L, fit = c("dmrseq", "global"), perms = NULL)
}
\arguments{
  \item{bs}{a sorted \code{BSseq} object with raw counts. The counts are read into
    memory: 8 bytes per CpG and sample.}
  \item{testCovariate}{name or column number in \code{pData(bs)} of the covariate
    of interest. With two values: two groups. A factor or character column with more
    values: a global test among several groups. A numeric column with more values: a
    continuous covariate. See Details for the last two.}
  \item{adjustCovariate}{optional names or column numbers of covariates to adjust for.}
  \item{cutoff}{cutoff on the smoothed difference for candidate regions, above 0 and
    at most 1. Its meaning with several groups and a continuous covariate is in Details.}
  \item{minNumRegion}{smallest number of CpGs in a region, at least 3.}
  \item{smooth}{logical: smooth the differences.}
  \item{bpSpan}{smoothing window in base pairs.}
  \item{minInSpan}{about how many CpGs the smoothing window contains at least.}
  \item{maxGapSmooth}{CpGs further apart than this are smoothed separately.}
  \item{maxGap}{CpGs further apart than this are not in the same region.}
  \item{verbose}{logical: report progress.}
  \item{maxPerms}{largest number of permutations of the labels.}
  \item{matchCovariate}{optional name or column number of a covariate, for two groups
    only. Permutations that do not have at least one value of it in both groups are
    skipped.}
  \item{stat}{statistic for significance: one of \code{"stat"}, \code{"L"},
    \code{"area"}, \code{"beta"} or \code{"avg"}. \code{"beta"} needs one coefficient,
    so several groups can't use it.}
  \item{threads}{number of threads. Results do not depend on it.}
  \item{fit}{how the correlation parameter of regions of less than 40 CpGs is found.
    \code{"dmrseq"} (default) reproduces what \code{nlme::gls()} does inside
    \code{dmrseq()}, see Details. \code{"global"} takes the global optimum of the
    restricted likelihood.}
  \item{perms}{optional matrix with one row per sample and one column per permutation,
    instead of drawing permutations. For two groups it holds the 0 or 1 of each sample,
    otherwise each column is a permutation of the samples.}
}
\details{
The statistical model and the defaults are those of \code{dmrseq()}. Candidate regions are identical up to rounding errors, including the interpolation
of the local regression. For regions of at least 40 CpGs, the test statistics agree
to about 1e-4.

Regions of less than 40 CpGs use a correlation that decays with the distance in
base pairs, and \code{gls()} starts its optimization at a correlation parameter of
0.2. When the CpGs of a region are about 8 bp or more apart, the likelihood is flat
there and the optimizer stops at once, so these regions keep the start value, which
means almost no correlation. In a small band of cases rounding errors decide whether
it stops: the same region can give a very different statistic after the samples are
reordered. \code{fit = "dmrseq"} applies a rule on the likelihood at the start value,
and the statistics of about 98\% of these regions agree with \code{dmrseq()} to 1e-3.
The remaining difference can't be removed. \code{fit = "global"} is the true restricted
maximum likelihood estimate; it gives smaller statistics for the regions that
\code{dmrseq()} leaves at the start value.

With several groups or a continuous covariate, the effect of each CpG doesn't come from
pooled proportions: the methylation proportions, with values of 0 and 1 moved inward by a
small value for each chromosome, are regressed on the design on the logit scale. As in
\code{dmrseq()}, where the code marks this path as experimental, what differs from two
groups is:
\itemize{
  \item With several groups the effect is the largest minus the smallest fitted proportion
    among the groups. It is never negative, so that regions have no direction, and noise
    alone can reach the default cutoff when groups are small: with 2 samples in each of 4
    groups, 86\% of 6,000 simulated CpGs fell into 10 regions. Use a larger cutoff then.
    \code{stat} is the square root of the sequential F statistic of the groups, and
    there is one coefficient column \code{beta_<level>} per level except the first.
  \item With a continuous covariate the effect is the change of the fitted proportion between
    the values 0 and 1 of the covariate. \code{cutoff} therefore depends on both the units and
    the origin of the covariate, and a value of 0 far outside the data means an effect at an
    extrapolated value. Center and scale the covariate so that this change is meaningful.
    For a given region, \code{stat} is the t statistic of the slope and doesn't depend on the
    units, and \code{beta} is the slope itself.
  \item Permutations reorder the samples without restriction, and adjustment covariates are
    reordered with them for the region model. There are no refits for extreme coefficients.
  \item A permutation can leave a group without coverage at a CpG, even if the observed labels
    don't. Such a CpG has no effect of its own in that pass, and the smoother passes over it.
    \code{dmrseq()} usually stops with an error here, and also when a chromosome has only one
    pattern of samples without coverage. To run both, keep CpGs that lack coverage in fewer
    samples than the smallest group has.
}

Random numbers come from the R generator: use \code{set.seed()} for reproducible
permutations and jitter. Blocks are not supported.
}
\value{
A \code{GRanges} object with the same columns as \code{dmrseq()}, sorted by p-value:
\code{L}, the number of CpGs; \code{area}, the absolute sum of the smoothed differences;
\code{beta}, the coefficient of the second level of the covariate of interest, or of the
continuous covariate, on the arcsine scale (with several groups \code{beta_<level>}, one
column per level except the first); \code{stat}, its test statistic; \code{pval}, the permutation p-value;
\code{qval}, the Benjamini-Hochberg q-value; and \code{index}, the rows of \code{bs}.
\code{NULL} if there are no candidate regions. \code{metadata(result)$dmrseqx} holds the
labels of the permutations and the seed of the run.
}
\references{
Korthauer K, Chakraborty S, Benjamini Y, Irizarry RA (2019). Detection and accurate
false discovery rate control of differentially methylated regions from whole genome
bisulfite sequencing. Biostatistics 20(3), 367-383.
}
\examples{
\dontrun{
regions <- dmrseqx(bs, testCovariate = "group", threads = 8)
regions[which(regions$qval < 0.05)]

# a global test among three or more groups: a larger cutoff if groups are small
regions <- dmrseqx(bs, testCovariate = "tissue", cutoff = 0.2, threads = 8)

# a continuous covariate: the effect is the change between the values 0 and 1
bsseq::pData(bs)$age10 <- (bsseq::pData(bs)$age - 50) / 10
regions <- dmrseqx(bs, testCovariate = "age10", cutoff = 0.05, threads = 8)
}
}
