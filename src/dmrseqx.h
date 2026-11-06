#ifndef DMRSEQX_H
#define DMRSEQX_H

#include <stddef.h>
#include <stdint.h>

/* No fused multiply-add contraction, so that results do not depend on the
 * compiler and CPU more than they have to. */
#ifdef __clang__
#pragma clang fp contract(off)
#endif

/* Everything in this header runs inside worker threads and must not touch the
 * R API. */

/* Counts of all CpGs: nTotal x nSamples, column-major. */
typedef struct {
    const int *M, *Cov, *pos;
    size_t nTotal;
    int nSamples;
} dataset_t;

typedef struct {
    double cutoff, bpSpan, minInSpan, bpSpan2, maxGapSmooth, maxGap, covQ75;
    int minNumRegion, smooth, fitGlobal;
    int twoGroup;         /* two groups: pooled effects, refit of extreme coefficients */
    int factor;           /* the covariate of interest is a factor */
    int nTest;            /* columns of the covariate of interest */
    int nDesign;          /* columns of the design: intercept, covariate of interest, others */
    const double *tCrit;  /* qt(0.995, df) for df = 1, ..., nTCrit */
    int nTCrit;
    double grubbsT[11];   /* qt(1 - 0.01 / n, n - 2) for n = 3, ..., 10 */
} params_t;

/* ------------------------------------------------------------------ random */

typedef struct { uint64_t s[2]; } rng_t;

void rng_seed(rng_t *r, uint64_t seed, uint64_t a, uint64_t b);
double rng_unif(rng_t *r);

/* ---------------------------------------------------------------- smoothing */

/* Local quadratic regression of one smoothing block on an adaptive tree of
 * vertices. xFit is sorted, yFit and wFit are finite. Fitted values are
 * written for the sorted positions xEval. Returns 0 if out of memory. */
int smooth_block(const double *xFit, const double *yFit, const double *wFit, long nFit,
                 double nn, double h, double cut,
                 const double *xEval, long nEval, double *fitted);

/* --------------------------------------------------------------- region fit */

typedef struct fit_ws fit_ws_t;

fit_ws_t *fit_ws_new(void);
void fit_ws_free(fit_ws_t *ws);

typedef struct {
    double stat, phi;
    double *beta;  /* nTest coefficients, set by the caller */
    int constant;
} fit_out_t;

/* Fits one region: L rows of the dataset, rows offset + ix[0], ...; X holds the
 * covariates of each sample without intercept (nSamples x p), the P->nTest columns
 * of the covariate of interest first. With several columns the statistic is the
 * square root of the sequential F statistic of these columns. weighted = 0 is the
 * refit without variance weights. Returns 0 if out of memory. */
int fit_region(const dataset_t *ds, const params_t *P, size_t offset, const int *ix, int L,
               const double *X, int p, int weighted, rng_t *rng, fit_ws_t *ws, fit_out_t *out);

/* ---------------------------------------------------------------- one task */

typedef struct {
    int start, end, indexStart, indexEnd, L;
    double area, stat;
} region_t;

typedef struct {
    region_t *rows;
    double *beta;  /* n x nTest, by rows */
    size_t n, cap;
} regions_t;

typedef struct chrom_ws chrom_ws_t;

chrom_ws_t *chrom_ws_new(int nSamples, int maxCpG, int nDesign, int nTest);
void chrom_ws_free(chrom_ws_t *ws);

/* Candidate regions and their fit on one chromosome (rows offset, ..., offset +
 * n - 1) for one labelling: design is the design matrix with intercept (nSamples x
 * P->nDesign). Returns 0 if out of memory. */
int process_chromosome(const dataset_t *ds, const params_t *P, size_t offset, int n,
                       const double *design, const double *X, int p, uint64_t seed, uint64_t pass,
                       uint64_t chrom, chrom_ws_t *ws, regions_t *out);

#endif
