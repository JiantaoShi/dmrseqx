# dmrseqx() against a complete run of unmodified dmrseq 1.32.0: chr21 of the ESCC data, 3 vs 3 samples,
# so that all 9 relabellings are used. Usage: Rscript check_original.R <bsseq.rds> <results of dmrseq> [lib]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({ library(dmrseqx); library(bsseq); library(data.table); library(GenomicRanges) })
options(width = 170, digits = 4)
bs <- readRDS(args[1])
bs <- bs[as.character(seqnames(bs)) == "chr21"]
control <- c("ctrl_N1", "ctrl_N2", "ctrl_N6"); case <- c("case_T1", "case_T2", "case_T6")
bs <- bs[, c(control, case)]
pData(bs)$group <- factor(c(rep("control", 3), rep("case", 3)), levels = c("control", "case"))
cov <- as.matrix(getCoverage(bs, type = "Cov"))
bs <- bs[rowSums(cov == 0) == 0]
cat("CpGs", nrow(bs), "\n")
set.seed(1)
t0 <- Sys.time()
res <- dmrseqx(bs, testCovariate = "group", cutoff = 0.1, threads = 1L, verbose = FALSE)
cat("dmrseqx", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), "s, permutations", length(metadata(res)$dmrseqx$labels) - 1, "\n")
A <- as.data.table(as.data.frame(res))[, .(chr = as.character(seqnames), start, end, L, area, beta, stat, pval, qval)]
ref <- fread(file.path(args[2], "regions_ranked.tsv.gz"))
print(fread(file.path(args[2], "runtime.tsv")))
m <- merge(ref, A, by = c("chr", "start", "end"))
cat("regions dmrseq", nrow(ref), "dmrseqx", nrow(A), "same coordinates", nrow(m), "\n")
m[, dt := abs(stat.y - stat.x)]
print(m[, .(n = .N, max_dt = max(dt, na.rm = TRUE), within_1e3 = mean(dt < 1e-3, na.rm = TRUE), max_dp = max(abs(pval.x - pval.y), na.rm = TRUE), max_dq = max(abs(padj - qval), na.rm = TRUE)), keyby = .(long = L >= 40)])
a <- m$padj <= 0.05; b <- m$qval <= 0.05
cat(sprintf("calls q <= 0.05: dmrseq %d, dmrseqx %d, both %d, differ %.2f%%\n", sum(a), sum(b), sum(a & b), 100 * sum(a != b) / sum(a | b)))
cat("ORIGINAL_DONE\n")
