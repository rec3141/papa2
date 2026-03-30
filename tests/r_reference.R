#!/usr/bin/env Rscript
# =============================================================================
# r_reference.R
# Run every major dada2 function on paired-end test data and persist all
# outputs so they can be used as ground-truth references for GPU reimplementation.
# =============================================================================

library(dada2)

set.seed(42)

# ---- paths ------------------------------------------------------------------
test_dir   <- "/data/papa2/tests"
data_dir   <- file.path(test_dir, "data")
out_dir    <- file.path(test_dir, "r_outputs")
filt_dir   <- file.path(out_dir, "filtered")

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
fnFs <- sort(list.files(data_dir, pattern = "F\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(data_dir, pattern = "R\\.fastq\\.gz$", full.names = TRUE))

sample_names <- sub("F\\.fastq\\.gz$", "", basename(fnFs))   # "sam1", "sam2"

filtFs <- file.path(filt_dir, paste0(sample_names, "F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample_names, "R_filt.fastq.gz"))

cat("Input forward files:", fnFs, "\n")
cat("Input reverse files:", fnRs, "\n")

# =============================================================================
# 1. filterAndTrim
# =============================================================================
cat("\n===== 1. filterAndTrim =====\n")
filt_result <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen  = c(240, 160),
  maxEE     = c(2, 2),
  truncQ    = 2,
  rm.phix   = TRUE,
  compress  = TRUE,
  multithread = FALSE
)
cat("Filter result:\n"); print(filt_result)
write.csv(filt_result, file.path(out_dir, "filter_result.csv"))

# =============================================================================
# 2. derepFastq
# =============================================================================
cat("\n===== 2. derepFastq =====\n")
derepFs <- lapply(filtFs, derepFastq)
derepRs <- lapply(filtRs, derepFastq)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

save_derep <- function(drp, tag) {
  seqs   <- drp$uniques |> names()
  abunds <- drp$uniques |> unname()
  quals  <- drp$quals

  writeLines(seqs,          file.path(out_dir, paste0("derep_", tag, "_seqs.txt")))
  writeLines(as.character(abunds), file.path(out_dir, paste0("derep_", tag, "_abunds.txt")))
  write.csv(quals,          file.path(out_dir, paste0("derep_", tag, "_quals.csv")))
  cat(sprintf("  %s: %d uniques, qual matrix %dx%d\n",
              tag, length(seqs), nrow(quals), ncol(quals)))
}

for (s in sample_names) {
  save_derep(derepFs[[s]], paste0(s, "F"))
  save_derep(derepRs[[s]], paste0(s, "R"))
}

# =============================================================================
# 3. learnErrors
# =============================================================================
cat("\n===== 3. learnErrors =====\n")
set.seed(42)
errF <- learnErrors(filtFs, multithread = FALSE, randomize = FALSE)
cat("Error matrix dimensions:", dim(errF$err_out), "\n")

# err_out is 16 rows (transitions) x ncol (quality scores)
err_mat <- errF$err_out
write.csv(err_mat, file.path(out_dir, "err_fwd.csv"))
saveRDS(errF, file.path(out_dir, "errF.rds"))

# =============================================================================
# 4. dada
# =============================================================================
cat("\n===== 4. dada =====\n")
set.seed(42)
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
# Also need reverse for mergePairs
set.seed(42)
errR <- learnErrors(filtRs, multithread = FALSE, randomize = FALSE)
set.seed(42)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

# If only one sample, dada returns a single dada object, not a list
if (inherits(dadaFs, "dada")) dadaFs <- list(dadaFs)
if (inherits(dadaRs, "dada")) dadaRs <- list(dadaRs)
names(dadaFs) <- sample_names
names(dadaRs) <- sample_names

for (s in sample_names) {
  dd <- dadaFs[[s]]
  cl_seqs   <- dd$sequence      # character vector of cluster seqs
  cl_abunds <- dd$denoised      # named integer vector
  cl_map    <- dd$map            # integer vector: unique-idx -> cluster-idx

  writeLines(cl_seqs,                     file.path(out_dir, paste0("dada_", s, "F_seqs.txt")))
  writeLines(as.character(cl_abunds),     file.path(out_dir, paste0("dada_", s, "F_abunds.txt")))
  writeLines(as.character(cl_map),        file.path(out_dir, paste0("dada_", s, "F_map.txt")))
  cat(sprintf("  %s: %d clusters, %d uniques mapped\n",
              s, length(cl_seqs), length(cl_map)))
}

# =============================================================================
# 5. mergePairs
# =============================================================================
cat("\n===== 5. mergePairs =====\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
if (is.data.frame(mergers)) mergers <- list(mergers)
names(mergers) <- sample_names

for (s in sample_names) {
  mg <- mergers[[s]]
  write.csv(mg[, c("sequence", "abundance", "forward", "reverse",
                    "nmatch", "nmismatch", "nindel", "accept")],
            file.path(out_dir, paste0("merge_", s, ".csv")),
            row.names = FALSE)
  cat(sprintf("  %s: %d merged pairs\n", s, nrow(mg)))
}

# =============================================================================
# 6. makeSequenceTable
# =============================================================================
cat("\n===== 6. makeSequenceTable =====\n")
seqtab <- makeSequenceTable(mergers)
cat("Sequence table dimensions:", dim(seqtab), "\n")
write.csv(seqtab, file.path(out_dir, "seqtab.csv"))

# =============================================================================
# 7. removeBimeraDenovo
# =============================================================================
cat("\n===== 7. removeBimeraDenovo =====\n")
set.seed(42)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
cat("After chimera removal:", dim(seqtab_nochim), "\n")
write.csv(seqtab_nochim, file.path(out_dir, "seqtab_nochim.csv"))

# =============================================================================
# 8. nwalign
# =============================================================================
cat("\n===== 8. nwalign =====\n")
aln <- nwalign("ACGTACGT", "ACGAACGT")
cat("Alignment:\n", aln[1], "\n", aln[2], "\n")
writeLines(aln, file.path(out_dir, "nwalign.txt"))

# =============================================================================
# 9. rc (reverse complement)
# =============================================================================
cat("\n===== 9. rc =====\n")
rc_result <- rc("ACGTACGTNNNN")
cat("RC:", rc_result, "\n")
writeLines(rc_result, file.path(out_dir, "rc.txt"))

# =============================================================================
# 10. nwhamming
# =============================================================================
cat("\n===== 10. nwhamming =====\n")
hamming <- nwhamming("ACGTACGT", "ACGAACGT")
cat("Hamming distance:", hamming, "\n")
writeLines(as.character(hamming), file.path(out_dir, "nwhamming.txt"))

# =============================================================================
# 11. isPhiX
# =============================================================================
cat("\n===== 11. isPhiX =====\n")
first5_seqs <- names(derepFs[[1]]$uniques)[1:5]
phix_result <- isPhiX(first5_seqs)
cat("isPhiX:", phix_result, "\n")
writeLines(as.character(phix_result), file.path(out_dir, "isphix.txt"))

# =============================================================================
# 12. seqComplexity
# =============================================================================
cat("\n===== 12. seqComplexity =====\n")
complexity <- seqComplexity(first5_seqs)
cat("Complexity:", complexity, "\n")
writeLines(as.character(complexity), file.path(out_dir, "complexity.txt"))

# =============================================================================
# 13. collapseNoMismatch
# =============================================================================
cat("\n===== 13. collapseNoMismatch =====\n")
collapsed <- collapseNoMismatch(seqtab_nochim)
cat("Collapsed table dimensions:", dim(collapsed), "\n")
write.csv(collapsed, file.path(out_dir, "collapse.csv"))

# =============================================================================
# 14. loessErrfun
# =============================================================================
cat("\n===== 14. loessErrfun =====\n")
# The dada object contains $trans — the observed transition counts matrix.
# loessErrfun expects a matrix with 16 rows (one per transition) and
# columns corresponding to quality scores.
trans_mat <- dadaFs[[1]]$trans
loess_err <- loessErrfun(trans_mat)
cat("loessErrfun output dimensions:", dim(loess_err), "\n")
write.csv(loess_err, file.path(out_dir, "loess_err.csv"))

# =============================================================================
# 15. isShiftDenovo
# =============================================================================
cat("\n===== 15. isShiftDenovo =====\n")

# isShiftDenovo expects a derep object (not a dada result).
# Per dada2 docs: it identifies sequences that are shifted versions of more
# abundant sequences in the dereplicated set.
shift_result <- tryCatch({
  res <- dada2:::isShiftDenovo(derepFs[[1]])
  cat("isShiftDenovo result (first 10):", head(res, 10), "\n")
  writeLines(as.character(res), file.path(out_dir, "isshift.txt"))
  res
}, error = function(e) {
  # Some versions expect a dada object or different input. Try alternatives.
  cat("  Trying with character vector of sequences...\n")
  seqs <- names(derepFs[[1]]$uniques)
  res <- dada2:::isShiftDenovo(seqs)
  cat("isShiftDenovo result (first 10):", head(res, 10), "\n")
  writeLines(as.character(res), file.path(out_dir, "isshift.txt"))
  res
})

# =============================================================================
cat("\n===== DONE =====\n")
cat("All outputs saved to:", out_dir, "\n")
