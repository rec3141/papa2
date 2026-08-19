#!/usr/bin/env Rscript
# Per-step benchmark: R dada2, multithreaded. Args: raw_dir out_dir tf tr eeF eeR silva
suppressMessages(library(dada2))
args <- commandArgs(trailingOnly=TRUE)
raw <- args[1]; out <- args[2]
tf <- as.integer(args[3]); tr <- as.integer(args[4])
eeF <- as.numeric(args[5]); eeR <- as.numeric(args[6])
silva <- args[7]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
filt <- file.path(out, "filtered"); dir.create(filt, showWarnings=FALSE)

accs <- readLines(file.path(raw, "accs.txt"))
fnFs <- file.path(raw, paste0(accs, "_1.fastq.gz"))
fnRs <- file.path(raw, paste0(accs, "_2.fastq.gz"))
filtFs <- file.path(filt, paste0(accs, "_F.fastq.gz"))
filtRs <- file.path(filt, paste0(accs, "_R.fastq.gz"))

tt <- function(name, expr) {
  t0 <- proc.time()[["elapsed"]]; res <- expr
  cat(sprintf("TIME %s %.2f\n", name, proc.time()[["elapsed"]]-t0)); flush.console(); res
}

ft <- tt("filter", filterAndTrim(fnFs, filtFs, fnRs, filtRs,
        truncLen=c(tf,tr), maxEE=c(eeF,eeR), truncQ=2, rm.phix=TRUE,
        compress=TRUE, multithread=TRUE))
write.csv(ft, file.path(out, "filter_result.csv"))
keep <- ft[,2] > 0
accs <- accs[keep]; filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
writeLines(accs, file.path(out, "kept.txt"))

derepFs <- tt("derep", { a <- lapply(filtFs, derepFastq); b <- lapply(filtRs, derepFastq); list(a,b) })
derepRs <- derepFs[[2]]; derepFs <- derepFs[[1]]

set.seed(42)
errF <- tt("learnErrorsF", learnErrors(filtFs, multithread=TRUE, randomize=FALSE))
set.seed(42)
errR <- tt("learnErrorsR", learnErrors(filtRs, multithread=TRUE, randomize=FALSE))
write.csv(errF$err_out, file.path(out, "err_fwd.csv"))

set.seed(42)
dadaFs <- tt("dadaF", dada(derepFs, err=errF, multithread=TRUE))
set.seed(42)
dadaRs <- tt("dadaR", dada(derepRs, err=errR, multithread=TRUE))
if (inherits(dadaFs, "dada")) { dadaFs <- list(dadaFs); dadaRs <- list(dadaRs) }

mergers <- tt("merge", mergePairs(dadaFs, derepFs, dadaRs, derepRs))
if (is.data.frame(mergers)) mergers <- list(mergers)
seqtab <- tt("seqtab", makeSequenceTable(mergers))
rownames(seqtab) <- accs
nochim <- tt("chimera", removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE))
cat("dims:", dim(seqtab), "->", dim(nochim), "\n")
write.csv(seqtab, file.path(out, "seqtab.csv"))
write.csv(nochim, file.path(out, "seqtab_nochim.csv"))
for (i in seq_along(accs)) writeLines(dadaFs[[i]]$sequence, file.path(out, paste0("dada_", accs[i], "_F.txt")))

set.seed(100)
taxa <- tt("taxonomy", assignTaxonomy(colnames(nochim), silva, multithread=TRUE,
           taxLevels=c("Kingdom","Phylum","Class","Order","Family","Genus")))
write.csv(taxa, file.path(out, "taxonomy.csv"))
