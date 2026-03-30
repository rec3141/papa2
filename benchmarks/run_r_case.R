args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    stop(paste("Missing required argument", flag))
  }
  args[[idx + 1]]
}

manifest_path <- get_arg("--manifest")
r_repo <- if ("--r-repo" %in% args) get_arg("--r-repo") else NA_character_

manifest <- jsonlite::fromJSON(manifest_path)

if (!is.na(r_repo) && nzchar(r_repo)) {
  configure_path <- file.path(r_repo, "configure")
  if (file.exists(configure_path)) {
    Sys.chmod(configure_path, mode = "0755")
  }
  pkgload::load_all(r_repo, quiet = TRUE, helpers = FALSE, export_all = FALSE)
} else {
  suppressPackageStartupMessages(library(dada2))
}

files <- manifest$files
learn_nbases <- as.integer(manifest$learn_nbases)

err <- learnErrors(files, nbases = learn_nbases, verbose = FALSE)
results <- lapply(files, function(f) {
  drp <- derepFastq(f, verbose = FALSE)
  dada(drp, err = err, verbose = FALSE)
})

n_asv_total <- sum(vapply(results, function(x) length(getUniques(x)), integer(1)))
n_reads_assigned <- sum(vapply(results, function(x) sum(getUniques(x)), numeric(1)))

cat(jsonlite::toJSON(
  list(
    engine = "r",
    case_id = manifest$case_id,
    sample_count = length(files),
    learn_nbases = learn_nbases,
    n_asv_total = unname(as.integer(n_asv_total)),
    n_reads_assigned = unname(as.integer(n_reads_assigned))
  ),
  auto_unbox = TRUE,
  pretty = TRUE
))
