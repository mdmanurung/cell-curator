# Train SingleR on the v8 training cells and classify a v8 query.
#
# Purpose: the SingleR comparator for protocol v8 (Q6), with v6's settings and rules:
#          wilcox marker detection, pruned calls are Unknown, uncertainty =
#          1 - delta.next / max(delta.next over the query). SingleR can abstain: pruning
#          sets pruned.labels to NA when a cell's best score is not clearly above the rest.
# Run:     Rscript benchmarks/scripts/v8_singler.R TRAIN_EXPORT QUERY_EXPORT OUT_TSV [--check]
#          QUERY_EXPORT and OUT_TSV may be comma-separated lists of equal length: the
#          training runs once and each query is classified in turn.
#          (sbatch benchmarks/scripts/v8_singler.sbatch)
# Inputs:
#     TRAIN_EXPORT, QUERY_EXPORT  v8_export_for_r.py folders; TRAIN_EXPORT has labels.tsv
# Outputs:
#     OUT_TSV   barcode, predicted_label, uncertainty
#     OUT_TSV.sessionInfo.txt

# Parameters --------------------------------------------------------------------------
DE_METHOD <- "wilcox" # v6 (v5 route): trainSingleR(de.method = "wilcox")
SEED <- 0 # v6
N_CORES <- 4 # speed only; the result does not change (4 in the 4-CPU allocation)

args <- commandArgs(trailingOnly = TRUE)
TRAIN <- args[1]
QUERY <- strsplit(args[2], ",")[[1]]
OUT <- strsplit(args[3], ",")[[1]]
stopifnot(
  "train export missing" = dir.exists(TRAIN),
  "query export missing" = all(dir.exists(QUERY)),
  "one output per query" = length(QUERY) == length(OUT),
  "output exists" = !any(file.exists(OUT))
)
if ("--check" %in% args) {
  cat("check ok:", TRAIN, QUERY, OUT, "\n")
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleR)
  library(BiocParallel)
})

# Read the exports ---------------------------------------------------------------------

read_export <- function(folder) {
  size <- as.integer(readLines(file.path(folder, "dim.txt")))
  p <- readBin(file.path(folder, "p.bin"), "integer", n = size[2] + 1, size = 4)
  i <- readBin(
    file.path(folder, "i.bin"),
    "integer",
    n = p[size[2] + 1],
    size = 4
  )
  x <- readBin(
    file.path(folder, "x.bin"),
    "numeric",
    n = p[size[2] + 1],
    size = 4
  )
  new(
    "dgCMatrix",
    i = i,
    p = p,
    x = x,
    Dim = size,
    Dimnames = list(
      readLines(file.path(folder, "genes.tsv")),
      readLines(file.path(folder, "barcodes.tsv"))
    )
  )
}

train <- read_export(TRAIN)
labels <- read.delim(file.path(TRAIN, "labels.tsv"))
stopifnot(
  "labels match training cells" = identical(labels$barcode, colnames(train))
)
queries <- lapply(QUERY, read_export)

# Set up -------------------------------------------------------------------------------

# Each query is classified on the genes it shares with the training set (as with one
# query); queries with the same shared genes reuse one training.
set.seed(SEED)
param <- MulticoreParam(N_CORES)
trained_for <- list()

# Train, classify and write ------------------------------------------------------------

for (k in seq_along(queries)) {
  shared <- intersect(rownames(train), rownames(queries[[k]]))
  key <- paste(length(shared), digest::digest(shared))
  if (is.null(trained_for[[key]])) {
    trained_for[[key]] <- trainSingleR(
      train[shared, ],
      labels = labels$cell_type,
      de.method = DE_METHOD,
      BPPARAM = param
    )
  }
  result <- classifySingleR(
    queries[[k]][shared, ],
    trained_for[[key]],
    BPPARAM = param
  )
  predicted <- ifelse(
    is.na(result$pruned.labels),
    "Unknown",
    result$pruned.labels
  )
  delta <- result$delta.next
  write.table(
    data.frame(
      barcode = rownames(result),
      predicted_label = predicted,
      uncertainty = 1 - delta / max(delta)
    ),
    OUT[k],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
writeLines(capture.output(sessionInfo()), paste0(OUT[1], ".sessionInfo.txt"))
