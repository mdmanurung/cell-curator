# Train SingleR on the v8 training cells and classify a v8 query.
#
# Purpose: the SingleR comparator for protocol v8 (Q6), with v6's settings and rules:
#          wilcox marker detection, pruned calls are Unknown, uncertainty =
#          1 - delta.next / max(delta.next over the query). SingleR can abstain: pruning
#          sets pruned.labels to NA when a cell's best score is not clearly above the rest.
# Run:     Rscript benchmarks/scripts/v8_singler.R TRAIN_EXPORT QUERY_EXPORT OUT_TSV [--check]
#          (sbatch benchmarks/scripts/v8_singler.sbatch)
# Inputs:
#     TRAIN_EXPORT, QUERY_EXPORT  v8_export_for_r.py folders; TRAIN_EXPORT has labels.tsv
# Outputs:
#     OUT_TSV   barcode, predicted_label, uncertainty
#     OUT_TSV.sessionInfo.txt

# Parameters --------------------------------------------------------------------------
DE_METHOD <- "wilcox" # v6 (v5 route): trainSingleR(de.method = "wilcox")
SEED <- 0 # v6
N_CORES <- 8 # speed only; the result does not change

args <- commandArgs(trailingOnly = TRUE)
TRAIN <- args[1]
QUERY <- args[2]
OUT <- args[3]
stopifnot(
  "train export missing" = dir.exists(TRAIN),
  "query export missing" = dir.exists(QUERY),
  "output exists" = !file.exists(OUT)
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
  i <- readBin(file.path(folder, "i.bin"), "integer", n = p[size[2] + 1], size = 4)
  x <- readBin(file.path(folder, "x.bin"), "numeric", n = p[size[2] + 1], size = 4)
  new(
    "dgCMatrix",
    i = i, p = p, x = x, Dim = size,
    Dimnames = list(
      readLines(file.path(folder, "genes.tsv")),
      readLines(file.path(folder, "barcodes.tsv"))
    )
  )
}

train <- read_export(TRAIN)
query <- read_export(QUERY)
labels <- read.delim(file.path(TRAIN, "labels.tsv"))
stopifnot("labels match training cells" = identical(labels$barcode, colnames(train)))
shared <- intersect(rownames(train), rownames(query))

# Train and classify -------------------------------------------------------------------

set.seed(SEED)
param <- MulticoreParam(N_CORES)
trained <- trainSingleR(
  train[shared, ], labels = labels$cell_type, de.method = DE_METHOD, BPPARAM = param
)
result <- classifySingleR(query[shared, ], trained, BPPARAM = param)

# Write --------------------------------------------------------------------------------

predicted <- ifelse(is.na(result$pruned.labels), "Unknown", result$pruned.labels)
delta <- result$delta.next
write.table(
  data.frame(
    barcode = rownames(result),
    predicted_label = predicted,
    uncertainty = 1 - delta / max(delta)
  ),
  OUT,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), paste0(OUT, ".sessionInfo.txt"))
