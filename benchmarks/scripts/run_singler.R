# Run SingleR trained on v5 split A and classify split B.
#
# Purpose: the SingleR comparator for protocol v5, with the confidence and abstention
#          rules the protocol declares.
# Run:     Rscript benchmarks/scripts/run_singler.R SPLITS_DIR OUT_TSV
# Inputs:
#     SPLITS_DIR/A/{lognorm.mtx,genes.tsv,barcodes.tsv,labels.tsv}  [export_splits.py]
#     SPLITS_DIR/B/{lognorm.mtx,genes.tsv,barcodes.tsv}             [export_splits.py]
# Outputs:
#     OUT_TSV   barcode, predicted_label, uncertainty
#               Pruned calls are Unknown; uncertainty = 1 - delta.next / max(delta.next on B).

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleR)
})

args <- commandArgs(trailingOnly = TRUE)
SPLITS <- args[1]
OUT <- args[2]
stopifnot(!file.exists(OUT))

# Read splits --------------------------------------------------------------------

read_split <- function(name) {
  folder <- file.path(SPLITS, name)
  values <- as(readMM(file.path(folder, "lognorm.mtx")), "CsparseMatrix")
  dimnames(values) <- list(
    readLines(file.path(folder, "genes.tsv")),
    readLines(file.path(folder, "barcodes.tsv"))
  )
  values
}

train <- read_split("A")
test <- read_split("B")
labels <- read.delim(file.path(SPLITS, "A", "labels.tsv"))
stopifnot(identical(labels$barcode, colnames(train)))

# Train and classify ------------------------------------------------------------------

set.seed(0)
trained <- trainSingleR(train, labels = labels$cell_type, de.method = "wilcox")
result <- classifySingleR(test, trained)

# Write ------------------------------------------------------------------------------

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
