#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
local_lib <- file.path(root, ".r-lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

# Formal DESeq2 analysis for the airway RNA-seq mini-project.
# Requires:
#   Rscript scripts/install_r_dependencies.R

suppressPackageStartupMessages({
  library(DESeq2)
})

counts_path <- file.path(root, "data", "raw", "airway_counts.csv")
metadata_path <- file.path(root, "data", "raw", "airway_colData.csv")
out_dir <- file.path(root, "results", "tables")
results_dir <- file.path(root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

counts <- counts[, rownames(coldata)]
coldata$dex <- relevel(factor(coldata$dex), ref = "untrt")
coldata$cell <- factor(coldata$cell)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = coldata,
  design = ~ cell + dex
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("dex", "trt", "untrt"))
res <- res[order(res$padj), ]

write.csv(as.data.frame(res), file.path(out_dir, "deseq2_results.csv"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))

message("Wrote: ", file.path(out_dir, "deseq2_results.csv"))
message("Wrote: ", file.path(results_dir, "sessionInfo.txt"))
