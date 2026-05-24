# Data Source

## Dataset

This project uses the public Bioconductor `airway` RNA-seq dataset.

The dataset contains RNA-seq measurements from human airway smooth muscle cells treated with dexamethasone (`trt`) and untreated controls (`untrt`). This mini-project starts from a prepared gene-level count matrix and sample metadata derived from that dataset.

## Analysis scope

This repository performs exploratory visualization and formal differential expression analysis from the prepared count matrix. It does not process raw FASTQ files, perform read quality control, trimming, alignment, or transcript quantification.

## Reproducibility note

The formal differential-expression analysis is implemented in `scripts/deseq2_airway.R` using DESeq2 with design `~ cell + dex`. The R runtime and package versions used for the generated results are recorded in `results/sessionInfo.txt`.

## References

- Bioconductor `airway` experiment data package
- DESeq2 documentation and Bioconductor RNA-seq workflows
