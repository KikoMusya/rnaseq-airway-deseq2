# DESeq2 results summary

## Question

Which genes change expression after dexamethasone treatment in airway smooth muscle cells?

## Method

I used DESeq2 with design `~ cell + dex`, controlling for donor/cell-line effects and testing treated vs untreated samples.

## Main result

Detected glucocorticoid-response marker genes among the top DESeq2 rows: DUSP1, PER1, KLF15

## Interpretation

The increased expression of FKBP5, DUSP1, KLF15, PER1, TSC22D3, ZBTB16, and CRISPLD2 supports that the dataset captures a glucocorticoid-response transcriptional program after dexamethasone treatment.

## Limitations

This project starts from a prepared count matrix and does not include raw FASTQ QC, trimming, alignment, or quantification.
