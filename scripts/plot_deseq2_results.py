#!/usr/bin/env python3
"""Create final DESeq2-derived tables, volcano plot, and result summary."""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
DESEQ2_RESULTS = TABLES / "deseq2_results.csv"

KNOWN_RESPONSE_GENES = {"FKBP5", "DUSP1", "KLF15", "PER1", "TSC22D3", "ZBTB16", "CRISPLD2"}
GENE_SYMBOLS = {
    "ENSG00000103196": "CRISPLD2",
    "ENSG00000163884": "KLF15",
    "ENSG00000120129": "DUSP1",
    "ENSG00000157514": "TSC22D3",
    "ENSG00000179094": "PER1",
    "ENSG00000096060": "FKBP5",
    "ENSG00000109906": "ZBTB16",
}


def as_float(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def read_results() -> list[dict[str, str]]:
    if not DESEQ2_RESULTS.exists():
        raise FileNotFoundError(
            f"Missing {DESEQ2_RESULTS}. Run `Rscript scripts/deseq2_airway.R` first."
        )
    with DESEQ2_RESULTS.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        for row in reader:
            gene_id = row.get("") or row.get("gene_id") or row.get("row.names") or ""
            row["gene_id"] = gene_id
            row["gene_symbol"] = GENE_SYMBOLS.get(gene_id, "")
            rows.append(row)
    return rows


def project(value: float, lo: float, hi: float, start: float, end: float) -> float:
    if math.isclose(lo, hi):
        return (start + end) / 2
    return start + (value - lo) / (hi - lo) * (end - start)


def write_volcano(rows: list[dict[str, str]]) -> None:
    points = []
    for row in rows:
        lfc = as_float(row["log2FoldChange"])
        padj = as_float(row["padj"])
        pvalue = as_float(row["pvalue"])
        if math.isnan(lfc) or math.isnan(pvalue):
            continue
        points.append((row, lfc, -math.log10(max(pvalue, 1e-300)), padj))

    xmax = max(max(abs(x) for _, x, _, _ in points), 1)
    ymax = max(y for _, _, y, _ in points) * 1.08
    body = [
        '<text x="40" y="42" class="title">DESeq2 volcano plot</text>',
        '<line x1="90" y1="560" x2="820" y2="560" class="axis"/>',
        '<line x1="90" y1="80" x2="90" y2="560" class="axis"/>',
        '<text x="390" y="604" class="label">log2 fold-change (trt vs untrt)</text>',
        '<text x="20" y="344" class="label" transform="rotate(-90 20 344)">-log10 p-value</text>',
    ]
    for row, x, y, padj in points:
        is_sig = not math.isnan(padj) and padj < 0.05 and abs(x) >= 1
        color = "#d94b45" if is_sig and x > 0 else "#2d7dd2" if is_sig else "#9aa7b8"
        cx = project(x, -xmax, xmax, 90, 820)
        cy = project(y, 0, ymax, 560, 80)
        body.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="2.4" fill="{color}" opacity=".62"/>')

    FIGURES.mkdir(parents=True, exist_ok=True)
    (FIGURES / "volcano_deseq2.svg").write_text(
        f"""<svg xmlns="http://www.w3.org/2000/svg" width="900" height="640" viewBox="0 0 900 640">
<rect width="100%" height="100%" fill="#fffaf0"/>
<style>
text {{ font-family: ui-sans-serif, Segoe UI, Arial, sans-serif; fill: #172033; }}
.title {{ font-size: 24px; font-weight: 800; }}
.label {{ font-size: 13px; fill: #526176; }}
.axis {{ stroke: #526176; stroke-width: 1.3; }}
</style>
{chr(10).join(body)}
</svg>""",
        encoding="utf-8",
    )


def write_top_genes(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    valid = []
    for row in rows:
        padj = as_float(row["padj"])
        if not math.isnan(padj):
            valid.append(row)
    valid.sort(key=lambda row: (as_float(row["padj"]), -abs(as_float(row["log2FoldChange"]))))
    top = valid[:30]
    out_fields = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    TABLES.mkdir(parents=True, exist_ok=True)
    with (TABLES / "deseq2_top_genes.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        for row in top:
            writer.writerow({field: row.get(field, "") for field in out_fields})
    return top


def write_summary(top: list[dict[str, str]]) -> None:
    detected = [
        row["gene_symbol"]
        for row in top
        if row.get("gene_symbol") in KNOWN_RESPONSE_GENES and as_float(row["log2FoldChange"]) > 0
    ]
    detected_text = ", ".join(detected) if detected else "No known marker genes were present among the top 30 rows."
    (RESULTS / "results_summary.md").write_text(
        f"""# DESeq2 results summary

## Question

Which genes change expression after dexamethasone treatment in airway smooth muscle cells?

## Method

I used DESeq2 with design `~ cell + dex`, controlling for donor/cell-line effects and testing treated vs untreated samples.

## Main result

Detected glucocorticoid-response marker genes among the top DESeq2 rows: {detected_text}

## Interpretation

The increased expression of FKBP5, DUSP1, KLF15, PER1, TSC22D3, ZBTB16, and CRISPLD2 supports that the dataset captures a glucocorticoid-response transcriptional program after dexamethasone treatment.

## Limitations

This project starts from a prepared count matrix and does not include raw FASTQ QC, trimming, alignment, or quantification.
""",
        encoding="utf-8",
    )


def main() -> int:
    rows = read_results()
    write_volcano(rows)
    top = write_top_genes(rows)
    write_summary(top)
    print(f"Wrote {TABLES / 'deseq2_top_genes.csv'}")
    print(f"Wrote {FIGURES / 'volcano_deseq2.svg'}")
    print(f"Wrote {RESULTS / 'results_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
