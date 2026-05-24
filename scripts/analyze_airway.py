#!/usr/bin/env python3
"""Exploratory RNA-seq analysis from the prepared airway count matrix.

Creates quick exploratory outputs. Formal differential-expression claims
should use scripts/deseq2_airway.R.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import median

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data" / "raw" / "airway_counts.csv"
METADATA = ROOT / "data" / "raw" / "airway_colData.csv"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

GENE_SYMBOLS = {
    "ENSG00000103196": "CRISPLD2",
    "ENSG00000163884": "KLF15",
    "ENSG00000120129": "DUSP1",
    "ENSG00000157514": "TSC22D3",
    "ENSG00000179094": "PER1",
    "ENSG00000096060": "FKBP5",
    "ENSG00000109906": "ZBTB16",
    "ENSG00000152583": "SPARCL1",
}


def read_counts() -> tuple[list[str], list[str], np.ndarray]:
    with COUNTS.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        genes = []
        values = []
        for row in reader:
            genes.append(row[0])
            values.append([float(x) for x in row[1:]])
    return genes, samples, np.array(values, dtype=float)


def read_metadata() -> dict[str, dict[str, str]]:
    with METADATA.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        out = {}
        for row in reader:
            sample = row.get("") or row.get("Run")
            if sample:
                out[sample] = row
        return out


def count_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as f:
        return max(0, sum(1 for _ in f) - 1)


def size_factors(counts: np.ndarray) -> np.ndarray:
    valid = (counts > 0).all(axis=1)
    safe = counts[valid]
    geo = np.exp(np.mean(np.log(safe), axis=1))
    ratios = safe / geo[:, None]
    return np.array([median(ratios[:, i]) for i in range(ratios.shape[1])])


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda i: p_values[i], reverse=True)
    adjusted = [1.0] * n
    running = 1.0
    for rank_from_end, idx in enumerate(order, start=1):
        rank = n - rank_from_end + 1
        running = min(running, p_values[idx] * n / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def write_svg(path: Path, body: str, width: int = 900, height: int = 620) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#fffaf0"/>
<style>
text {{ font-family: ui-sans-serif, Segoe UI, Arial, sans-serif; fill: #172033; }}
.title {{ font-size: 24px; font-weight: 800; }}
.label {{ font-size: 13px; fill: #526176; }}
.small {{ font-size: 11px; fill: #334155; }}
.axis {{ stroke: #526176; stroke-width: 1.3; }}
</style>
{body}
</svg>""",
        encoding="utf-8",
    )


def project(value: float, lo: float, hi: float, start: float, end: float) -> float:
    if math.isclose(lo, hi):
        return (start + end) / 2
    return start + (value - lo) / (hi - lo) * (end - start)


def pca_plot(log_counts: np.ndarray, samples: list[str], metadata: dict[str, dict[str, str]]) -> None:
    x = log_counts.T
    x = x - x.mean(axis=0, keepdims=True)
    _, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = x @ vt[:2].T
    explained = (s ** 2) / max(1, x.shape[0] - 1)
    explained = explained / explained.sum()
    xmin, xmax = scores[:, 0].min(), scores[:, 0].max()
    ymin, ymax = scores[:, 1].min(), scores[:, 1].max()
    body = [
        '<text x="40" y="42" class="title">PCA: airway samples</text>',
        '<line x1="90" y1="540" x2="820" y2="540" class="axis"/>',
        '<line x1="90" y1="80" x2="90" y2="540" class="axis"/>',
        f'<text x="390" y="590" class="label">PC1 ({explained[0] * 100:.1f}% variance)</text>',
        f'<text x="20" y="330" class="label" transform="rotate(-90 20 330)">PC2 ({explained[1] * 100:.1f}% variance)</text>',
    ]
    for sample, (pc1, pc2) in zip(samples, scores):
        dex = metadata[sample]["dex"]
        color = "#d94b45" if dex == "trt" else "#2d7dd2"
        cx = project(float(pc1), xmin, xmax, 130, 780)
        cy = project(float(pc2), ymin, ymax, 500, 110)
        body.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="10" fill="{color}" opacity=".88"/>')
        body.append(f'<text x="{cx + 13:.1f}" y="{cy + 4:.1f}" class="small">{sample} ({dex})</text>')
    write_svg(FIGURES / "pca_airway.svg", "\n".join(body))


def volcano_plot(rows: list[dict[str, str]]) -> None:
    xs = [float(r["log2FC_trt_vs_untrt"]) for r in rows]
    ys = [-math.log10(max(float(r["p_value"]), 1e-300)) for r in rows]
    xmax = max(max(abs(x) for x in xs), 1)
    ymax = max(ys) * 1.08
    body = [
        '<text x="40" y="42" class="title">Approximate exploratory volcano plot</text>',
        '<line x1="90" y1="560" x2="820" y2="560" class="axis"/>',
        '<line x1="90" y1="80" x2="90" y2="560" class="axis"/>',
        '<text x="390" y="600" class="label">log2 fold-change</text>',
        '<text x="20" y="340" class="label" transform="rotate(-90 20 340)">-log10 p-value</text>',
    ]
    for row, x, y in zip(rows, xs, ys):
        padj = float(row["padj_bh"])
        color = "#d94b45" if padj < 0.05 and x > 1 else "#2d7dd2" if padj < 0.05 and x < -1 else "#9aa7b8"
        cx = project(x, -xmax, xmax, 90, 820)
        cy = project(y, 0, ymax, 560, 80)
        body.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="2.2" fill="{color}" opacity=".58"/>')
    write_svg(FIGURES / "volcano_airway.svg", "\n".join(body), 900, 640)


def heatmap(log_counts: np.ndarray, genes: list[str], samples: list[str], metadata: dict[str, dict[str, str]], top_ids: list[str]) -> None:
    top_idx = [genes.index(g) for g in top_ids if g in genes][:24]
    labels = [GENE_SYMBOLS.get(genes[i], genes[i].replace("ENSG00000", "E")) for i in top_idx]
    mat = log_counts[top_idx]
    z = (mat - mat.mean(axis=1, keepdims=True)) / np.where(mat.std(axis=1, keepdims=True) == 0, 1, mat.std(axis=1, keepdims=True))
    cell_w, cell_h = 70, 18
    body = ['<text x="34" y="40" class="title">Heatmap: top exploratory genes</text>']
    for j, sample in enumerate(samples):
        color = "#d94b45" if metadata[sample]["dex"] == "trt" else "#2d7dd2"
        body.append(f'<rect x="{170+j*cell_w}" y="56" width="{cell_w-2}" height="12" fill="{color}"/>')
        body.append(f'<text x="{170+j*cell_w+5}" y="76" class="small" transform="rotate(-35 {170+j*cell_w+5} 76)">{sample}</text>')
    for i, label in enumerate(labels):
        y = 90 + i * cell_h
        body.append(f'<text x="20" y="{y+12}" class="small">{label}</text>')
        for j in range(len(samples)):
            val = max(-2.5, min(2.5, float(z[i, j])))
            color = "#ffb38a" if val > 0 else "#9db7ff"
            body.append(f'<rect x="{170+j*cell_w}" y="{y}" width="{cell_w-2}" height="{cell_h-2}" fill="{color}" opacity="{0.35 + min(abs(val) / 2.5, 1) * 0.65:.2f}"/>')
    write_svg(FIGURES / "heatmap_top_genes.svg", "\n".join(body), 820, 620)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    genes, samples, counts = read_counts()
    metadata = read_metadata()
    keep_samples = [sample for sample in samples if sample in metadata]
    idx = [samples.index(sample) for sample in keep_samples]
    samples = keep_samples
    counts = counts[:, idx]
    keep = counts.sum(axis=1) >= 10
    genes = [gene for gene, keep_gene in zip(genes, keep) if keep_gene]
    counts = counts[keep]
    factors = size_factors(counts)
    normalized = counts / factors
    log_counts = np.log2(normalized + 1)
    trt = [i for i, sample in enumerate(samples) if metadata[sample]["dex"] == "trt"]
    ctl = [i for i, sample in enumerate(samples) if metadata[sample]["dex"] == "untrt"]

    rows = []
    p_values = []
    for i, gene in enumerate(genes):
        a = log_counts[i, trt]
        b = log_counts[i, ctl]
        lfc = float(a.mean() - b.mean())
        se = math.sqrt(float(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))
        z = lfc / se if se else 0.0
        p = math.erfc(abs(z) / math.sqrt(2))
        p_values.append(p)
        rows.append({
            "gene_id": gene,
            "gene_symbol": GENE_SYMBOLS.get(gene, ""),
            "baseMean": f"{float(normalized[i].mean()):.3f}",
            "log2FC_trt_vs_untrt": f"{lfc:.4f}",
            "p_value": f"{p:.6g}",
        })
    for row, padj in zip(rows, bh_adjust(p_values)):
        row["padj_bh"] = f"{padj:.6g}"
    rows.sort(key=lambda row: (float(row["padj_bh"]), -abs(float(row["log2FC_trt_vs_untrt"]))))

    with (TABLES / "de_results_approx.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with (TABLES / "normalized_counts_log2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id", "gene_symbol", *samples])
        for gene, values in zip(genes, log_counts):
            writer.writerow([gene, GENE_SYMBOLS.get(gene, ""), *[f"{v:.4f}" for v in values]])

    pca_plot(log_counts, samples, metadata)
    volcano_plot(rows)
    heatmap(log_counts, genes, samples, metadata, [row["gene_id"] for row in rows[:24]])
    (RESULTS / "analysis_summary.json").write_text(json.dumps({
        "dataset": "Bioconductor airway",
        "input": "prepared gene-level count matrix",
        "comparison": "dexamethasone treated vs untreated",
        "n_genes_in_matrix": count_rows(COUNTS),
        "n_genes_after_filter_sum_ge_10": len(genes),
        "n_samples": len(samples),
        "formal_de_note": "Run scripts/deseq2_airway.R before making final DE claims.",
    }, indent=2), encoding="utf-8")
    (RESULTS / "report.md").write_text(
        """# Exploratory RNA-seq report: airway dexamethasone response

This report summarizes exploratory Python results. Formal DESeq2 results are generated by `scripts/deseq2_airway.R` and visualized by `scripts/plot_deseq2_results.py`.

## Exploratory outputs

- `results/figures/pca_airway.svg`
- `results/figures/pca_airway.png`
- `results/figures/heatmap_top_genes.svg`
- `results/figures/heatmap_top_genes.png`
- `results/figures/volcano_airway.svg`
- `results/figures/volcano_airway.png`
- `results/tables/de_results_approx.csv`
- `results/tables/normalized_counts_log2.csv`

## Method note

The current report summarizes exploratory Python results. Formal DESeq2 results are generated by `scripts/deseq2_airway.R` and visualized by `scripts/plot_deseq2_results.py`.
""",
        encoding="utf-8",
    )
    print(f"Wrote exploratory outputs under {RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
