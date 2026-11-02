#!/usr/bin/env python3
"""Build synchronized, evidence-bounded Supplementary Figures S1 and S2.

The western-blot composites are deliberately labelled as archival screening
records A and B. Surviving records do not permit a cell-line assignment for
those composites. Quantitative phenotype panels remain cell-line specific and
use independent culture/transfection experiments as the inferential unit.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


COLORS = {
    "NC": "#4A4A4A",
    "si-IFI44L-2": "#2A9D8F",
    "si-IFI44L-3": "#D95F4B",
}
MARKERS = {"NC": "o", "si-IFI44L-2": "s", "si-IFI44L-3": "^"}
ORDER = ["NC", "si-IFI44L-2", "si-IFI44L-3"]


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def panel_label(ax: plt.Axes, label: str, x: float = -0.14, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top")


def save_all(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(out_dir / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(
        out_dir / f"{stem}.tiff",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )


def mean_sd(ax: plt.Axes, x: float, values: np.ndarray, color: str) -> None:
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    ax.errorbar(
        x,
        mean,
        yerr=sd,
        fmt="_",
        ms=13,
        mew=1.3,
        color=color,
        ecolor=color,
        capsize=3,
        elinewidth=1,
        zorder=5,
    )


def paired_endpoint(ax: plt.Axes, df: pd.DataFrame, ylabel: str, title: str, normalize: bool = False) -> None:
    pivot = df.pivot(index="experiment_id", columns="treatment", values="mean_value")[ORDER]
    if normalize:
        pivot = pivot.div(pivot["NC"], axis=0)
    xs = np.arange(len(ORDER), dtype=float)
    for _, row in pivot.iterrows():
        ax.plot(xs, row.values, color="#B8B8B8", lw=0.8, zorder=1)
    for x, treatment in zip(xs, ORDER):
        values = pivot[treatment].to_numpy(float)
        jitter = np.linspace(-0.045, 0.045, len(values))
        ax.scatter(
            np.repeat(x, len(values)) + jitter,
            values,
            s=25,
            marker=MARKERS[treatment],
            color=COLORS[treatment],
            edgecolor="white",
            linewidth=0.5,
            zorder=4,
        )
        mean_sd(ax, x, values, COLORS[treatment])
    ax.set_xticks(xs, ORDER, rotation=25, ha="right", rotation_mode="anchor")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold", loc="center")
    ax.margins(x=0.16)


def time_course(ax: plt.Axes, df: pd.DataFrame, ylabel: str, title: str) -> None:
    time_order = ["24h", "48h", "72h"] if "72h" in set(df["endpoint"]) else ["24h", "48h"]
    xs = np.arange(len(time_order), dtype=float)
    for treatment in ORDER:
        sub = df[df["treatment"] == treatment]
        pivot = sub.pivot(index="experiment_id", columns="endpoint", values="mean_value")[time_order]
        for _, row in pivot.iterrows():
            ax.plot(xs, row.values, color=COLORS[treatment], alpha=0.22, lw=0.7, zorder=1)
        means = pivot.mean(axis=0).to_numpy(float)
        sds = pivot.std(axis=0, ddof=1).to_numpy(float)
        ax.errorbar(
            xs,
            means,
            yerr=sds,
            marker=MARKERS[treatment],
            color=COLORS[treatment],
            label=treatment,
            lw=1.4,
            ms=4.5,
            capsize=2.5,
            zorder=3,
        )
        for index, endpoint in enumerate(time_order):
            values = pivot[endpoint].to_numpy(float)
            ax.scatter(
                np.repeat(xs[index], len(values)),
                values,
                s=11,
                marker=MARKERS[treatment],
                color=COLORS[treatment],
                alpha=0.72,
                zorder=2,
            )
    ax.set_xticks(xs, [value.replace("h", " h") for value in time_order])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold", loc="center")


def image_axis(ax: plt.Axes, path: Path, title: str | None = None) -> None:
    ax.imshow(Image.open(path).convert("RGB"))
    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=6.5, pad=2)


def condition_strip(
    spec,
    paths: list[Path],
    labels: list[str],
    panel: str,
    row_labels: list[str] | None = None,
) -> None:
    rows = len(row_labels) if row_labels else 1
    cols = len(paths) // rows
    grid = spec.subgridspec(rows, cols, wspace=0.025, hspace=0.035)
    for index, (path, label) in enumerate(zip(paths, labels)):
        ax = plt.subplot(grid[index // cols, index % cols])
        image_axis(ax, path, label if index < cols else None)
        if index == 0:
            panel_label(ax, panel, x=-0.18, y=1.17)
        if row_labels and index % cols == 0:
            ax.text(
                -0.08,
                0.5,
                row_labels[index // cols],
                transform=ax.transAxes,
                rotation=90,
                rotation_mode="anchor",
                ha="right",
                va="center",
                fontsize=6.5,
            )


def build_s1(data: pd.DataFrame, wb_a: Path, wb_b: Path, colony_dir: Path, out_dir: Path) -> None:
    fig = plt.figure(figsize=(7.09, 8.45), constrained_layout=False)
    outer = fig.add_gridspec(
        5,
        2,
        height_ratios=[1.0, 0.18, 1.05, 0.95, 1.1],
        hspace=0.40,
        wspace=0.34,
    )
    wb_paths = [wb_a, wb_b]
    for col, (title, wb_path) in enumerate(
        zip(["Archival siRNA-screening blot A", "Archival siRNA-screening blot B"], wb_paths)
    ):
        ax = fig.add_subplot(outer[0, col])
        image_axis(ax, wb_path)
        ax.set_title(title, fontsize=7.2, fontweight="bold", pad=3)
        panel_label(ax, ["a", "b"][col], x=-0.10, y=1.08)

    note_ax = fig.add_subplot(outer[1, :])
    note_ax.set_axis_off()
    note_ax.text(
        0.5,
        0.56,
        "Cell-line assignment for panels a-b was not recoverable from surviving records; qualitative documentation only.",
        ha="center",
        va="center",
        fontsize=6.7,
        color="#4A4A4A",
    )

    for col, cell_line in enumerate(["SUM159", "MDA-MB-231"]):
        ax = fig.add_subplot(outer[2, col])
        sub = data[(data["assay"] == "CCK-8") & (data["cell_line"] == cell_line)]
        time_course(ax, sub, "CCK-8 absorbance", cell_line)
        panel_label(ax, ["c", "d"][col])
        if col == 0:
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.0, 1.02),
                ncol=3,
                fontsize=6,
                handletextpad=0.4,
                columnspacing=0.9,
            )

        colony_names = (
            ["sum159-NC.png", "sum159-Si2.png", "sum159-si3.png"]
            if cell_line == "SUM159"
            else ["MDAMB231-NC.png", "MDAMB231-si2.png", "MDAMB231-si3.png"]
        )
        condition_strip(
            outer[3, col],
            [colony_dir / name for name in colony_names],
            ["NC", "si-IFI44L-2", "si-IFI44L-3"],
            ["e", "f"][col],
        )

        ax = fig.add_subplot(outer[4, col])
        sub = data[(data["assay"] == "Clonogenic assay") & (data["cell_line"] == cell_line)]
        paired_endpoint(ax, sub, "Relative colony formation", cell_line)
        panel_label(ax, ["g", "h"][col])

    fig.suptitle(
        "Archival IFI44L siRNA-screening blots and in vitro growth phenotypes",
        x=0.06,
        y=0.995,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    save_all(fig, out_dir, "Supplementary_Figure_S1_v8")
    plt.close(fig)


def build_s2(data: pd.DataFrame, wound_dir: Path, transwell_dir: Path, out_dir: Path) -> None:
    fig = plt.figure(figsize=(7.09, 8.6), constrained_layout=False)
    outer = fig.add_gridspec(4, 2, height_ratios=[1.65, 1.05, 1.05, 1.1], hspace=0.50, wspace=0.34)
    for col, cell_line in enumerate(["SUM159", "MDA-MB-231"]):
        prefix = "sum159" if cell_line == "SUM159" else "MDAMB231"
        wound_paths = [
            wound_dir / f"{prefix}-{time}-{group}.png"
            for time in ["0h", "48h"]
            for group in ["NC", "si2", "si3"]
        ]
        condition_strip(
            outer[0, col],
            wound_paths,
            ["NC", "si-IFI44L-2", "si-IFI44L-3"] * 2,
            ["a", "b"][col],
            row_labels=["0 h", "48 h"],
        )

        ax = fig.add_subplot(outer[1, col])
        sub = data[(data["assay"] == "Wound-healing assay") & (data["cell_line"] == cell_line)]
        time_course(ax, sub, "Wound closure (%)", cell_line)
        panel_label(ax, ["c", "d"][col])
        if col == 0:
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.0, 1.02),
                ncol=3,
                fontsize=6,
                handletextpad=0.4,
                columnspacing=0.9,
            )

        transwell_names = [f"{prefix}-{group}.png" for group in ["NC", "si2", "si3"]]
        condition_strip(
            outer[2, col],
            [transwell_dir / name for name in transwell_names],
            ["NC", "si-IFI44L-2", "si-IFI44L-3"],
            ["e", "f"][col],
        )

        ax = fig.add_subplot(outer[3, col])
        sub = data[(data["assay"] == "Transwell invasion") & (data["cell_line"] == cell_line)]
        paired_endpoint(ax, sub, "Invasive cells per field", cell_line)
        panel_label(ax, ["g", "h"][col])

    fig.suptitle(
        "Representative migration and invasion fields with experiment-level quantification",
        x=0.06,
        y=0.995,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    save_all(fig, out_dir, "Supplementary_Figure_S2_v8")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument(
        "--wb-a",
        required=True,
        type=Path,
        help="Path to archival siRNA-screening blot composite A; no cell-line assignment is implied.",
    )
    parser.add_argument(
        "--wb-b",
        required=True,
        type=Path,
        help="Path to archival siRNA-screening blot composite B; no cell-line assignment is implied.",
    )
    parser.add_argument("--colony-dir", required=True, type=Path)
    parser.add_argument("--wound-dir", required=True, type=Path)
    parser.add_argument("--transwell-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = pd.read_csv(args.data, sep="\t")
    build_s1(data, args.wb_a, args.wb_b, args.colony_dir, args.out_dir)
    build_s2(data, args.wound_dir, args.transwell_dir, args.out_dir)
    print(args.out_dir.resolve())


if __name__ == "__main__":
    main()
