#!/usr/bin/env python3
"""Build the final spatial main figure with a prespecified representative section."""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import pandas as pd
from PIL import Image


BLUE = "#2B6CA3"
TEAL = "#23877D"
RED = "#B64E55"
GREY = "#747474"
LIGHT_GREY = "#D8D8D8"
BLACK = "#222222"


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6.6,
        "ytick.labelsize": 6.6,
        "legend.fontsize": 6.6,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def save_bundle(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_analysis_module(project: Path):
    path = project / "scripts" / "09_gse210616_spatial_analysis.py"
    spec = importlib.util.spec_from_file_location("spatial_analysis", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def spatial_panel(ax: plt.Axes, image: np.ndarray, frame: pd.DataFrame, scale: float,
                  column: str, title: str, cmap: str, scalebar_px: float) -> None:
    ax.imshow(image)
    x = frame["pxl_col_fullres"].to_numpy() * scale
    y = frame["pxl_row_fullres"].to_numpy() * scale
    values = frame[column].to_numpy(float)
    vmax = np.nanquantile(values, 0.98)
    scatter = ax.scatter(x, y, c=values, s=6.2, cmap=cmap, vmin=0, vmax=vmax,
                         alpha=0.82, linewidths=0)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_axis_off()
    cbar = ax.figure.colorbar(scatter, ax=ax, fraction=0.04, pad=0.015)
    cbar.ax.tick_params(labelsize=6.6, length=2)
    cbar.outline.set_linewidth(0.4)
    cbar.set_label("Normalized expression", fontsize=6.6)
    # Full-resolution pixels per nominal 55-um Visium capture spot come from
    # the original Space Ranger scalefactors, then scale to the plotted image.
    x0, y0 = image.shape[1] * 0.12, image.shape[0] * 0.91
    scale_line = ax.plot([x0, x0 + scalebar_px], [y0, y0], color="white", lw=2.0)[0]
    scale_line.set_path_effects([path_effects.Stroke(linewidth=3.1, foreground="black"), path_effects.Normal()])
    scale_text = ax.text(x0 + scalebar_px / 2, y0 - image.shape[0] * 0.025, "500 μm",
        color="white", ha="center", va="bottom", fontsize=6.6)
    scale_text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground="black"), path_effects.Normal()])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    set_style()
    project = args.project
    source = project / "results" / "spatial" / "GSE210616"
    data_dir = project / "source_metadata" / "spatial" / "GSE210616"
    out_dir = project / "figures" / "manuscript_nature_style"
    sections = pd.read_csv(source / "section_qc_and_IFI44L_spatial_autocorrelation.tsv", sep="\t")
    corr = pd.read_csv(source / "patient_IFI44L_signature_correlations.tsv", sep="\t")
    partial = pd.read_csv(source / "patient_partial_spearman.tsv", sep="\t")

    median_detection = sections["IFI44L_detection_fraction"].median()
    representative = sections.iloc[(sections["IFI44L_detection_fraction"] - median_detection).abs().argsort().iloc[0]]
    gsm = str(representative["gsm"])
    sample_dir = data_dir / "files" / gsm
    module = load_analysis_module(project)
    frame, _ = module.load_section(sample_dir)
    image_path = module.find_one(sample_dir, "*_tissue_hires_image.png.gz")
    scale_path = module.find_one(sample_dir, "*_scalefactors_json.json.gz")
    with gzip.open(image_path, "rb") as handle:
        image = np.asarray(Image.open(handle).convert("RGB"))
    with gzip.open(scale_path, "rt") as handle:
        scalefactors = json.load(handle)
    scale = float(scalefactors["tissue_hires_scalef"])
    spot_diameter_fullres_px = float(scalefactors["spot_diameter_fullres"])
    nominal_spot_diameter_um = 55.0
    scalebar_um = 500.0
    scalebar_px = scalebar_um / nominal_spot_diameter_um * spot_diameter_fullres_px * scale

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.02, 0.98], width_ratios=[1.04, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    ax_f = fig.add_subplot(gs[1, 2])

    spatial_panel(ax_a, image, frame, scale, "IFI44L", "IFI44L", "magma", scalebar_px)
    spatial_panel(ax_b, image, frame, scale, "IFN_response_excluding_IFI44L", "IFN response", "viridis", scalebar_px)
    spatial_panel(ax_c, image, frame, scale, "Epithelial", "Epithelial module", "cividis", scalebar_px)
    ax_a.text(0.01, -0.06, f"Prespecified section: {representative['patient_id']} / {gsm}",
              transform=ax_a.transAxes, fontsize=6.6, color=GREY, ha="left", va="top")
    panel_label(ax_a, "a", x=-0.08)
    panel_label(ax_b, "b", x=-0.08)
    panel_label(ax_c, "c", x=-0.08)

    patient_order = sections.groupby("patient_id")["IFI44L_moran_I"].median().sort_values().index.tolist()
    for i, patient in enumerate(patient_order):
        values = sections.loc[sections["patient_id"] == patient].sort_values("section")
        ax_d.plot(values["IFI44L_moran_I"], [i] * len(values), color=LIGHT_GREY, lw=1)
        for _, row in values.iterrows():
            marker = "o" if int(row["section"]) == 1 else "s"
            colour = BLUE if row["moran_fdr_bh"] < 0.05 and row["IFI44L_moran_I"] > 0 else GREY
            ax_d.scatter(row["IFI44L_moran_I"], i, marker=marker, s=20, color=colour,
                         edgecolor="white", linewidth=0.4, zorder=2)
    ax_d.axvline(0, color=GREY, lw=0.7, ls="--")
    ax_d.set_yticks(range(len(patient_order)), patient_order)
    ax_d.set_xlabel("IFI44L Moran's I")
    ax_d.set_title("Spatial autocorrelation (28 sections)", loc="left", fontweight="bold")
    panel_label(ax_d, "d", x=-0.22)

    order = ["IFN_response_excluding_IFI44L", "Myeloid", "T_cell", "Cytotoxicity", "Checkpoint", "Endothelial", "Epithelial"]
    labels = ["IFN response", "Myeloid", "T cell", "Cytotoxicity", "Checkpoint", "Endothelial", "Epithelial"]
    for i, sig in enumerate(order):
        values = corr.loc[corr["signature"] == sig, "spearman_rho"].to_numpy(float)
        values = values[np.isfinite(values)]
        colour = RED if sig == "Epithelial" else (TEAL if sig == "IFN_response_excluding_IFI44L" else BLUE)
        offsets = np.linspace(-0.10, 0.10, len(values))
        ax_e.scatter(values, i + offsets, color=colour, s=13,
                     alpha=0.72, edgecolor="white", linewidth=0.3)
        ax_e.plot([np.median(values)] * 2, [i - 0.18, i + 0.18], color=BLACK, lw=1.25)
    ax_e.axvline(0, color=GREY, lw=0.7, ls="--")
    ax_e.set_yticks(range(len(order)), labels)
    ax_e.set_xlabel("Within-tumour Spearman ρ")
    ax_e.set_title("Patient-level co-localization", loc="left", fontweight="bold")
    panel_label(ax_e, "e", x=-0.32)

    specs = [
        ("IFN_adjusted_for_lineage_modules", "IFN adjusted\nfor lineage", TEAL),
        ("Checkpoint_adjusted_for_IFN_and_lineage_modules", "Checkpoint adjusted\nfor IFN and lineage", GREY),
    ]
    for i, (analysis, label, colour) in enumerate(specs):
        values = partial.loc[partial["analysis"] == analysis, "partial_spearman_rho"].to_numpy(float)
        values = values[np.isfinite(values)]
        ax_f.scatter(i + np.linspace(-0.06, 0.06, len(values)), values, color=colour, s=18,
                     alpha=0.78, edgecolor="white", linewidth=0.35)
        ax_f.plot([i - 0.16, i + 0.16], [np.median(values)] * 2, color=BLACK, lw=1.35)
    ax_f.axhline(0, color=GREY, lw=0.7, ls="--")
    ax_f.set_xticks([0, 1], [x[1] for x in specs])
    ax_f.set_ylabel("Partial Spearman ρ")
    ax_f.set_title("Composition-adjusted associations", loc="left", fontweight="bold")
    panel_label(ax_f, "f", x=-0.17)

    base = out_dir / "Figure4_spatial_IFI44L_ecology"
    save_bundle(fig, base)
    audit = {
        "figure": "Figure 4 final",
        "base": str(base),
        "n_patients": 14,
        "n_sections": 28,
        "n_spots_after_qc": int(sections["n_spots_after_qc"].sum()),
        "representative_section": gsm,
        "representative_patient": str(representative["patient_id"]),
        "selection_rule": "section closest to cohort median IFI44L detection fraction",
        "image_processing": "unaltered H&E background with normalized-expression spot overlay",
        "inference_unit": "patient after averaging two sections",
        "scale_bar": {
            "length_um": scalebar_um,
            "nominal_Visium_spot_diameter_um": nominal_spot_diameter_um,
            "spot_diameter_fullres_pixels": spot_diameter_fullres_px,
            "tissue_hires_scalef": scale,
            "rendered_bar_pixels": scalebar_px,
            "formula": "length_um / 55 * spot_diameter_fullres * tissue_hires_scalef",
            "source_scalefactors": scale_path.name,
        },
    }
    (out_dir / "figure_build_audit_figure4_final.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
