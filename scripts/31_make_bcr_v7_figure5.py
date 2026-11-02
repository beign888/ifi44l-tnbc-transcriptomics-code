#!/usr/bin/env python3
"""Create the patient-free conceptual Figure 5 for the BCR package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BLUE = "#2B6CA3"
TEAL = "#23877D"
RED = "#B64E55"
GREY = "#66727A"
BLACK = "#20262B"
BLUE_BG = "#F1F6FA"
TEAL_BG = "#EEF8F6"
RED_BG = "#FCF1F2"


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    lines: list[str],
    edge: str,
    face: str,
    title_size: float = 7.1,
    body_size: float = 6.0,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.014",
        linewidth=1.0,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height - 0.050,
        title,
        ha="center",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        color=edge,
    )
    ax.text(
        x + 0.022,
        y + height - 0.120,
        "\n".join(lines),
        ha="left",
        va="top",
        fontsize=body_size,
        linespacing=1.34,
        color=BLACK,
    )


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], colour: str, dashed: bool = False) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=1.0,
            linestyle="--" if dashed else "-",
            color=colour,
            connectionstyle="arc3,rad=0",
        )
    )


def build(
    out_dir: Path,
    heading: str = "Public multiscale transcriptomic TNBC evidence",
    version: str = "v7_no_local_human_cohort",
) -> Path:
    set_style()
    fig, ax = plt.subplots(figsize=(7.15, 4.45), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.045, 0.965, heading, fontsize=8.2, fontweight="bold", color=TEAL, va="top")
    ax.text(0.785, 0.965, "Auxiliary evidence", fontsize=8.2, fontweight="bold", color=RED, va="top")

    box(
        ax,
        0.035,
        0.58,
        0.215,
        0.275,
        "Bulk cohorts",
        [
            "IFI44L tracks IFN, antigen-",
            "presentation and T-cell",
            "programmes; favourable",
            "outcome association in",
            "GSE58812",
        ],
        BLUE,
        BLUE_BG,
    )
    box(
        ax,
        0.285,
        0.58,
        0.215,
        0.275,
        "Single-cell source",
        [
            "Highest in myeloid,",
            "endothelial and T-cell",
            "compartments; lower in",
            "malignant epithelium",
        ],
        TEAL,
        TEAL_BG,
    )
    box(
        ax,
        0.535,
        0.58,
        0.215,
        0.275,
        "Spatial organization",
        [
            "Recurrent IFI44L clustering",
            "within IFN- and immune-",
            "enriched multicellular",
            "niches",
        ],
        TEAL,
        TEAL_BG,
    )
    box(
        ax,
        0.785,
        0.58,
        0.18,
        0.275,
        "Cell perturbation",
        [
            "Two TNBC lines",
            "Qualitative immunoblot",
            "Reduced growth,",
            "migration and invasion",
            "after knockdown",
        ],
        RED,
        RED_BG,
        body_size=5.7,
    )

    box(
        ax,
        0.205,
        0.19,
        0.50,
        0.235,
        "Compartment-aware interpretation",
        [
            "At tissue level, IFI44L primarily marks an organized IFN-inflamed",
            "microenvironment rather than a tumour-cell-only programme.",
            "Cell-intrinsic effects may coexist at a distinct biological scale.",
        ],
        BLUE,
        BLUE_BG,
        title_size=7.5,
        body_size=6.3,
    )
    box(
        ax,
        0.755,
        0.19,
        0.21,
        0.235,
        "Evidence boundary",
        [
            "n = 3 independent",
            "experiments per assay",
            "No rescue or downstream",
            "mechanism tested",
        ],
        RED,
        RED_BG,
        body_size=5.8,
    )

    arrow(ax, (0.142, 0.58), (0.31, 0.425), BLUE)
    arrow(ax, (0.392, 0.58), (0.43, 0.425), TEAL)
    arrow(ax, (0.642, 0.58), (0.58, 0.425), TEAL)
    arrow(ax, (0.875, 0.58), (0.86, 0.425), RED)
    arrow(ax, (0.755, 0.305), (0.705, 0.305), RED, dashed=True)

    ax.text(
        0.49,
        0.085,
        "Scale-dependent synthesis: public tumour data define the tissue-level signal; cell-line perturbation remains supportive.",
        ha="center",
        va="center",
        fontsize=6.3,
        fontweight="bold",
        color=BLACK,
    )
    ax.text(
        0.49,
        0.040,
        "Arrows summarize the evidence logic and do not imply a causal pathway.",
        ha="center",
        va="center",
        fontsize=5.6,
        color=GREY,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_dir / "Figure5_compartment_dependent_model_no_local_human_cohort"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    contract = {
        "figure": "Figure 5",
        "version": version,
        "core_conclusion": "Tissue-level IFI44L primarily marks IFN-inflamed multicellular niches, while cell-line perturbation supports a distinct auxiliary cell-intrinsic phenotype.",
        "question": "How can tissue-level immune enrichment and an isolated tumour-cell phenotype coexist?",
        "archetype": "schematic-led synthesis",
        "backend": "Python/matplotlib",
        "size": "double-column, approximately 170 mm",
        "causal_claim": False,
        "human_biospecimen_component": False,
        "evidence_boundaries": [
            "WB is qualitative",
            "three independent experiments per assay",
            "no rescue experiment",
            "no downstream mediator tested",
        ],
        "files": [str(base.with_suffix(ext)) for ext in [".pdf", ".svg", ".png", ".tiff"]],
    }
    contract_path = out_dir / "Figure5_contract_v7.json"
    contract_path.write_text(json.dumps(contract, indent=2), encoding="utf-8")
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--heading", default="Public multiscale transcriptomic TNBC evidence")
    parser.add_argument("--version", default="v7_no_local_human_cohort")
    args = parser.parse_args()
    if args.output_dir is not None:
        output_dir = args.output_dir
    elif args.project is not None:
        output_dir = args.project / "figures/v7_no_local_human_cohort"
    else:
        parser.error("one of --project or --output-dir is required")
    print(build(output_dir, heading=args.heading, version=args.version))
