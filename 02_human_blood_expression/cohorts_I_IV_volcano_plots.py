#!/usr/bin/env python3
"""Reproduce the four discovery-cohort miRNA volcano plots."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = (
    SCRIPT_DIR.parent.parent
    / "source_data"
    / "cohorts_I_IV_miRNA_volcano_plot_data.xlsx"
)

COHORT_CONFIGS = (
    {
        "sheet": "Cohort I",
        "title": "Cohort 1",
        "output": "cohort_I_volcano.pdf",
        "xticks_raw": [-4, -2, 0, 2, 4, 6],
        "ymax_tick": 10,
        "yticks_step": 1,
        "ylabel_y": 7.5,
        "internal_title_x": 0.77,
        "legend_x": 0.80,
    },
    {
        "sheet": "Cohort II",
        "title": "Cohort 2",
        "output": "cohort_II_volcano.pdf",
        "xticks_raw": [-4, -2, -1, 0, 1, 2, 4],
        "ymax_tick": 8,
        "yticks_step": 1,
        "ylabel_y": 6.0,
        "internal_title_x": 0.74,
        "legend_x": 0.75,
    },
    {
        "sheet": "Cohort III",
        "title": "Cohort 3",
        "output": "cohort_III_volcano.pdf",
        "xticks_raw": [-4, -2, -1, 0, 1, 2, 4],
        "ymax_tick": 15,
        "yticks_step": 2,
        "ylabel_y": 11.0,
        "internal_title_x": 0.74,
        "legend_x": 0.75,
    },
    {
        "sheet": "Cohort IV",
        "title": "Cohort 4",
        "output": "cohort_IV_volcano.pdf",
        "xticks_raw": [-4, -2, -1, 0, 1, 2, 4],
        "ymax_tick": 8,
        "yticks_step": 1,
        "ylabel_y": 6.0,
        "internal_title_x": 0.74,
        "legend_x": 0.75,
    },
)

HIGHLIGHT_GENES = ("hsa-let-7b-5p", "hsa-miR-93-5p", "hsa-miR-4508")
HIGHLIGHT_COLORS = ("#F47F72", "#FFB347", "#A37ACC")
UP_COLOR = "#8DD2C5"
DOWN_COLOR = "#7FB2D5"
NONSIGNIFICANT_COLOR = "#EBEBEB"
P_THRESHOLD = -np.log10(0.05)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Workbook containing sheets Cohort I through Cohort IV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "volcano_outputs",
        help="Directory for the four PDF files.",
    )
    parser.add_argument(
        "--font",
        type=Path,
        default=None,
        help="Optional Arial-compatible TTF file.",
    )
    return parser.parse_args()


def resolve_font(font_path: Path | None) -> font_manager.FontProperties:
    if font_path is not None:
        if not font_path.is_file():
            raise FileNotFoundError(f"Font file not found: {font_path}")
        return font_manager.FontProperties(fname=str(font_path))

    for family in ("Arial", "Helvetica", "DejaVu Sans"):
        try:
            resolved = font_manager.findfont(family, fallback_to_default=False)
            return font_manager.FontProperties(fname=resolved)
        except ValueError:
            continue
    return font_manager.FontProperties()


def stretch_fc(values: pd.Series | np.ndarray, scale: float = 3.0):
    values = np.asarray(values, dtype=float)
    return np.sign(values) * np.log2(np.abs(values) + 1) * scale


def load_cohort(workbook: Path, sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(workbook, sheet_name=sheet)
    required = {"miRNA", "log2 fold change", "-log10(P value)"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{sheet} is missing columns: {sorted(missing)}")

    frame = frame.loc[
        :,
        ["miRNA", "log2 fold change", "-log10(P value)"],
    ].rename(
        columns={
            "miRNA": "gene",
            "log2 fold change": "log2FC",
            "-log10(P value)": "-log10(P)",
        }
    )
    frame["log2FC"] = pd.to_numeric(frame["log2FC"], errors="raise")
    frame["-log10(P)"] = pd.to_numeric(frame["-log10(P)"], errors="raise")
    if frame.isna().any().any():
        raise ValueError(f"{sheet} contains missing plotting values.")
    if frame["gene"].duplicated().any():
        duplicates = frame.loc[frame["gene"].duplicated(), "gene"].tolist()
        raise ValueError(f"{sheet} contains duplicated miRNAs: {duplicates}")
    return frame


def plot_cohort(
    frame: pd.DataFrame,
    config: dict,
    output_path: Path,
    font: font_manager.FontProperties,
) -> None:
    xtick_fontsize = 28
    ytick_fontsize = 28
    label_fontsize = 28
    annot_fontsize = 28
    point_size = 60
    highlight_point_size = 200
    alpha = 0.95
    axis_linewidth = 2.5
    hline_linewidth = 2.5
    p05_text_dx = 2.5
    p05_text_dy = 0.10
    ylabel_x = 0.3
    ylabel_rotation = 90
    ymin_fixed = 0.8
    xlim_margin = 0.5
    internal_title_y_ratio = 0.75
    legend_point_size = 100
    legend_up_y_ratio = 0.68
    legend_down_y_ratio = 0.60
    legend_fontsize = 28

    frame = frame.copy()
    frame["log2FC_stretch"] = stretch_fc(frame["log2FC"])
    frame["color"] = NONSIGNIFICANT_COLOR
    frame.loc[
        (frame["-log10(P)"] >= P_THRESHOLD) & (frame["log2FC"] > 0),
        "color",
    ] = UP_COLOR
    frame.loc[
        (frame["-log10(P)"] >= P_THRESHOLD) & (frame["log2FC"] < 0),
        "color",
    ] = DOWN_COLOR
    for gene, color in zip(HIGHLIGHT_GENES, HIGHLIGHT_COLORS):
        frame.loc[frame["gene"] == gene, "color"] = color

    fig, ax = plt.subplots(figsize=(7, 8))
    ax.scatter(
        frame["log2FC_stretch"],
        frame["-log10(P)"],
        c=frame["color"],
        s=point_size,
        linewidths=0,
        alpha=alpha,
        zorder=1,
    )
    for gene, color in zip(HIGHLIGHT_GENES, HIGHLIGHT_COLORS):
        selected = frame.loc[frame["gene"] == gene]
        if not selected.empty:
            row = selected.iloc[0]
            ax.scatter(
                [row["log2FC_stretch"]],
                [row["-log10(P)"]],
                c=color,
                s=highlight_point_size,
                linewidths=0,
                alpha=alpha,
                zorder=2,
            )

    ax.axhline(
        y=P_THRESHOLD,
        color="black",
        linestyle="dashed",
        linewidth=hline_linewidth,
        zorder=4,
    )
    ax.text(
        frame["log2FC_stretch"].max() + p05_text_dx,
        P_THRESHOLD + p05_text_dy,
        r"$P = 0.05$",
        fontproperties=font,
        fontsize=annot_fontsize,
        va="bottom",
        ha="right",
        style="italic",
        fontweight="bold",
        zorder=5,
    )

    ax.spines["left"].set_position(("data", 0))
    ax.spines["left"].set_linewidth(axis_linewidth)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_linewidth(axis_linewidth)

    xticks_raw = np.asarray(config["xticks_raw"])
    xticks_stretch = stretch_fc(xticks_raw)
    ax.set_xticks(xticks_stretch)
    ax.set_xticklabels(
        [str(value) for value in xticks_raw],
        fontproperties=font,
        fontsize=xtick_fontsize,
        fontweight="bold",
    )
    right_lim = max(
        frame["log2FC_stretch"].max() + xlim_margin,
        xticks_stretch.max() + xlim_margin,
    )
    left_lim = frame["log2FC_stretch"].min() - xlim_margin
    ax.set_xlim(left_lim, right_lim)

    yticks = np.arange(
        1,
        config["ymax_tick"] + 1,
        config["yticks_step"],
    )
    ax.set_yticks(yticks)
    ax.set_yticklabels(
        [str(value) for value in yticks],
        fontproperties=font,
        fontsize=ytick_fontsize,
        fontweight="bold",
    )
    ax.set_ylim(ymin_fixed, config["ymax_tick"])
    ax.set_xlabel(
        r"Stretched log$_2$FC",
        fontproperties=font,
        fontsize=label_fontsize,
        fontweight="bold",
    )
    ax.text(
        ylabel_x,
        config["ylabel_y"],
        r"$-\log_{10}(P)$",
        fontproperties=font,
        fontsize=label_fontsize,
        va="bottom",
        ha="left",
        rotation=ylabel_rotation,
        fontweight="bold",
        zorder=10,
    )

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    title_x = x_min + (x_max - x_min) * config["internal_title_x"]
    title_y = y_min + (y_max - y_min) * internal_title_y_ratio
    ax.text(
        title_x,
        title_y,
        config["title"],
        fontproperties=font,
        fontsize=label_fontsize,
        fontweight="bold",
        ha="left",
        va="center",
        color="black",
        zorder=20,
    )

    legend_x = x_min + (x_max - x_min) * config["legend_x"]
    legend_up_y = y_min + (y_max - y_min) * legend_up_y_ratio
    legend_down_y = y_min + (y_max - y_min) * legend_down_y_ratio
    for y_value, color, label in (
        (legend_up_y, UP_COLOR, "Up"),
        (legend_down_y, DOWN_COLOR, "Down"),
    ):
        ax.scatter(
            [legend_x],
            [y_value],
            c=color,
            s=legend_point_size,
            linewidth=0,
            alpha=alpha,
            zorder=21,
        )
        ax.text(
            legend_x + 0.3,
            y_value,
            label,
            fontproperties=font,
            fontsize=legend_fontsize,
            fontweight="bold",
            ha="left",
            va="center",
            color="black",
            zorder=21,
        )

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input workbook not found: {input_path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams["mathtext.fontset"] = "custom"
    matplotlib.rcParams["mathtext.rm"] = "sans"
    matplotlib.rcParams["mathtext.it"] = "sans:italic"
    matplotlib.rcParams["mathtext.bf"] = "sans:bold"
    font = resolve_font(args.font)

    for config in COHORT_CONFIGS:
        cohort = load_cohort(input_path, config["sheet"])
        output_path = output_dir / config["output"]
        plot_cohort(cohort, config, output_path, font)
        print(f"{config['sheet']}: {len(cohort)} points -> {output_path}")


if __name__ == "__main__":
    main()
