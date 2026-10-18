#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from scipy.stats import pearsonr, shapiro, spearmanr


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Analyze the cohort VII IL6 and let-7b-5p correlation."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=(
            script_dir.parent.parent
            / "source_data"
            / "cohort_VII_IL6_let7b_expression_data.xlsx"
        ),
        help="Input workbook containing the Expression data worksheet.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "output",
        help="Directory for the PDF and statistics table.",
    )
    return parser.parse_args()


def font_properties(size: int, weight: str) -> fm.FontProperties:
    return fm.FontProperties(family="sans-serif", size=size, weight=weight)


def add_trendline(
    ax: plt.Axes, x_values: np.ndarray, y_values: np.ndarray
) -> None:
    valid = np.isfinite(x_values) & np.isfinite(y_values)
    fit_x = np.asarray(x_values[valid], dtype=float)
    fit_y = np.asarray(y_values[valid], dtype=float)
    if len(fit_x) < 3 or np.unique(fit_x).size < 2:
        return

    x_grid = np.linspace(fit_x.min(), fit_x.max(), 200)
    model = sm.OLS(fit_y, sm.add_constant(fit_x)).fit()
    prediction = model.get_prediction(sm.add_constant(x_grid))
    summary = prediction.summary_frame(alpha=0.05)
    line_color = "#555555"
    dash_style = (0, (1.2, 2.4))

    ax.plot(
        x_grid,
        summary["mean"].to_numpy(dtype=float),
        color=line_color,
        linewidth=2.2,
        alpha=0.95,
        zorder=1,
    )
    for column in ("mean_ci_lower", "mean_ci_upper"):
        ax.plot(
            x_grid,
            summary[column].to_numpy(dtype=float),
            color=line_color,
            linewidth=1.8,
            linestyle=dash_style,
            alpha=0.9,
            zorder=1,
        )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(args.input, sheet_name="Expression data").rename(
        columns={
            "IL6 log2(CPM + 1)": "IL6 log2(CPM+1)",
            "let-7b-5p log2(CPM + 1)": "let-7b-5p log2(CPM+1)",
        }
    )
    required = {
        "Group",
        "IL6 log2(CPM+1)",
        "let-7b-5p log2(CPM+1)",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df["Group"] = (
        df["Group"]
        .astype(str)
        .str.upper()
        .replace(
            {
                "SCHIZOPHRENIA": "SCZ",
                "SZ": "SCZ",
                "CASE": "SCZ",
                "HC": "CTR",
                "CONTROL": "CTR",
            }
        )
    )
    unexpected_groups = sorted(set(df["Group"]).difference({"SCZ", "CTR"}))
    if unexpected_groups:
        raise ValueError(f"Unexpected group labels: {unexpected_groups}")

    x = df["let-7b-5p log2(CPM+1)"].to_numpy(dtype=float)
    y = df["IL6 log2(CPM+1)"].to_numpy(dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Expression columns contain missing or non-finite values.")

    shapiro_il6 = shapiro(y)
    shapiro_let7b = shapiro(x)
    if shapiro_il6.pvalue > 0.05 and shapiro_let7b.pvalue > 0.05:
        method = "Pearson"
        coefficient_symbol = "r"
        coefficient, p_value = pearsonr(x, y)
    else:
        method = "Spearman"
        coefficient_symbol = r"$\rho$"
        coefficient, p_value = spearmanr(x, y)

    p_text = r"$P$ < 0.001" if p_value < 0.001 else rf"$P$ = {p_value:.3f}"
    title = (
        f"{method} {coefficient_symbol} = {coefficient:.3f}; {p_text}"
    )

    rng = np.random.default_rng(42)
    jitter_x = (x.max() - x.min()) * 0.004 if len(np.unique(x)) < len(x) else 0
    jitter_y = (y.max() - y.min()) * 0.004 if len(np.unique(y)) < len(y) else 0
    df["_x"] = x + rng.normal(0, jitter_x, size=len(df))
    df["_y"] = y + rng.normal(0, jitter_y, size=len(df))

    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
    ]
    matplotlib.rcParams["mathtext.fontset"] = "custom"
    matplotlib.rcParams["mathtext.rm"] = "sans"
    matplotlib.rcParams["mathtext.it"] = "sans:italic"
    matplotlib.rcParams["mathtext.bf"] = "sans:bold"
    label_font = font_properties(20, "bold")
    tick_font = font_properties(20, "normal")
    title_font = font_properties(20, "bold")
    legend_title_font = font_properties(18, "bold")
    legend_label_font = font_properties(18, "bold")

    fig, ax = plt.subplots(figsize=(5, 4.5))
    add_trendline(ax, x, y)
    sns.scatterplot(
        data=df,
        x="_x",
        y="_y",
        hue="Group",
        hue_order=["SCZ", "CTR"],
        s=90,
        alpha=0.9,
        edgecolor="white",
        linewidth=0.6,
        palette={"SCZ": "#F47F72", "CTR": "#7FB2D5"},
        zorder=3,
        ax=ax,
    )

    ax.set_xlabel(r"Let-7b-5p log$_2$(CPM+1)", fontproperties=label_font)
    ax.set_ylabel(r"$\mathit{IL6}$ log$_2$(CPM+1)", fontproperties=label_font)
    fig.suptitle(title, fontproperties=title_font, x=0.5, y=0.955)

    legend = ax.legend(
        title="Group",
        frameon=True,
        loc="upper right",
        bbox_to_anchor=(1, 1),
    )
    if legend is not None:
        legend.set_title("Group", prop=legend_title_font)
        for text in legend.get_texts():
            text.set_fontproperties(legend_label_font)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(2.3)

    ax.tick_params(axis="both", which="major", labelsize=20, width=2.3)
    for label in (*ax.get_xticklabels(), *ax.get_yticklabels()):
        label.set_fontproperties(tick_font)

    x_pad = (x.max() - x.min()) * 0.08 if x.max() > x.min() else 0.5
    y_pad = (y.max() - y.min()) * 0.08 if y.max() > y.min() else 0.5
    ax.set_xlim(x.min() - x_pad, x.max() + x_pad)
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    output_pdf = args.output_dir / "Fig2f_cohort_VII_IL6_let7b_correlation.pdf"
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        [
            {
                "n": len(df),
                "method": method,
                "coefficient": coefficient,
                "p_value": p_value,
                "shapiro_IL6_statistic": shapiro_il6.statistic,
                "shapiro_IL6_p_value": shapiro_il6.pvalue,
                "shapiro_let7b_statistic": shapiro_let7b.statistic,
                "shapiro_let7b_p_value": shapiro_let7b.pvalue,
            }
        ]
    ).to_csv(
        args.output_dir / "Fig2f_cohort_VII_IL6_let7b_statistics.csv",
        index=False,
    )

    print(
        f"{method} correlation: coefficient={coefficient:.12g}, "
        f"p={p_value:.12g}, n={len(df)}"
    )
    print(f"Figure: {output_pdf}")


if __name__ == "__main__":
    main()
