#!/usr/bin/env python3

import argparse
import os
import random
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from matplotlib import font_manager
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


FILES = ("ko3.txt", "ko11.txt", "con2.txt", "con10.txt")
ID_COLUMN = "cell.name"
TYPE_COLUMN = "celltype"
ALLOWED_TYPES = (
    "ExciNeuN",
    "InhiNeuN",
    "Astrocyte",
    "Microglia",
    "Endo",
)
TYPE_LABELS = {
    "ExciNeuN": "ExtNeu",
    "InhiNeuN": "InbNeu",
    "Astrocyte": "Astro",
    "Microglia": "MiGlia",
    "Endo": "Endoth",
}
COLORS = {
    "ExciNeuN": "#377eb8",
    "InhiNeuN": "#a6cee3",
    "Astrocyte": "#e41a1c",
    "Microglia": "#4daf4a",
    "Endo": "#984ea3",
}
ZORDER_PRIORITY = ("Endo", "Microglia", "Astrocyte", "InhiNeuN", "ExciNeuN")
EXPECTED_CELL_COUNT = 17554
RANDOM_SEED = 2025

warnings.filterwarnings(
    "ignore",
    message="n_jobs value .* overridden .*random_state",
    category=UserWarning,
)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_input = Path(
        os.environ.get(
            "IISS_SOURCE_DATA_DIR",
            script_dir.parent.parent / "source_data" / "IISS",
        )
    )
    default_output = Path(
        os.environ.get("IISS_OUTPUT_DIR", script_dir / "output")
    )
    parser = argparse.ArgumentParser(
        description="Generate the five-cell-type IISS UMAP for Figure 3k."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input,
        help="Directory containing con2, con10, ko3 and ko11 tables.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help="Directory for the UMAP PDF and Source Data.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def configure_font() -> None:
    for font_family in ("Arial", "Helvetica", "DejaVu Sans"):
        try:
            font_manager.findfont(font_family, fallback_to_default=False)
            plt.rcParams["font.family"] = font_family
            return
        except ValueError:
            continue


def load_samples(input_dir: Path) -> pd.DataFrame:
    samples = []
    for file_name in FILES:
        file_path = input_dir / file_name
        if not file_path.is_file():
            raise FileNotFoundError(f"Missing IISS input: {file_path}")
        sample = pd.read_csv(file_path, sep=None, engine="python")
        sample.columns = [column.strip() for column in sample.columns]
        missing = {ID_COLUMN, TYPE_COLUMN}.difference(sample.columns)
        if missing:
            raise ValueError(
                f"{file_name} is missing columns: {sorted(missing)}"
            )
        sample_name = file_path.stem
        sample["sample"] = sample_name
        sample["group"] = (
            "Syn-control"
            if sample_name.startswith("con")
            else "Syn-let-7b KO"
        )
        samples.append(sample)
    return pd.concat(samples, axis=0, ignore_index=True)


def filter_cell_types(data: pd.DataFrame) -> pd.DataFrame:
    labels = data[TYPE_COLUMN].astype("string")
    astrocyte_subtypes = labels.str.endswith("_Astrocyte", na=False)
    data.loc[astrocyte_subtypes, TYPE_COLUMN] = "Astrocyte"
    data = data[data[TYPE_COLUMN].isin(ALLOWED_TYPES)].copy()
    if data.empty:
        raise ValueError("No cells remained after cell-type filtering.")
    data[TYPE_COLUMN] = pd.Categorical(
        data[TYPE_COLUMN],
        categories=ALLOWED_TYPES,
        ordered=True,
    )
    observed_types = set(data[TYPE_COLUMN].dropna().unique())
    if observed_types != set(ALLOWED_TYPES):
        raise ValueError(
            f"Expected exactly five cell types, found: {sorted(observed_types)}"
        )
    if len(data) != EXPECTED_CELL_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_CELL_COUNT} cells, found {len(data)}."
        )
    return data


def transform_features(data: pd.DataFrame) -> np.ndarray:
    metadata = {ID_COLUMN, TYPE_COLUMN, "sample", "group"}
    features = [column for column in data.columns if column not in metadata]
    matrix = (
        data[features]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
    )
    retained = (matrix > 0).sum(axis=0) >= 10
    matrix = matrix.loc[:, retained.values]
    if matrix.shape[1] == 0:
        raise ValueError("No features passed the nonzero-cell filter.")
    return np.log1p(matrix).to_numpy(dtype=float)


def class_balanced_weights(labels: np.ndarray) -> np.ndarray:
    values, counts = np.unique(labels, return_counts=True)
    counts_by_label = dict(zip(values, counts))
    weights = np.array(
        [1.0 / counts_by_label[label] for label in labels],
        dtype=float,
    )
    weights *= len(labels) / weights.sum()
    return weights


def class_balanced_indices(
    labels: np.ndarray,
    cap_per_class: int,
    rng: np.random.RandomState,
) -> np.ndarray:
    all_indices = np.arange(len(labels))
    selected = []
    for cell_type in pd.unique(labels):
        class_indices = all_indices[labels == cell_type]
        if len(class_indices) > cap_per_class:
            class_indices = rng.choice(
                class_indices,
                size=cap_per_class,
                replace=False,
            )
        selected.append(class_indices)
    return np.concatenate(selected)


def calculate_embedding(data: pd.DataFrame) -> np.ndarray:
    labels = data[TYPE_COLUMN].astype(str).to_numpy()
    class_codes = data[TYPE_COLUMN].cat.codes.to_numpy()
    matrix = transform_features(data)
    matrix = StandardScaler(with_mean=True, with_std=True).fit_transform(matrix)

    weights = class_balanced_weights(labels)
    centered = matrix - np.average(matrix, axis=0, weights=weights)
    n_components = min(50, centered.shape[1])
    pca_matrix = PCA(
        n_components=n_components,
        random_state=RANDOM_SEED,
    ).fit_transform(centered)

    reducer = umap.UMAP(
        n_neighbors=45,
        min_dist=0.06,
        metric="cosine",
        n_components=2,
        random_state=RANDOM_SEED,
        target_metric="categorical",
        target_weight=0.7,
        densmap=True,
        repulsion_strength=1.0,
        verbose=False,
    )
    rng = np.random.RandomState(RANDOM_SEED)
    balanced_indices = class_balanced_indices(labels, 1500, rng)
    remaining_indices = np.setdiff1d(
        np.arange(len(labels)),
        balanced_indices,
        assume_unique=False,
    )
    fit_order = np.concatenate([balanced_indices, remaining_indices])
    ordered_embedding = reducer.fit_transform(
        pca_matrix[fit_order],
        y=class_codes[fit_order],
    )
    inverse_order = np.empty_like(fit_order)
    inverse_order[fit_order] = np.arange(len(fit_order))
    embedding = ordered_embedding[inverse_order]

    center = embedding.mean(axis=0, keepdims=True)
    embedding = center + (embedding - center) * 0.2
    jitter_rng = np.random.RandomState(RANDOM_SEED)
    embedding = embedding + jitter_rng.normal(
        0,
        0.003,
        size=embedding.shape,
    )
    return embedding


def add_umap_axes(fig: plt.Figure, ax: plt.Axes) -> None:
    box = ax.get_position()
    inset_width = box.width * 0.28
    inset_height = inset_width
    x_position = max(0.02, box.x0 - inset_width - 0.06)
    y_position = max(0.02, box.y0 - inset_height - 0.06)
    inset = fig.add_axes(
        [x_position, y_position, inset_width, inset_height]
    )
    inset.set_axis_off()
    inset.set_xlim(0, 1)
    inset.set_ylim(0, 1)
    arrow_properties = {
        "arrowstyle": "-|>",
        "lw": 4.0,
        "color": "black",
        "shrinkA": 0,
        "shrinkB": 0,
    }
    inset.annotate(
        "",
        xy=(0.85, 0.15),
        xytext=(0.15, 0.15),
        arrowprops=arrow_properties,
    )
    inset.annotate(
        "",
        xy=(0.15, 0.85),
        xytext=(0.15, 0.15),
        arrowprops=arrow_properties,
    )
    inset.text(
        0.90,
        0.12,
        "UMAP 1",
        fontsize=30,
        ha="left",
        va="center",
    )
    inset.text(
        0.12,
        0.90,
        "UMAP 2",
        fontsize=30,
        ha="center",
        va="bottom",
        rotation=90,
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    set_seed(RANDOM_SEED)
    configure_font()
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["font.size"] = 25

    data = filter_cell_types(load_samples(args.input_dir))
    labels = data[TYPE_COLUMN].astype(str).to_numpy()
    embedding = calculate_embedding(data)

    source_data = pd.DataFrame(
        {
            "cell_id": data[ID_COLUMN].astype(str).to_numpy(),
            "sample": data["sample"].astype(str).to_numpy(),
            "group": data["group"].astype(str).to_numpy(),
            "cell_type": labels,
            "cell_type_label": [TYPE_LABELS[label] for label in labels],
            "UMAP_1": embedding[:, 0],
            "UMAP_2": embedding[:, 1],
        }
    )
    source_file = args.output_dir / "Fig3k_IISS_UMAP_source_data.csv"
    source_data.to_csv(source_file, index=False)

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=300)
    for z_order, cell_type in enumerate(ZORDER_PRIORITY, start=1):
        selected = labels == cell_type
        if not np.any(selected):
            continue
        ax.scatter(
            embedding[selected, 0],
            embedding[selected, 1],
            s=9,
            c=COLORS[cell_type],
            alpha=1,
            edgecolors="none",
            linewidths=0.0,
            label=TYPE_LABELS[cell_type],
            zorder=z_order,
            rasterized=True,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.tight_layout()
    add_umap_axes(fig, ax)

    output_pdf = args.output_dir / "Fig3k_IISS_UMAP.pdf"
    fig.savefig(output_pdf)
    plt.close(fig)
    print(f"Figure: {output_pdf}")
    print(f"Source Data: {source_file}")


if __name__ == "__main__":
    main()
