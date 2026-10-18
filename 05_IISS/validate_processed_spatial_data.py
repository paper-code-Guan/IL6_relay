#!/usr/bin/env python3
"""Validate processed IISS spatial files against the released count matrices.

The CellGene tables are streamed one row at a time. Only one sample and its
sparse cell-by-gene counts are held in memory.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


SAMPLES = {
    "con2": "con2.txt",
    "con10": "con10.txt",
    "KO3": "ko3.txt",
    "KO11": "ko11.txt",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("source_data/IISS"),
        help="Directory containing the released IISS matrices.",
    )
    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=None,
        help="Directory containing per-sample processed spatial files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/IISS/processed_spatial_validation.csv"),
        help="CSV validation summary.",
    )
    return parser.parse_args()


def normalize_identifier(value: str) -> str:
    return value.strip()


def normalize_gene(value: str) -> str:
    gene = value.strip()
    if gene.startswith("IL6"):
        return "Il6" + gene[3:]
    return gene


def validate_tiff_header(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty image file: {path}")
    with path.open("rb") as handle:
        signature = handle.read(4)
    if signature not in {b"II*\x00", b"MM\x00*"}:
        raise ValueError(f"Invalid TIFF header: {path}")


def read_released_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Released matrix not found: {path}")
    matrix: dict[str, dict[str, int]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell.name", "celltype"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain cell.name and celltype columns")
        genes = [field for field in reader.fieldnames if field not in required]
        if not genes:
            raise ValueError(f"No gene columns found in {path}")
        for row in reader:
            cell_id = normalize_identifier(row["cell.name"])
            if not cell_id:
                raise ValueError(f"Blank cell.name in {path}")
            if cell_id in matrix:
                raise ValueError(f"Duplicated cell.name {cell_id!r} in {path}")
            counts: dict[str, int] = {}
            for gene in genes:
                value = float(row[gene])
                if value < 0 or not value.is_integer():
                    raise ValueError(
                        f"Non-integer count for {cell_id}, {gene} in {path}: {value}"
                    )
                counts[gene] = int(value)
            matrix[cell_id] = counts
    return genes, matrix


def stream_cellgene(
    path: Path,
) -> tuple[dict[str, Counter[str]], set[str], int, int]:
    if not path.is_file():
        raise FileNotFoundError(f"CellGene table not found: {path}")
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    genes: set[str] = set()
    rows = 0
    invalid_coordinates = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"gene", "y", "x", "cell"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain gene, y, x and cell columns")
        for row in reader:
            rows += 1
            gene = normalize_gene(row["gene"])
            cell_id = normalize_identifier(row["cell"])
            if not gene or not cell_id:
                raise ValueError(f"Blank gene or cell identifier in {path}, row {rows + 1}")
            genes.add(gene)
            counts[cell_id][gene] += 1
            try:
                float(row["x"])
                float(row["y"])
            except ValueError:
                invalid_coordinates += 1
    return counts, genes, rows, invalid_coordinates


def read_centroid_ids(path: Path) -> tuple[set[str], int]:
    if not path.is_file():
        raise FileNotFoundError(f"CellCentroid table not found: {path}")
    identifiers: set[str] = set()
    invalid_coordinates = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"row", "column", "cell"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain row, column and cell columns")
        for record in reader:
            cell_id = normalize_identifier(record["cell"])
            if not cell_id or cell_id in identifiers:
                raise ValueError(f"Blank or duplicated cell identifier in {path}: {cell_id!r}")
            identifiers.add(cell_id)
            try:
                float(record["row"])
                float(record["column"])
            except ValueError:
                invalid_coordinates += 1
    return identifiers, invalid_coordinates


def read_geojson_ids(path: Path) -> set[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Cell-boundary GeoJSON not found: {path}")
    with path.open("r", encoding="utf-8-sig") as handle:
        data = json.load(handle)
    if data.get("type") != "FeatureCollection":
        raise ValueError(f"Expected a FeatureCollection in {path}")
    identifiers: set[str] = set()
    for feature in data.get("features", []):
        cell_id = normalize_identifier(str(feature.get("properties", {}).get("name", "")))
        if not cell_id or cell_id in identifiers:
            raise ValueError(f"Blank or duplicated feature name in {path}: {cell_id!r}")
        identifiers.add(cell_id)
    return identifiers


def validate_sample(input_dir: Path, processed_dir: Path, sample: str) -> dict[str, object]:
    sample_dir = processed_dir / sample
    cells_dir = sample_dir / "cells"
    dapi_path = sample_dir / "DAPI.tif"
    mask_path = cells_dir / "masks.tif"
    for image_path in (dapi_path, mask_path):
        validate_tiff_header(image_path)

    genes, released = read_released_matrix(input_dir / SAMPLES[sample])
    decoded, decoded_genes, decoded_rows, invalid_signal_coordinates = stream_cellgene(
        cells_dir / "CellGene.csv"
    )
    centroid_ids, invalid_centroid_coordinates = read_centroid_ids(
        cells_dir / "CellCentroid.csv"
    )
    boundary_ids = read_geojson_ids(cells_dir / "Cells.geojson")

    if centroid_ids != boundary_ids:
        missing_boundaries = len(centroid_ids - boundary_ids)
        missing_centroids = len(boundary_ids - centroid_ids)
    else:
        missing_boundaries = 0
        missing_centroids = 0

    missing_genes = set(genes) - decoded_genes
    mismatched_cells = 0
    mismatched_values = 0
    for cell_id, expected in released.items():
        observed = decoded.get(cell_id, Counter())
        cell_has_mismatch = False
        for gene in genes:
            if observed.get(gene, 0) != expected[gene]:
                mismatched_values += 1
                cell_has_mismatch = True
        if cell_has_mismatch:
            mismatched_cells += 1

    released_ids = set(released)
    result = {
        "sample": sample,
        "released_matrix": SAMPLES[sample],
        "released_cells": len(released),
        "released_genes": len(genes),
        "decoded_rows": decoded_rows,
        "decoded_cells": len(decoded),
        "centroid_cells": len(centroid_ids),
        "boundary_cells": len(boundary_ids),
        "released_cells_missing_from_decoded": len(released_ids - set(decoded)),
        "released_cells_missing_from_centroids": len(released_ids - centroid_ids),
        "released_cells_missing_from_boundaries": len(released_ids - boundary_ids),
        "processed_cells_not_in_released_matrix": len(set(decoded) - released_ids),
        "decoded_cells_missing_from_centroids": len(set(decoded) - centroid_ids),
        "missing_released_genes_in_cellgene": len(missing_genes),
        "mismatched_cells": mismatched_cells,
        "mismatched_gene_values": mismatched_values,
        "centroids_without_boundaries": missing_boundaries,
        "boundaries_without_centroids": missing_centroids,
        "invalid_signal_coordinates": invalid_signal_coordinates,
        "invalid_centroid_coordinates": invalid_centroid_coordinates,
        "dapi_bytes": dapi_path.stat().st_size,
        "auxiliary_mask_bytes": mask_path.stat().st_size,
    }
    result["status"] = "PASS" if all(
        result[field] == 0
        for field in (
            "released_cells_missing_from_decoded",
            "released_cells_missing_from_centroids",
            "released_cells_missing_from_boundaries",
            "missing_released_genes_in_cellgene",
            "decoded_cells_missing_from_centroids",
            "mismatched_cells",
            "mismatched_gene_values",
            "centroids_without_boundaries",
            "boundaries_without_centroids",
            "invalid_signal_coordinates",
            "invalid_centroid_coordinates",
        )
    ) else "FAIL"
    return result


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    processed_dir = (
        args.processed_dir.resolve()
        if args.processed_dir is not None
        else input_dir / "processed_spatial"
    )
    results = [validate_sample(input_dir, processed_dir, sample) for sample in SAMPLES]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    failed = [result["sample"] for result in results if result["status"] != "PASS"]
    if failed:
        raise SystemExit(f"Validation failed for: {', '.join(failed)}")
    print(f"Validated {len(results)} samples; report written to {args.output}")


if __name__ == "__main__":
    main()
