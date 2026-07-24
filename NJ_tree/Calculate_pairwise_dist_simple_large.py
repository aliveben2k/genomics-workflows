#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd


def calculate_pairwise_distances(input_file):
    input_file = Path(input_file)
    print(f"Reading input file: {input_file}")

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")
    if "tmp.txt" not in input_file.name:
        raise ValueError(
            f"Input filename must contain 'tmp.txt': {input_file}"
        )
    if input_file.stat().st_size == 0:
        raise ValueError(f"Input file is empty: {input_file}")

    try:
        data = pd.read_csv(input_file, sep="\t", header=0)
    except pd.errors.EmptyDataError as error:
        raise ValueError(
            f"Input file contains no readable data: {input_file}"
        ) from error

    if data.empty:
        raise ValueError(f"Input file contains no data rows: {input_file}")
    if data.shape[1] < 7:
        raise ValueError(
            "Input table must contain five SNP-information columns "
            "and at least two sample columns"
        )

    sample_names = data.columns[5:].astype(str).tolist()
    region_start_pos = str(data.iloc[0, 0])
    data_use = (
        data.iloc[:, 5:]
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy(dtype=float)
    )

    print(f"Calculating pairwise distances for {input_file.name}")
    n_samples = data_use.shape[1]
    available_values = []
    difference_values = []

    for i in range(n_samples - 1):
        for j in range(i + 1, n_samples):
            absolute_difference = np.abs(data_use[:, i] - data_use[:, j])
            available_values.append(np.count_nonzero(~np.isnan(absolute_difference)))
            difference_values.append(np.nansum(absolute_difference))

    available_values = np.asarray(available_values, dtype=int)
    difference_values = np.asarray(difference_values, dtype=float)
    upper_rows, upper_columns = np.triu_indices(n_samples, k=1)

    available_matrix = np.zeros((n_samples, n_samples), dtype=int)
    available_matrix[upper_rows, upper_columns] = available_values
    available_matrix[upper_columns, upper_rows] = available_values

    difference_matrix = np.zeros((n_samples, n_samples), dtype=float)
    difference_matrix[upper_rows, upper_columns] = difference_values
    difference_matrix[upper_columns, upper_rows] = difference_values

    return (
        sample_names,
        available_matrix,
        difference_matrix,
        region_start_pos,
    )


def process_chunk(input_file, output_file):
    output_file = Path(output_file)
    (
        sample_names,
        available_matrix,
        difference_matrix,
        region_start_pos,
    ) = calculate_pairwise_distances(input_file)

    np.savez_compressed(
        output_file,
        sample_names=np.asarray(sample_names),
        available_matrix=available_matrix,
        difference_matrix=difference_matrix,
        region_start_pos=np.asarray(region_start_pos),
    )
    print(f"Saved: {output_file}")
    return output_file


def task_from_manifest(input_manifest, task_id):
    input_manifest = Path(input_manifest)
    if not input_manifest.is_file():
        raise FileNotFoundError(f"Manifest does not exist: {input_manifest}")

    with open(input_manifest, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"task_id", "chunk_file", "matrix_file"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Manifest must contain columns {sorted(required)}: "
                f"{input_manifest}"
            )
        for row in reader:
            if int(row["task_id"]) == task_id:
                return Path(row["chunk_file"]), Path(row["matrix_file"])

    raise ValueError(f"Task ID {task_id} is not present in {input_manifest}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate pairwise distances for one genotype-table chunk."
    )
    parser.add_argument("-i", "--input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "-m",
        "--manifest",
        type=Path,
        help="Manifest created by vcf2table_large.py prepare",
    )
    parser.add_argument(
        "-t",
        "--task-id",
        type=int,
        help="One-based task ID from the manifest",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    direct_mode = args.input is not None or args.output is not None
    manifest_mode = args.manifest is not None or args.task_id is not None

    if direct_mode and manifest_mode:
        raise ValueError(
            "Use either --input/--output or --manifest/--task-id, not both"
        )
    if direct_mode:
        if args.input is None or args.output is None:
            raise ValueError("Both --input and --output are required")
        input_file, output_file = args.input, args.output
    elif manifest_mode:
        if args.manifest is None or args.task_id is None:
            raise ValueError("Both --manifest and --task-id are required")
        if args.task_id <= 0:
            raise ValueError("--task-id must be greater than zero")
        input_file, output_file = task_from_manifest(
            args.manifest,
            args.task_id,
        )
    else:
        raise ValueError(
            "Provide --input/--output or --manifest/--task-id"
        )

    process_chunk(input_file, output_file)


if __name__ == "__main__":
    main()
