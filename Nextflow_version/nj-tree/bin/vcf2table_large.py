#!/usr/bin/env python3

import argparse
import csv
import gzip
import pickle
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from Bio import Phylo
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.tree import nj

from Calculate_pairwise_dist_simple_large import process_chunk


def positive_integer(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return number


def chunk_output_path(output_prefix, chunk_number):
    output_name = str(output_prefix)
    if output_name.endswith(".txt.gz"):
        prefix = output_name.removesuffix(".txt.gz")
    elif output_name.endswith(".gz"):
        prefix = output_name.removesuffix(".gz")
    elif output_name.endswith(".txt"):
        prefix = output_name.removesuffix(".txt")
    else:
        prefix = output_name
    return Path(f"{prefix}.{chunk_number}.tmp.txt.gz")


def matrix_output_path(chunk_file):
    name = chunk_file.name
    if not name.endswith(".tmp.txt.gz"):
        raise ValueError(f"Unexpected chunk filename: {chunk_file}")
    return chunk_file.with_name(
        name.removesuffix(".tmp.txt.gz") + ".tmp.npz"
    )


def open_output_chunk(output_prefix, chunk_number, output_names):
    chunk_file = chunk_output_path(output_prefix, chunk_number)
    output_handle = gzip.open(
        chunk_file,
        mode="wt",
        encoding="utf-8",
        newline="",
    )
    writer = csv.writer(
        output_handle,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writerow(["SNP", "chr", "pos", "ref", "alt"] + output_names)
    return chunk_file, output_handle, writer


def genotype_value(genotype):
    if genotype is None or len(genotype) != 2:
        return None
    allele1, allele2 = genotype
    if allele1 is None or allele2 is None:
        return None
    if allele1 not in (0, 1) or allele2 not in (0, 1):
        return None
    return (allele1 + allele2) / 2


def read_population_file(pop_file):
    population_ids = {}
    with open(pop_file, "rt", encoding="utf-8") as population_handle:
        lines = [
            line.rstrip("\r\n")
            for line in population_handle
            if line.strip()
        ]

    if lines and re.match(
        r"(?:ind\B|id|taxa|samp\B)",
        lines[0],
        flags=re.IGNORECASE,
    ):
        lines.pop(0)

    for line_number, line in enumerate(lines, start=1):
        elements = line.split("\t")
        if len(elements) < 2:
            raise ValueError(
                f"Population line {line_number} has fewer than two columns"
            )
        sample_id, population = elements[0], elements[1]
        population_ids.setdefault(population, []).append(sample_id)

    if not population_ids:
        raise ValueError(f"Population file contains no sample records: {pop_file}")
    return population_ids


def population_indices(population_ids, sample_names):
    sample_index = {
        sample_name: index
        for index, sample_name in enumerate(sample_names)
    }
    return {
        population: [
            sample_index[sample_id]
            for sample_id in population_ids[population]
            if sample_id in sample_index
        ]
        for population in sorted(population_ids)
    }


def population_frequencies(genotype_values, pop_freq, populations):
    frequencies = []
    for population in populations:
        valid_values = [
            genotype_values[index]
            for index in pop_freq[population]
            if genotype_values[index] is not None
        ]
        if valid_values:
            frequencies.append(f"{sum(valid_values) / len(valid_values):.3f}")
        else:
            frequencies.append("NA")
    return frequencies


def vcf_to_table_chunks(
    input_file,
    output_prefix,
    pop_file=None,
    chunk_size=10_000,
):
    if not input_file.is_file():
        raise FileNotFoundError(f"Input VCF does not exist: {input_file}")
    if pop_file is not None and not pop_file.is_file():
        raise FileNotFoundError(
            f"Population file does not exist: {pop_file}"
        )

    population_ids = (
        read_population_file(pop_file)
        if pop_file is not None
        else None
    )
    chunk_files = []

    with pysam.VariantFile(input_file) as vcf:
        sample_names = list(vcf.header.samples)
        if population_ids is not None:
            pop_freq = population_indices(population_ids, sample_names)
            output_names = sorted(pop_freq)
        else:
            pop_freq = None
            output_names = sample_names

        chunk_number = 0
        rows_in_chunk = 0
        output_handle = None
        writer = None

        try:
            for record in vcf:
                if record.alts is None or len(record.alts) != 1:
                    continue

                genotype_values = [
                    genotype_value(record.samples[name]["GT"])
                    for name in sample_names
                ]
                observed_values = {
                    value for value in genotype_values if value is not None
                }
                if len(observed_values) <= 1:
                    continue

                if output_handle is None:
                    chunk_number += 1
                    chunk_file, output_handle, writer = open_output_chunk(
                        output_prefix,
                        chunk_number,
                        output_names,
                    )
                    chunk_files.append(chunk_file)

                if pop_freq is not None:
                    output_values = population_frequencies(
                        genotype_values,
                        pop_freq,
                        output_names,
                    )
                else:
                    output_values = [
                        "NA" if value is None else f"{value:g}"
                        for value in genotype_values
                    ]

                writer.writerow(
                    [
                        f"{record.chrom}-{record.pos}",
                        record.chrom,
                        record.pos,
                        record.ref,
                        record.alts[0],
                    ]
                    + output_values
                )
                rows_in_chunk += 1

                if rows_in_chunk == chunk_size:
                    output_handle.close()
                    output_handle = None
                    writer = None
                    rows_in_chunk = 0
        finally:
            if output_handle is not None:
                output_handle.close()

    if not chunk_files:
        raise ValueError("The VCF contains no retained polymorphic variants")
    return chunk_files


def manifest_path(output_prefix):
    return Path(f"{output_prefix}.chunks.tsv")


def write_manifest(output_prefix, chunk_files):
    output_manifest = manifest_path(output_prefix)
    with open(output_manifest, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["task_id", "chunk_file", "matrix_file"])
        for task_id, chunk_file in enumerate(chunk_files, start=1):
            writer.writerow(
                [
                    task_id,
                    chunk_file.resolve(),
                    matrix_output_path(chunk_file).resolve(),
                ]
            )
    return output_manifest


def read_manifest(input_manifest):
    input_manifest = Path(input_manifest)
    if not input_manifest.is_file():
        raise FileNotFoundError(f"Manifest does not exist: {input_manifest}")

    tasks = []
    with open(input_manifest, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"task_id", "chunk_file", "matrix_file"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Manifest must contain columns {sorted(required)}: {input_manifest}"
            )
        for row in reader:
            tasks.append((Path(row["chunk_file"]), Path(row["matrix_file"])))

    if not tasks:
        raise ValueError(f"Manifest contains no tasks: {input_manifest}")
    return tasks


def run_chunk_workers(tasks, jobs):
    if jobs == 1:
        return [process_chunk(*task) for task in tasks]

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(process_chunk, input_file, output_file)
            for input_file, output_file in tasks
        ]
        return [future.result() for future in futures]


def prepare_pipeline(input_file, output_prefix, pop_file, chunk_size):
    chunk_files = vcf_to_table_chunks(
        input_file,
        output_prefix,
        pop_file,
        chunk_size,
    )
    output_manifest = write_manifest(output_prefix, chunk_files)
    print(f"Prepared {len(chunk_files)} chunk(s)")
    print(f"Manifest: {output_manifest}")
    return output_manifest


def join_distance_matrices(matrix_files):
    sample_names_all = None
    available_matrix_all = None
    difference_matrix_all = None

    for matrix_file in matrix_files:
        print(f"Loading: {matrix_file}")
        with np.load(matrix_file) as result:
            required = {
                "sample_names",
                "available_matrix",
                "difference_matrix",
            }
            missing = required.difference(result.keys())
            if missing:
                raise ValueError(
                    f"{matrix_file} is missing objects: {sorted(missing)}"
                )

            sample_names = result["sample_names"]
            available_matrix = result["available_matrix"]
            difference_matrix = result["difference_matrix"]

            if sample_names_all is None:
                sample_names_all = sample_names.copy()
                available_matrix_all = available_matrix.copy()
                difference_matrix_all = difference_matrix.copy()
                continue

            if not np.array_equal(sample_names, sample_names_all):
                raise ValueError(
                    f"Sample names or order differ in {matrix_file}"
                )
            if available_matrix.shape != available_matrix_all.shape:
                raise ValueError(f"Matrix dimensions differ in {matrix_file}")

            available_matrix_all += available_matrix
            difference_matrix_all += difference_matrix

    distance_matrix = np.full(
        difference_matrix_all.shape,
        np.nan,
        dtype=float,
    )
    np.divide(
        difference_matrix_all,
        available_matrix_all,
        out=distance_matrix,
        where=available_matrix_all > 0,
    )
    np.fill_diagonal(distance_matrix, 0.0)

    distance_table = pd.DataFrame(
        distance_matrix,
        index=sample_names_all,
        columns=sample_names_all,
    )
    return (
        distance_table,
        sample_names_all,
        available_matrix_all,
        difference_matrix_all,
        distance_matrix,
    )


def save_distance_outputs(
    output_prefix,
    distance_table,
    sample_names,
    available_matrix,
    difference_matrix,
    distance_matrix,
):
    matrix_file = Path(f"{output_prefix}.mtx")
    distance_table.to_csv(matrix_file, sep="\t", index=True, header=True)

    archive_file = Path(f"{output_prefix}_final.npz")
    np.savez_compressed(
        archive_file,
        sample_names=sample_names,
        available_matrix=available_matrix,
        difference_matrix=difference_matrix,
        distance_matrix=distance_matrix,
    )
    return matrix_file, archive_file


def create_tree_and_pcoa(distance_matrix, sample_names, output_prefix):
    if not np.isfinite(distance_matrix).all():
        raise ValueError(
            "The final distance matrix contains missing or infinite values"
        )

    names = [str(name) for name in sample_names]
    my_dist = DistanceMatrix(distance_matrix, ids=names)
    my_pcoa = pcoa(my_dist)
    my_tree = nj(my_dist, neg_as_zero=False)

    coordinates = my_pcoa.samples
    coordinates.to_csv(
        f"{output_prefix}_pcoa.txt",
        sep="\t",
        index=True,
        index_label="sample",
    )

    figure, axis = plt.subplots(figsize=(6, 6))
    axis.scatter(
        coordinates["PC1"],
        coordinates["PC2"],
        s=30,
        color="black",
        alpha=0.5,
    )
    axis.set_xlabel("PCoA1")
    axis.set_ylabel("PCoA2")
    figure.tight_layout()
    figure.savefig(f"{output_prefix}_pcoa.pdf")
    plt.close(figure)

    newick_file = f"{output_prefix}_tree.newick"
    nexus_file = f"{output_prefix}_tree.nex"
    my_tree.write(newick_file, format="newick")
    Phylo.convert(newick_file, "newick", nexus_file, "nexus")

    with open(f"{output_prefix}_tree_pcoa.pkl", "wb") as handle:
        pickle.dump(
            {
                "sample_names": names,
                "pcoa": my_pcoa,
                "tree": my_tree,
            },
            handle,
        )


def finalize_pipeline(output_prefix, input_manifest):
    tasks = read_manifest(input_manifest)
    missing_outputs = [
        matrix_file
        for _, matrix_file in tasks
        if not matrix_file.is_file()
    ]
    if missing_outputs:
        preview = "\n".join(str(path) for path in missing_outputs[:10])
        remainder = len(missing_outputs) - 10
        if remainder > 0:
            preview += f"\n... and {remainder} more"
        raise FileNotFoundError(
            "Cannot finalize because worker outputs are missing:\n" + preview
        )

    matrix_files = [matrix_file for _, matrix_file in tasks]
    merged = join_distance_matrices(matrix_files)
    save_distance_outputs(output_prefix, *merged)
    create_tree_and_pcoa(merged[-1], merged[1], output_prefix)
    print("Finalization complete.")


def finalize_matrix_files(output_prefix, matrix_files):
    missing_outputs = [
        Path(matrix_file)
        for matrix_file in matrix_files
        if not Path(matrix_file).is_file()
    ]
    if missing_outputs:
        raise FileNotFoundError(
            "Cannot finalize because worker outputs are missing:\n"
            + "\n".join(str(path) for path in missing_outputs)
        )

    merged = join_distance_matrices([Path(path) for path in matrix_files])
    save_distance_outputs(output_prefix, *merged)
    create_tree_and_pcoa(merged[-1], merged[1], output_prefix)
    print("Finalization complete.")


def add_prepare_arguments(parser):
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("-p", "--pop", type=Path)
    parser.add_argument(
        "-w",
        "--window",
        type=positive_integer,
        default=10_000,
        help="Retained variants per table chunk (default: 10000)",
    )


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Convert a VCF to chunked genotype tables, calculate pairwise "
            "distances, merge them, and generate a PCoA and NJ tree."
        )
    )
    subparsers = parser.add_subparsers(dest="stage", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Create genotype-table chunks and a worker manifest",
    )
    add_prepare_arguments(prepare_parser)

    finalize_parser = subparsers.add_parser(
        "finalize",
        help="Merge completed worker outputs and create PCoA/tree outputs",
    )
    finalize_parser.add_argument("-o", "--output", type=Path, required=True)
    finalize_parser.add_argument(
        "-m",
        "--manifest",
        type=Path,
        help="Chunk manifest (default: OUTPUT.chunks.tsv)",
    )

    finalize_files_parser = subparsers.add_parser(
        "finalize-files",
        help="Finalize from an explicit list of worker matrix files",
    )
    finalize_files_parser.add_argument(
        "-o", "--output", type=Path, required=True
    )
    finalize_files_parser.add_argument(
        "--matrices", type=Path, nargs="+", required=True
    )

    run_parser = subparsers.add_parser(
        "run",
        help="Run every stage locally",
    )
    add_prepare_arguments(run_parser)
    run_parser.add_argument(
        "-j",
        "--jobs",
        type=positive_integer,
        default=1,
        help="Concurrent pairwise-distance workers (default: 1)",
    )
    return parser.parse_args(argv)


def main():
    argv = sys.argv[1:]
    known_stages = {"prepare", "finalize", "finalize-files", "run"}
    if argv and argv[0] not in known_stages:
        # Preserve the previous command line by treating it as local run mode.
        argv.insert(0, "run")
    args = parse_arguments(argv)

    if args.stage == "prepare":
        prepare_pipeline(
            args.input,
            args.output,
            args.pop,
            args.window,
        )
        return

    if args.stage == "finalize":
        input_manifest = args.manifest or manifest_path(args.output)
        finalize_pipeline(args.output, input_manifest)
        return

    if args.stage == "finalize-files":
        finalize_matrix_files(args.output, args.matrices)
        return

    print("Stage 1/3: Preparing chunks and manifest")
    output_manifest = prepare_pipeline(
        args.input,
        args.output,
        args.pop,
        args.window,
    )
    tasks = read_manifest(output_manifest)

    print(f"Stage 2/3: Processing {len(tasks)} chunk(s)")
    run_chunk_workers(tasks, args.jobs)

    print("Stage 3/3: Merging matrices and generating PCoA/tree")
    finalize_pipeline(args.output, output_manifest)
    print("Done.")


if __name__ == "__main__":
    main()
