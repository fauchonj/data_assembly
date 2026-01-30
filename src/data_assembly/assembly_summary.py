"""Utility function to manipulate assembly summary files."""

import csv
from pathlib import Path
from data_assembly.config import INPUT_PATH, OUTPUT_PATH


def filter_assembly_summary(assembly_summary: Path, output_path: Path) -> None:
    """Filter an assembly summary to keep onlly GCF lines.

    When downloading assembly summary from NCBI their is multiple line for
    the same genomes. This function rewrite the assembly summary to keep
    only one line by accession nuber.
    """
    with assembly_summary.open("r") as f_i, output_path.open("w+") as f_w:
        csv_reader = csv.reader(f_i, delimiter="\t")
        next(csv_reader)
        csv_writter = csv.writer(f_w, delimiter="\t")
        for line in csv_reader:
            if line[1][0:3] == "GCA" and line[2] != "":
                continue
            else:
                csv_writter.writerow(line)


def get_column_from_tsv(tsv_path: Path, column_id: int) -> list[str]:
    """Get all genome accession from a tsv file."""
    with tsv_path.open("r") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        return [line[1] for line in csv_reader]


if __name__ == "__main__":
    # filter_assembly_summary(
    #     Path(__file__).parents[2] / Path("input_data/thermococcales.tsv"),
    #     Path(__file__).parents[2] / Path("input_data/thermococcales_filtered.tsv"),
    # )
    # filter_assembly_summary(
    #     Path(__file__).parents[2] / Path("input_data/alteromonadales.tsv"),
    #     Path(__file__).parents[2] / Path("input_data/alteromonadales_filtered.tsv"),
    # )
