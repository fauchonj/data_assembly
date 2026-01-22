"""Utility function to manipulate assembly summary files."""


import csv
from pathlib import Path
from config import INPUT_PATH, OUTPUT_PATH

def filter_assembly_summary(assembly_summary: Path, output_path: Path, genome_accesion: list[str]) -> None:
    """Filter an assembly summary according to some GCF/GCA ID."""
    # print(gcf)
    with assembly_summary.open("r") as f_i, output_path.open("w+") as f_w:
        for line in f_i.readlines()[2:]:
            # print(line.split("\t")[0])
            if line.split("\t")[0] in genome_accesion:
                print(line)

def get_genome_accession_from_tsv(tsv_path: Path) -> list[str]:
    """Get all genome accession from a tsv file."""
    with tsv_path.open("r") as f:
        csv_reader = csv.reader(f, delimiter="\t")
        next(csv_reader)
        return [line[1] for line in csv_reader]

if __name__ == "__main__":
    gcf = get_genome_accession_from_tsv(INPUT_PATH/Path("ncbi_dataset.tsv"))
    print(gcf)
    filter_assembly_summary(INPUT_PATH / Path("archaea_assembly_summary.txt"), OUTPUT_PATH / Path("thermococcales_assembly_summary.txt"), gcf)
