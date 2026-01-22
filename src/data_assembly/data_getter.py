"""Utilities function to get data from ifferent online database."""

import subprocess
from pathlib import Path
import shutil
import logging
from assembly_summary import get_genome_accession_from_tsv

INPUT_PATH = Path(__file__).parents[2] / Path("input_data")
OUTPUT_PATH = Path(__file__).parents[2] / Path("output_data")

logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=(Path(__file__).parents[2] / Path("log/data_getter.log")),
    encoding="utf-8",
    level=logging.DEBUG,
)


def get_datasets(gcf_number: str):
    """Downloads dataset from NCBI."""
    gcf_number = gcf_number.replace(" ", "_")
    cmd = f"datasets download genome accession {gcf_number} --include genome,protein,seq-report"
    process = subprocess.Popen(cmd.split(" "), cwd="/tmp", text=False)
    process.wait()
    if process.returncode:
        logger.warning(
            f"Warning the download of the genome {gcf_number} did not succeed correctly."
        )
    else:
        logger.debug(f"Download of the genome {gcf_number} finished correctly.")
        process = subprocess.Popen(
            ["unzip", "/tmp/ncbi_dataset.zip", "-d", "/tmp/ncbi_dataset"]
        )
        process.wait()
        shutil.copytree(
            f"/tmp/ncbi_dataset/ncbi_dataset/data/{gcf_number}",
            OUTPUT_PATH / Path(f"thermococcales/{gcf_number}"),
        )
        shutil.rmtree("/tmp/ncbi_dataset")
        Path("/tmp/ncbi_dataset.zip").unlink()

        if process.returncode:
            logger.warning(
                f"Warning the unzip of the genome {gcf_number} did not succeed correctly."
            )
        else:
            logger.debug(f"Unzip of the genome {gcf_number} finished correctly.")


def get_genomes_prot(tsv_file: Path):
    """Download all genomes and its anotation if available in the tsv file."""
    gcfs = get_genome_accession_from_tsv(tsv_file)
    for gcf in gcfs:
        get_datasets(gcf)


if __name__ == "__main__":
    # get_genomes_prot(INPUT_PATH / Path("alteromonadales.tsv"))
    get_genomes_prot(INPUT_PATH / Path("thermococcales.tsv"))
