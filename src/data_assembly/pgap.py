"""Utilities to use pgap to annotate a genome."""

from pathlib import Path
import subprocess
import csv
import tqdm
from shlex import quote
import logging
from data_assembly.fasta import Fasta
from multiprocessing import Pool
from data_assembly.assembly_summary import get_column_from_tsv

logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=(Path(__file__).parents[2] / Path("log/pgap.log")),
    encoding="utf-8",
    level=logging.DEBUG,
)


def get_genome_file_from_accession(
    genome_path: Path, genome_accession: str
) -> Path | None:
    """Get the genome file from its genome accession.

    Genome accesion must be include in file name.
    """
    for path in genome_path.iterdir():
        if genome_accession in path.stem:
            return path
    return None


def get_pgap_inputs(
    dir_genomes: Path, dir_output: Path, tsv_path: Path
) -> list[tuple[Path, Path, str]]:
    """Get inputs of the pgap command from a tsv files and input to search genome in."""
    pgap_inputs = []
    with tsv_path.open("r") as tsv_file:
        for line in tsv_file.readlines():
            genome_accession = line[1]
            org_name = line[3]
            genome_path = get_genome_file_from_accession(dir_genomes, genome_accession)
            if genome_path:
                pgap_inputs.append(
                    (genome_path, (dir_output / Path(f"{genome_path.stem}")), org_name)
                )
    return pgap_inputs


def run_pgap(args: tuple[Path, Path, str]):
    """Run pgap."""
    (input_path, output_path, org_name) = args
    cmd = f"/home/pgap/pgap.py -n -d -g {str(input_path)} -o {str(output_path)} -s {quote(org_name)}"

    process = subprocess.Popen(cmd, text=True, shell=True)
    returncode = process.wait()

    if returncode:
        logger.warning(
            f"PGAP runned on the file {str(input_path)} but did not success."
        )
    else:
        logger.debug(f"PGAP runned on the file {str(input_path)} and success.")
    return 0


def parse_genome(args: tuple[Path, Path]):
    """Parse a genome and create its correct version for pgap."""
    (genome_path, parsed_dir) = args
    fasta_file = Fasta.from_fasta_file(genome_path)
    fasta_file.remove_first_last_n()
    fasta_file.remove_seq_too_short()
    fasta_file.to_fasta_file(
        (parsed_dir / Path(f"{genome_path.stem}{genome_path.suffix}"))
    )


if __name__ == "__main__":
    pgap_inputs = get_pgap_inputs(
        Path("/data/pgap/parsed_genomes/thermococcales"),
        Path("/data/pgap/proteomes/thermococcales"),
        Path(__file__).parents[2] / Path("input_data/thermococcales.tsv"),
    )
    with Pool(8) as p:
        tqdm.tqdm(
            p.imap(
                run_pgap,
                pgap_inputs,
            ),
            total=len(pgap_inputs),
        )
    # with Pool(8) as p:
    #     p.map(
    #         run_pgap,
    #         get_pgap_inputs(
    #             Path("/data/pgap/parsed_genomes/alteromonadales"),
    #             Path("/data/pgap/proteomes/alteromonadales"),
    #         ),
    #     )
