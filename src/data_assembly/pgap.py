"""Utilities to use pgap to annotate a genome."""

from pathlib import Path
import subprocess
import csv
import tqdm
from shlex import quote
import logging
from data_assembly.fasta import parse_genome
from multiprocessing import Pool
from data_assembly.assembly_summary import get_column_from_tsv
from parallelbar import progress_map
from ruyaml import YAML
from ruyaml.scalarstring import SingleQuotedScalarString

yaml = YAML()
yaml.default_flow_style = "'"
# yaml.default_style = "'"
yaml.preserve_quotes = True
yaml.explicit_start = True
yaml.explicit_end = True
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


def get_pgap_inputs(dir_genomes: Path, tsv_path: Path) -> list[tuple[Path, Path, str]]:
    """Get inputs of the pgap command from a tsv files and input to search genome in."""
    pgap_inputs = []
    with tsv_path.open("r") as tsv_file:
        for line in tsv_file.readlines():
            genome_accession = line[1]
            org_name = line[3]
            strain = line[7]
            genome_path = get_genome_file_from_accession(dir_genomes, genome_accession)
            if genome_path:
                pgap_inputs.append((genome_path, strain, org_name))
    return pgap_inputs


def run_pgap(output_path: Path, input_yaml: Path):
    """Run pgap."""
    cmd = f"/home/pgap/pgap.py -n -d -o {str(output_path)} {str(input_yaml)}"

    process = subprocess.Popen(cmd, text=True, shell=True)
    returncode = process.wait()

    if returncode:
        logger.warning(
            f"PGAP runned on the file {output_path.stem} but did not success."
        )
    else:
        logger.debug(f"PGAP runned on the file {output_path.stem} and success.")


def create_input_pgap(genome_path: Path, genus_species: str, strain: str):
    """Create PGAP input yaml file and its submol."""
    with (Path(__file__).parents[2] / Path("input_data/template_pgap.yaml")).open(
        "r"
    ) as f:
        input_yaml = yaml.load(f)
    with (Path(__file__).parents[2] / Path("input_data/template_submol.yaml")).open(
        "r"
    ) as f:
        submol_yaml = yaml.load(f)

    Path(f"/tmp/{genome_path.stem}").mkdir(exist_ok=False)
    genome_path.copy(
        Path(f"/tmp/{genome_path.stem}/{genome_path.stem}{genome_path.suffix}")
    )

    input_yaml["fasta"]["location"] = str(
        Path(f"/tmp/{genome_path.stem}/{genome_path.stem}{genome_path.suffix}")
    )
    input_yaml["submol"]["location"] = f"/tmp/{genome_path.stem}/submol.yaml"

    submol_yaml["organism"]["genus_species"] = SingleQuotedScalarString(genus_species)
    submol_yaml["organism"]["strain"] = SingleQuotedScalarString(strain)

    with Path(f"/tmp/{genome_path.stem}/submol.yaml").open("w+") as f:
        yaml.dump(submol_yaml, f)

    with Path(f"/tmp/{genome_path.stem}/input.yaml").open("w+") as f:
        yaml.dump(input_yaml, f)


def create_imput_and_run_pgap(args):
    """Create PGAP yaml input file and run it."""

    (genome_path, strain, genus_specied, output_path) = args
    create_input_pgap(
        genome_path=genome_path, genus_species=genus_specied, strain=strain
    )
    run_pgap(output_path, Path(f"/tmp/{genome_path.stem}/input.yaml"))


if __name__ == "__main__":
    # input_path = [
    #     (path, Path("/data/pgap/parsed_genomes/thermococcales/"))
    #     for path in Path("/data/pgap/input_genomes/thermococcales/").iterdir()
    # ]
    # progress_map(parse_genome, input_path)

    # input_path = [
    #     (path, Path("/data/pgap/parsed_genomes/alteromonadales/"))
    #     for path in Path("/data/pgap/input_genomes/alteromonadales/").iterdir()
    # ]
    # progress_map(parse_genome, input_path)

    pgap_inputs = get_pgap_inputs(
        Path("/data/pgap/parsed_genomes/thermococcales"),
        Path(__file__).parents[2] / Path("input_data/thermococcales.tsv"),
    )[0:16:]
    pgap_inputs = [
        (
            genome_path,
            strain,
            org_name,
            Path(f"/home/pgap/proteomes/thermococcales/{genome_path.stem}"),
        )
        for (genome_path, strain, org_name) in pgap_inputs
    ]
    print(pgap_inputs)
    progress_map(create_imput_and_run_pgap, pgap_inputs)
    # with Pool(8) as p:
    #     tqdm.tqdm(
    #         p.imap(
    #             run_pgap,
    #             pgap_inputs,
    #         ),
    #         total=len(pgap_inputs),
    #     )
    # with Pool(8) as p:
    #     p.map(
    #         run_pgap,
    #         get_pgap_inputs(
    #             Path("/data/pgap/parsed_genomes/alteromonadales"),
    #             Path("/data/pgap/proteomes/alteromonadales"),
    #         ),
    #     )
