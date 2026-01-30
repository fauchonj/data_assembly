"""Utilities to use pgap to annotate a genome."""

from pathlib import Path
import subprocess
import csv
import logging
from parallelbar import progress_map
import yaml

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
    paths = []
    with tsv_path.open("r") as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        next(tsv_reader)
        for line in tsv_reader:
            genome_accession = line[1]
            org_name = line[3]
            strain = line[7]
            genome_path = get_genome_file_from_accession(
                dir_genomes, "GCR_" + genome_accession.split("_")[1]
            )
            if genome_path and genome_path not in paths:
                pgap_inputs.append((genome_path, strain, org_name))
                paths.append(genome_path)
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
    with (Path(__file__).parents[2] / Path("templates/template_pgap.yaml")).open(
        "r"
    ) as f:
        input_yaml = yaml.safe_load(f)

    with (Path(__file__).parents[2] / Path("templates/template_submol.yaml")).open(
        "r"
    ) as f:
        submol_yaml = yaml.safe_load(f)

    Path(f"/tmp/{genome_path.stem}").mkdir(exist_ok=False)
    genome_path.copy(
        Path(f"/tmp/{genome_path.stem}/{genome_path.stem}{genome_path.suffix}")
    )

    input_yaml["fasta"]["location"] = str(
        Path(f"/tmp/{genome_path.stem}/{genome_path.stem}{genome_path.suffix}")
    )
    input_yaml["submol"]["location"] = f"/tmp/{genome_path.stem}/submol.yaml"

    submol_yaml["organism"]["genus_species"] = genus_species
    submol_yaml["organism"]["strain"] = strain

    with Path(f"/tmp/{genome_path.stem}/submol.yaml").open("w+") as f:
        yaml.safe_dump(submol_yaml, f)

    with Path(f"/tmp/{genome_path.stem}/input.yaml").open("w+") as f:
        yaml.safe_dump(input_yaml, f)


def create_imput_and_run_pgap(args):
    """Create PGAP yaml input file and run it."""

    (genome_path, strain, genus_specied, output_path) = args
    create_input_pgap(
        genome_path=genome_path, genus_species=genus_specied, strain=strain
    )
    run_pgap(output_path, Path(f"/tmp/{genome_path.stem}/input.yaml"))
    for path in Path(f"/tmp/{genome_path.stem}").iterdir():
        path.unlink()
    Path(f"/tmp/{genome_path.stem}").rmdir()


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
    )
    pgap_inputs = [
        (
            genome_path,
            strain,
            org_name,
            Path(f"/data/pgap/proteomes/thermococcales/{genome_path.stem}"),
        )
        for (genome_path, strain, org_name) in pgap_inputs
    ]
    progress_map(create_imput_and_run_pgap, pgap_inputs, n_cpu=8)

    pgap_inputs = get_pgap_inputs(
        Path("/data/pgap/parsed_genomes/alteromonadales"),
        Path(__file__).parents[2] / Path("input_data/alteromonadales.tsv"),
    )
    pgap_inputs = [
        (
            genome_path,
            strain,
            org_name,
            Path(f"/data/pgap/proteomes/alteromonadales/{genome_path.stem}"),
        )
        for (genome_path, strain, org_name) in pgap_inputs
    ]
    progress_map(create_imput_and_run_pgap, pgap_inputs, n_cpu=8)
