"""Module to manipulate fasta files."""

import logging
from pathlib import Path
from typing import Self

logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=(Path(__file__).parents[2] / Path("log/fasta.log")),
    encoding="utf-8",
    level=logging.DEBUG,
)


class Fasta:
    """Represent a FASTA file."""

    def __init__(self, sequences: list, titles: list, stem: str | None = None):
        self.sequences: list = sequences
        self.titles: list = titles
        self.stem: str = stem

    @classmethod
    def from_fasta_file(cls, fasta_path: Path) -> Self:
        """Parse a fasta file and create an instance of Self."""
        titles = []
        sequences = []
        seq = ""
        with fasta_path.open("r") as f:
            for line in f.readlines():
                if ">" in line:
                    if seq != "":
                        sequences.append(seq)
                    titles.append(line.replace("\n", "")[1:])
                    seq = ""
                else:
                    seq += line.replace("\n", "")
        if seq != "":
            sequences.append(seq)
        return cls(sequences, titles, fasta_path.stem)

    def remove_seq_too_short(self, threshold: int = 200):
        """Remove sequences with length less than `threshold`.

        PGAP input files need to have sequences with more than 200 nucleotides.
        """
        sequences = []
        titles = []
        for i, seq in enumerate(self.sequences):
            if len(seq) >= threshold:
                sequences.append(seq)
                titles.append(self.titles[i])
            else:
                logger.debug(
                    f"Remove sequence: ({self.titles[i]}, {seq}) because it's too short ({threshold})."
                )
        self.sequences = sequences
        self.titles = titles

    def remove_first_last_n(self):
        """Remove first and last nucleotides in sequences being N."""
        for i, seq in enumerate(self.sequences):
            if seq[0] == "N":
                self.sequences[i] = seq[1:]
            if seq[-1] == "N":
                self.sequences[i] = seq[:-1]

    def to_fasta_file(self, output_path: Path):
        """Write the fasta in a fasta file."""
        with output_path.open("w+") as f:
            for title, seq in zip(self.titles, self.sequences, strict=True):
                title = f">{title}"
                sequences = [
                    seq[i : min(i + 50, len(seq))] + "\n"
                    for i in range(0, len(seq), 50)
                ]
                f.writelines([title + "\n", *sequences])
