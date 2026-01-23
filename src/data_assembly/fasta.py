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
        new_seqs = []
        for i, seq in enumerate(self.sequences):
            n_id_min = 0
            n_id_max = len(seq) - 1
            while seq[n_id_min] == "N":
                n_id_min += 1
            while seq[n_id_max] == "N":
                n_id_max -= 1

            new_seq = seq[n_id_min : n_id_max + 1 :]
            if new_seq != 0:
                new_seqs.append(new_seq)
            else:
                print(i, len(new_seq))
                self.titles.pop(i)
        self.sequences = new_seqs

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

    @property
    def in_pgap_range(self) -> bool:
        """Get if the genome len is in ranges of PGAP limits.

        PGAP limits are 10e3 and 100e6.
        """
        genome_len = 0
        for seq in self.sequences:
            genome_len += len(seq)
        return 10e3 < genome_len and genome_len < 100e6

    def remove_all_n_seq(self):
        """Remove sequences with only n in it."""
        for i, seq in enumerate(self.sequences):
            remove = True
            for elem in seq:
                if elem != "N":
                    remove = False
                    break
            if remove:
                self.sequences.pop(i)
                self.titles.pop(i)

    def reduce_successives_n(self, limit=9):
        """Reduce the number of successives N in a sequence.

        The `limit` parameter define the limit when all N after this one will be removed.
        """
        n_id = 0
        new_sequences = []
        for seq in self.sequences:
            new_seq = ""
            for elem in seq:
                if elem == "N":
                    n_id += 1
                else:
                    n_id = 0
                if n_id <= limit:
                    new_seq += elem
            new_sequences.append(new_seq)
        self.sequences = new_sequences
