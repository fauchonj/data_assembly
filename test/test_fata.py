"""Test the module FASTA file."""

from pathlib import Path
from data_assembly.fasta import Fasta
import pytest


@pytest.mark.parametrize(
    "fasta_path, supposed_sequences, supposed_titles",
    [
        (
            Path(__file__).parent / Path("inputs/basic.fst"),
            ["AAAA"],
            [" Basic fasta file"],
        ),
        (
            Path(__file__).parent / Path("inputs/long_line_fasta.fst"),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            ],
            [" Long line fasta"],
        ),
        (
            Path(__file__).parent / Path("inputs/multi_line.fst"),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            ],
            [" Line 0", " Line 1", " Line 2"],
        ),
    ],
)
def test_parsing_fasta(
    fasta_path: Path, supposed_sequences: list, supposed_titles: list
):
    """Check if the classmethod `from_fasta` parse correctly a fasta file."""
    fasta_file = Fasta.from_fasta_file(fasta_path)
    assert fasta_file.sequences == supposed_sequences
    assert fasta_file.titles == supposed_titles


@pytest.mark.parametrize(
    "fasta_path, supposed_sequences, supposed_titles",
    [
        (
            Path(__file__).parent / Path("inputs/seq_n_start_n_end.fst"),
            ["AB", "AB", "ABC"],
            ["N start", "N end", "No N"],
        ),
        (
            Path(__file__).parent / Path("inputs/seq_too_short.fst"),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
            ],
            ["Keep", "Remove", "Keep", "Keep", "Remove"],
        ),
    ],
)
def test_n_first_last(
    fasta_path: Path, supposed_sequences: list, supposed_titles: list
):
    """Check the function removing N."""
    fasta_file = Fasta.from_fasta_file(fasta_path)
    fasta_file.remove_first_last_n()
    assert fasta_file.sequences == supposed_sequences
    assert fasta_file.titles == supposed_titles


@pytest.mark.parametrize(
    "fasta_path, supposed_sequences, supposed_titles",
    [
        (
            Path(__file__).parent / Path("inputs/seq_n_start_n_end.fst"),
            [],
            [],
        ),
        (
            Path(__file__).parent / Path("inputs/seq_too_short.fst"),
            [
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
                "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
            ],
            ["Keep", "Keep", "Keep"],
        ),
    ],
)
def test_lines_too_short(
    fasta_path: Path, supposed_sequences: list, supposed_titles: list
):
    """Check the function removing N."""
    fasta_file = Fasta.from_fasta_file(fasta_path)
    fasta_file.remove_seq_too_short()
    assert fasta_file.sequences == supposed_sequences
    assert fasta_file.titles == supposed_titles


@pytest.mark.parametrize(
    "input_path",
    [
        Path(__file__).parent / Path("inputs/basic.fst"),
        Path(__file__).parent / Path("inputs/long_line_fasta.fst"),
        Path(__file__).parent / Path("inputs/multi_line.fst"),
        Path(__file__).parent / Path("inputs/seq_n_start_n_end.fst"),
        Path(__file__).parent / Path("inputs/seq_too_short.fst"),
    ],
)
def test_write_faste(input_path: Path):
    """Test the writting of fasta in a file."""
    fasta_file = Fasta.from_fasta_file(input_path)
    fasta_file.to_fasta_file(
        Path(__file__).parent / Path(f"results/{input_path.stem}{input_path.suffix}")
    )
    assert 0 == 0
