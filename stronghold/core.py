import sys
from stronghold import nt_list
from utils import read_file


def count_dna_nucleotides(dna_string):
    """
    Count each nucleotide in a DNA string

    Parameters
    ----------
    dna_string: str
        The input DNA string.

    Returns
    -------
    dict
        The counts for A, C, T, G
    """
    counts = {nt: dna_string.count(nt) for nt in nt_list}
    return counts


def count_dna_nucleotides_file(filename):
    counts = count_dna_nucleotides(read_file(filename))
    return " ".join([str(counts[nt]) for nt in nt_list])


if __name__ == "__main__":
    sys.exit(count_dna_nucleotides_file(sys.argv[1]))
