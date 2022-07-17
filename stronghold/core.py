import sys
from stronghold import nt_list, nt_comp
from utils import read_file


def count_dna_nucleotides(dna_string: str) -> dict:
    """
    Count each nucleotide in a DNA string

    Parameters
    ----------
    dna_string: str
        The input DNA string.

    Returns
    -------
    dict
        The counts for A, C, T, G.
    """
    counts = {nt: dna_string.count(nt) for nt in nt_list}
    return counts


def count_dna_nucleotides_file(filename: str) -> str:
    """
    Count each nucleotide in a DNA string

    Parameters
    ----------
    filename: str
        The file containing the input DNA string.

    Returns
    -------
        A string containing the counts for A, C, T, G, in this order, space-delimited.
    """

    counts = count_dna_nucleotides(read_file(filename))
    return " ".join([str(counts[nt]) for nt in nt_list])


def transcribe_dna_into_rna(dna_string: str) -> str:
    """
    Transcribe a DNA string into an RNA string, i.e. replace 'T' by 'U'

    Parameters
    ----------
    dna_string: str
        The input DNA string.

    Returns
    -------
        The transcribed string.
    """

    return dna_string.replace('T', 'U')


def complement_dna_strand(dna_string: str) -> str:
    """
    Complement a DNA strand into an RNA string, i.e.
        1. reverse the string
        2. replace nucleotides by their complement bases (A <-> T, C <-> G)

    Parameters
    ----------
    dna_string: str
        The input DNA strand.

    Returns
    -------
        The complement strand.
    """

    table = str.maketrans(nt_comp)
    return dna_string[::-1].translate(table)


def rabbits_and_recurrence_relations_dynamic(f1, f2, n, k):
    print("What's up, Doc?")
    return fibo_dynamic(f1, f2, n, k)


def rabbits_and_recurrence_relations_recursive(f1, f2, n, k):
    print("What's up, Doc?")
    computed = {}
    return fibo_recursive(f1, f2, n, k, computed)


def fibo_recursive(f1, f2, n, k, computed):
    args = (n, k)
    if args in computed:
        return computed[args]

    if n == 1:
        res = f1
    elif n == 2:
        res = f2
    else:
        res = fibo_recursive(f1, f2, n-1, k, computed) + k * fibo_recursive(f1, f2, n-2, k, computed)
    computed[args] = res
    return res


def fibo_dynamic(f1, f2, n, k):
    total_n_2 = f1
    total_n_1 = f2
    total_n = 0
    for months in range(1, n-1):
        total_n = total_n_1 + k * total_n_2
        total_n_2 = total_n_1
        total_n_1 = total_n
    return total_n


if __name__ == "__main__":
    sys.exit(count_dna_nucleotides_file(sys.argv[1]))
