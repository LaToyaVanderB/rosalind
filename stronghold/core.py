import sys
from stronghold import nt_list, nt_comp
from utils import read_file
from Bio import SeqIO
from operator import ne
from scipy.special import comb


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

    lines = read_file(filename)
    counts = [count_dna_nucleotides(line) for line in lines]

    def flatten(dict_counts):
        return " ".join([str(dict_counts[nt]) for nt in nt_list])
    counts_flat = map(flatten, counts)

    return "\n".join(counts_flat)


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


def rabbits_and_recurrence_relations_recursive(f1, f2, n, k):
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


def fibonacci(n: int, k: int) -> int:
    """
    Fibonacci sequence:
        We start with one pair of rabbits.
        Rabbits need one month to become mature.
        Mature rabbit pairs produce k new rabbit pairs every month.
        Rabbits don´t die or stop reproducing.
    :param n: total number of months
    :param k: how many new pairs a mature pair produces each month
    :return: how many pairs after n months
    """
    total_n_2 = 1
    total_n_1 = 1
    total_n = 0
    for months in range(1, n-1):
        total_n = total_n_1 + k * total_n_2
        total_n_2 = total_n_1
        total_n_1 = total_n
    return total_n



def mortal_fibo(n: int, k: int, m: int) -> int:
    """
    Fibonacci sequence where rabbits die after m months:
        We start with one pair of rabbits.
        Rabbits need one month to become mature.
        Mature rabbit pairs produce k new rabbit pairs every month.
        Rabbits die after m months.

    :param n: total number of months
    :param k: how many new pairs a mature pair produces each month
    :param m: how many months a pair live
    :return: how many pairs after n months
    """

    t = 1
    pop = [1] + [0]*m #F1
    # print(t, pop, sum(pop[0:-1]))
    for t in range(2, n+1):
        pop = pop[:-1]
        new_b_t = sum(pop[1:])
        pop = [new_b_t] + pop
        # print(t, pop, sum(pop[0:-1]))
    return sum(pop[0:-1])


def gc_content(s):
    return (s.count('C') + s.count('G')) / len(s) * 100


def highest_gc_content_record(filename: str) -> str:
    it = SeqIO.parse(filename, "fasta")
    max_gc = 0
    max_gc_record = None
    for r in it:
        gc = gc_content(r.seq)
        if gc > max_gc:
            max_gc_record = r
            max_gc = gc
    max_gc_record.annotations['GC'] = max_gc
    return max_gc_record


def hamming_distance(s1, s2):
    return sum(map(ne, s1, s2))


def mendels_first_law(k, m, n):
    """
        Probability of at least a dominant gene at the next generation
        :param k: number of homozygous dominant individuals
        :param m: number of heterozygous individuals
        :param n: number of homozygous recessive individuals
        :return: number of offspring with a dominant allele

        THIS DOES NOT  RETURN WHAT ROSALIND WANTS
        """
    total_pop = k+m+m
    # total population

    total_pairs = comb(total_pop, 2)
    # Combinations of 2 among total population

    # Pairs with a dominant allele =
    #     Dom-Dom pairs: k choose 2
    #   + Dom-Het pairs: k*m
    #   + Dom-Rec pairs: k*n
    #   + Het-Het pairs with the allele: m choose 2 * .75
    #   + Het-Rec pairs with the allele: m*n * .5

    pairs_with_a_dominant_allele = comb(k, 2) + k*m + k*n + comb(m, 2)*.75 + m*n*.5
    return round(pairs_with_a_dominant_allele / total_pairs, 5)


def get_prob_of_dominant(k, m, n):
    # A - dominant factor
    # a - recessive factor
    # k - amount of organisms with AA factors (homozygous dominant)
    # m - amount of organisms with Aa factors (heterozygous)
    # n - amount of organisms with aa factors (homozygous recessive)
    events = ['AA+Aa', 'AA+aa', 'Aa+aa', 'AA+AA', 'Aa+Aa', 'aa+aa']

    # get the probability of dominant traits (set up Punnett square)
    punnett_probabilities = {
        'AA+Aa': 1,
        'AA+aa': 1,
        'Aa+aa': 1 / 2,
        'AA+AA': 1,
        'Aa+Aa': 3 / 4,
        'aa+aa': 0,
    }
    event_probabilities = {}
    totals = k + m + n

    # Event: AA+Aa -> P(X=k, Y=m) + P(X=m, Y=k):
    p_km = k / totals * m / (totals - 1)
    p_mk = m / totals * k / (totals - 1)
    event_probabilities['AA+Aa'] = p_km + p_mk

    # Event: AA+aa -> P(X=k, Y=n) + P(X=n, Y=k):
    p_kn = k / totals * n / (totals - 1)
    p_nk = n / totals * k / (totals - 1)
    event_probabilities['AA+aa'] = p_kn + p_nk

    # Event: Aa+aa -> P(X=m, Y=n) +P(X=n, Y=m):
    p_mn = m / totals * n / (totals - 1)
    p_nm = n / totals * m / (totals - 1)
    event_probabilities['Aa+aa'] = p_mn + p_nm

    # Event: AA+AA -> P(X=k, Y=k):
    p_kk = k / totals * (k - 1) / (totals - 1)
    event_probabilities['AA+AA'] = p_kk

    # Event: Aa+Aa -> P(X=m, Y=m):
    p_mm = m / totals * (m - 1) / (totals - 1)
    event_probabilities['Aa+Aa'] = p_mm

    # Event: aa+aa -> P(X=n, Y=n) + P(X=n, Y=n) = 0 (will be * 0, so just don't use it)
    event_probabilities['aa+aa'] = n / totals * (n-1) / (totals - 1)

    # Total probability is the sum of (prob of dominant factor * prob of the event)
    total_probability = 0
    for event in events:
        total_probability += punnett_probabilities[event] * event_probabilities[event]
    return round(total_probability, 5)


if __name__ == "__main__":
    # sys.exit(fibonacci(5, 3))
    sys.exit(mortal_fibo(100, 1, 16))
    # sys.exit(mortal_fibo(6, 1, 3))
