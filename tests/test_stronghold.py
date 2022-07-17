from unittest import TestCase, main

from stronghold import nt_list, nt_comp
from stronghold.core import count_dna_nucleotides, count_dna_nucleotides_file


class StrongholdTestCase(TestCase):
    def test_constants(self):
        self.assertEqual(nt_list,
                         ['A', 'C', 'G', 'T'],
                         "Nucleotide list should be ['A', 'C', 'G', 'T']")
        self.assertEqual(nt_comp,
                         {"A": "T", "C": "G", "G": "C", "T": "A"},
                         "Nucleotide complements should be {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}")

    def test_count_dna_nucleotides(self):
        self.assertEqual(count_dna_nucleotides('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'),
                         {'A': 20, 'C': 12, 'G': 17, 'T': 21},
                         "Nucleotide counts should be: {'A': 20, 'C': 12, 'G': 17, 'T': 21}")

    def test_count_dna_nucleotides_file(self):
        self.assertEqual(count_dna_nucleotides_file('../data/count_dna_nucleotides.txt'),
                         '20 12 17 21',
                         "Nucleotide counts should be: {'A': 20, 'C': 12, 'G': 17, 'T': 21}")


if __name__ == '__main__':
    main()
