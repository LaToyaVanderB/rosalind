from unittest import TestCase, main

from stronghold import nt_list, nt_comp
from stronghold.core import count_dna_nucleotides, count_dna_nucleotides_file, \
    transcribe_dna_into_rna, complement_dna_strand, rabbits_and_recurrence_relations_recursive, rabbits_and_recurrence_relations_dynamic


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

    def test_transcribe_dna_into_rna(self):
        self.assertEqual(transcribe_dna_into_rna('GATGGAACTTGACTACGTAAATT'),
                         'GAUGGAACUUGACUACGUAAAUU',
                         "Transcribed string should be: GAUGGAACUUGACUACGUAAAUU")

    def test_complement_dna_strand(self):
        self.assertEqual(complement_dna_strand('AAAACCCGGT'),
                         'ACCGGGTTTT',
                         "Complemented strand should be: ACCGGGTTTT")

    def test_rabbits_and_recurrence_relations_recursive(self):
        self.assertEqual(rabbits_and_recurrence_relations_recursive(1, 1, 5, 3),
                         19,
                         "fibonacci(1, 1, 5, 3) =  19")

    def test_rabbits_and_recurrence_relations_dynamic(self):
        self.assertEqual(rabbits_and_recurrence_relations_dynamic(1, 1, 5, 3),
                         19,
                         "fibonacci(1, 1, 5, 3) =  19")

if __name__ == '__main__':
    main()
