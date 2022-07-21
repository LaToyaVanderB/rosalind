from unittest import TestCase, main

from stronghold import nt_list, nt_comp
from stronghold.core import count_dna_nucleotides, count_dna_nucleotides_file, \
    transcribe_dna_into_rna, complement_dna_strand, rabbits_and_recurrence_relations_recursive, \
    highest_gc_content_record, hamming_distance, get_prob_of_dominant, \
    mendels_first_law, fibonacci, mortal_fibo


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

    def test_count_dna_nucleotides_file_multiple_lines(self):
        self.assertEqual(count_dna_nucleotides_file('../data/rosalind_hamm.txt'),
                         '12 12 12 11\n10 15 9 13',
                         "Nucleotide counts should be: '12 12 12 11\n10 15 9 13'")

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

    def test_fibonacci(self):
        self.assertEqual(fibonacci(5, 3),
                         19,
                         "fibonacci(5, 3) =  19")

    def test_mortal_fibo(self):
        self.assertEqual(mortal_fibo(6, 1, 3),
                         4,
                         "mortal_fibo(6, 1, 3) =  4")


    def test_highest_gc_content_record(self):
        record = highest_gc_content_record('../data/fasta.txt')
        self.assertEqual(record.id, 'Rosalind_0808')
        self.assertEqual(round(record.annotations['GC'], 6), 60.919540)

    def test_hamming_distance(self):
        self.assertEqual(hamming_distance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT'), 7)

    def test_mendels_first_law(self):
        self.assertEqual(mendels_first_law(2, 2, 2), 0.78333)

    def test_get_prob_of_dominant(self):
        self.assertEqual(get_prob_of_dominant(29, 28, 25), 0.77582)
        self.assertEqual(get_prob_of_dominant(2, 2, 2), 0.78333)


if __name__ == '__main__':
    main()
