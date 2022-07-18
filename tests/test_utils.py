import unittest

from utils import read_file


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(read_file('../data/count_dna_nucleotides.txt'),
                         ['AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'],
                         "File should contain: 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'")

        self.assertEqual(read_file('../data/rosalind_hamm.txt'),
                         ['AGGCCCTAATGTAACCTCTTCTGGCGCCGCGAAAGTTAGATTGACGA', 'CTCCCATGACATTAACTCCGTATTGGCCGCGAATGTTCCATGGACCT'],
                         "File should contain: 'AGGCCCTAATGTAACCTCTTCTGGCGCCGCGAAAGTTAGATTGACGA', 'AGGCCCTAATGTAACCTCTTCTGGCGCCGCGAAAGTTAGATTGACGA'")


if __name__ == '__main__':
    unittest.main()

