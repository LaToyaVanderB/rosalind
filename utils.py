from Bio import SeqIO

def read_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        return [line.strip() for line in lines]


