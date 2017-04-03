import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description = "reads one .fasta files")
parser.add_argument("--file", help = ".fasta or fastq file 1", required = True)
parser.add_argument("--output", help = ".output directory", required = True)


args = parser.parse_args()
filename = args.file
output = args.output

newfilename = os.path.basename(filename)
newfilename = newfilename[:-6]
newfilename += ".fasta"

filepath = os.path.join(output,newfilename)

with open(filename, "rU") as input_handle:
    with open(filepath, "w+") as output_handle:
        sequences = SeqIO.parse(input_handle, "fastq")
        count = SeqIO.write(sequences, output_handle, "fasta")
