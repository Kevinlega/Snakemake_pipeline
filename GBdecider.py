# ###############################################################################
#                                                                              #
# Copyright 2017 Megaprobe-Lab                                                 #
#                                                                              #
# This is software created by the megaprobe lab under the GPL3 license.        #
#                                                                              #
# This program reproduces the pipeline for de novo rna sequencing research. To #
# run the program, just move into the folder that contains the makefile and    #
# make sure the sources folder contains the source files for the compilation   #
# of Trinity, khmer, Trimmomatic and ________.                                 #
#                                                                              #
# Then utilize the command: "make" to run the desired operation.               #
# ###############################################################################
#
# aim for ~1G RAM per ~1M ~76 base


import os
import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description = "reads one .fasta files")
parser.add_argument("--file", help = ".fasta or fastq file 1", required = True)
parser.add_argument("--output", help = ".output directory", required = True)

args = parser.parse_args()

filename = args.file

output = args.output

def decide_ext(filenme):
    if filenme.endswith('.fasta'):
        # contig_dict, i = get_contig_dict(filenme)
        i = get_contig_dict(filenme)
        filepath = os.path.join(output,'GBFA.txt')
    elif filenme.endswith('.fastq'):
        filepath = os.path.join(output,'GBFQ.txt')
        # contig_dict, i = get_contig_dict2(filenme)
        i = get_contig_dict2(filenme)
    else:
        sys.exit("File extension not match, only fasta or fastq")
        exit()
    # return contig_dict,i
    return i



def get_contig_dict(fasta_file):
        contig_dict = {}
        i = 0
        with open(fasta_file, 'rU') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    k = len(record.seq)
                    # contig_dict[record.id] = k
                    i += k
        # return contig_dict , i
        return i
def get_contig_dict2(fastq_file):
    contig_dict = {}
    i = 0
    with open(fastq_file,'rU') as handle:
        for record in SeqIO.parse(handle,'fastq'):
            k = len(record.seq)
            # contig_dict[record.id] = k
            i +=k
        # return contig_dict, i
        return i


# arbol1, i = decide_ext(filename) #check file extension
i = decide_ext(filename) #check file extension

# dicoutput = os.path.join(output, filename+'dic.txt')
# GBneeded = os.path.join(output, filename+'GB.txt')



filename = os.path.basename(filename)

if not os.path.exists(output):
    os.makedirs(output)

gb = i /(1000000 * 76)

if os.path.exists(filepath) == False:

    f  = open(filepath, 'w+')
    lines = f.read()
    gb = str(gb)
    f.write(gb)
    f.close()
else:

    f  = open(filepath, 'r+')
    lines = f.read()
    line = float(lines) + float(gb)
    line = str(line)
    f.close()
    f  = open(filepath, 'w+')
    f.write(line)
    f.close()
