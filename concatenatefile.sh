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

# Trinity needs fasta and fastq files in different runs, so we will need to
#  separate them


#!/bin/bash


# fastq files
for file in ~/Files/*.fastq
do
  echo $file >> fq.txt | sed '$s/ $/\n/'
done


cat fq.txt | tr "\n" "," > allFQ.txt
cat allFQ.txt | head -c -1 >> allFQ.txt

rm fq.txt


# fasta files
for file in ~/Files/*.fasta
do
  echo $file >> fa.txt | sed '$s/ $/\n/'
done

cat fa.txt | tr "\n" "," > allFA.txt
cat allFA.txt | head -c -1 > allFA.txt
rm fa.txt
