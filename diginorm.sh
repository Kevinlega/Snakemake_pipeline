#!/bin/bash

# Apply digital normalization to the paired-end reads:
#
# /usr/local/share/khmer/scripts/normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savehash normC20k20.kh *.pe.qc.fq.gz
#
# and then to the single-end reads:
#
# /usr/local/share/khmer/scripts/normalize-by-median.py -C 20 --loadhash normC20k20.kh --savehash normC20k20.kh *.se.qc.fq.gz

# Now, run through all the reads and trim off low-abundance parts of high-coverage reads:
#
# /usr/local/share/khmer/scripts/filter-abund.py -V normC20k20.kh *.keep
#
# This will turn some reads into orphans, but that’s ok – their partner read was bad.

# for paired end
# normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savehash normC20k20.kh *.pe.qc.fq.gz





# normalize-by-median.py -C 20 -N 4 -x 15e9 --force_single Files/ZHour.fastq -o ZHournormalize.fastq





for filename in ./Files/* do
  normalize-by-median.py -C 20 -N 4 -x 15e9
done

for filename in ./Files/* do
  filter-abund.py -V normC20k20.kh *.keep
done
