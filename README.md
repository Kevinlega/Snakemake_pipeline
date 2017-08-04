 	        # ###############################################################################
		#                                                                              #
		# Copyright 2017 Megaprobe-Lab                                                 #
		#                                                                              #
		# This is software created by the megaprobe lab under the GPL3 license.        #
		#                                                                              #
		# This program reproduces the pipeline for de novo rna sequencing research. To #
		# run the program, just move into the folder that contains the Snakefile and   #
		# make sure the PATH has addeded Trinity, khmer, screed, Trimmomatic and       #
		# sourmash. For more explanation on what the program does read technical       #
		# report or to see usage read the README.md                                    #
		#                                                                              #
		# Then utilize the command: "snakemake" to run the desired operation.          #
		# ###############################################################################


# Programs:
<!-- Note: When creating enviorment use python3 not python2 -->
- Environment and khmer: http://khmer.readthedocs.io/en/v1.0.1/install.html 
- Trinity: https://github.com/trinityrnaseq/trinityrnaseq/releases, https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity (v2.4)
- Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic (v0.36)
- Sourmash: https://sourmash.readthedocs.io/en/latest/tutorials.html#installing-sourmash
- Screed: https://pypi.python.org/pypi/screed
- Biopython: https://anaconda.org/anaconda/biopython
- Anaconda: https://anaconda.org/
- Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
<!-- Once you have the enviorment you can do pip install <package> and add different packages -->

# App requirments:
<!-- For Trinity and Trimmomatic.  -->
- Java-1.8 (or higher) : sudo apt-get install default-jre
- gcc greater than 4.3: sudo apt-get install gcc
- zlib: sudo apt-get install zlib1g-dev
- bowtie2: sudo apt-get install bowtie2 or http://bowtie-bio.sourceforge.net/bowtie2/index.shtml


# Pip requirments:
- pip install PyYAML
- pip install screed
- pip install khmer or with conda install --name (name of enviorment) -c bioconda khmer=2.1
- pip install biopython
- pip install snakemake

# Export Programs to Path:
- Like so: PATH=$PATH:~/path/to/file
- Trimmomatic 
- Trinity
- Sourmash
- Screed

# Modify Config.yaml
- Add working directory
- Add path to Organisms 
- Add Comparison like so: "0 1" This is Organisms 1 and 2.


## Run the snakemake
- Have the Snakefile and config.yaml in the working directory
- run like so snakemake
- or snakemake <rule>
- for other uses use snakemake -h

## Command to see how many sequences:
- grep -o 'TRINITY_DN' Trinity_Organism0.fasta-in-Trinity_Organism1.fasta.txt | wc -l 

