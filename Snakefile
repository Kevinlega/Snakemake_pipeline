
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


# Helper functions
# gives all the unhidden files in directory
def listdir(path):
    x = os.listdir(path)
    tmp = []
    for e in x:
        if not e.startswith('.'):
            tmp.append(e)
    return tmp


# function to determine gb needed for Trinity and a dictonary not needed for this pipeline
def GBneeded(fasta_file):
    # contig_dict = {}
    i = 0
    with open(fasta_file, 'rU') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                k = len(record.seq)
                # contig_dict[record.id] = k
                i += k
    # return contig_dict , i
    return i


def compare_sequences(s1 , s2, output):
        
        seqfile1 = screed.open(s2)
        seqfile2 = screed.open(s1)
        

        #lists to hold all the sequences in the files
        f1 , f2 = [] , []

        #load the file into the list
        for read in seqfile1:
            E = sourmash_lib.Estimators(n=20, ksize=3)
            E.add_sequence(read.sequence)
            f1.append((read.name , read.sequence , E))

        #load the other file into the list
        for read in seqfile2:
            E = sourmash_lib.Estimators(n=20, ksize=3)
            E.add_sequence(read.sequence)
            f2.append((read.name , read.sequence , E))


        #output file
        s1filename = os.path.basename(s1)
        s2filename = os.path.basename(s2)

        outputfilename = output + "%s-in-%s.txt"%(s1filename , s2filename) 

        out = open(outputfilename, "w")

        #for every sequence in file A
        for seq1 in f1:

            #set start asuming nothing matches
            max = (None,None,0)

            #for every sequence in file B
            for seq2 in f2:

                #get the jaccard index
                j = seq2[2].jaccard(seq1[2])

                #if the index is larger than the previous maximum
                if j > max[2]:

                    #it becomes the new max secuence
                    max = (seq1[0] , seq2[0] , j )

            #write out the maximum sequence in the output
            out.write("\n%s\n%s\ncorrespondence:%s\n"%(max[0] , max[1] , max[2]))


# Snakemake pipeline start:


import os
import sys          
import re
import math
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator  
from os.path import expanduser
home = expanduser("~")    

# should have screed, sourmash, Trinity, Trimmomatic, TruSeq2-SE.fa (Trimmomatic) and khmer in path

# States the configuration file
configfile: "config.yaml"





Organisms = config['OrganismPATH']
Workdir = config['Workdir']
Workdir = Workdir[0]
ComparisonList = config["ComparisonList"]
TrimPath = config['TrimmomaticPath']
TrimPath = TrimPath[0]
MaxMemory = config['MaxMemory']
MaxMemory = MaxMemory[0]



if len(Organisms) == 1:

    finalout = 'trinity/'

elif len(Organisms) > 1:

    finalout = 'Comparison/'

    try:

        import screed
        import sourmash_lib

    except ImportError:

        sys.exit("screed and sourmash need to be installed and in the path")

else:

    sys.exit("At least one Organism must be present.")


rule all:
    input: finalout


# This rule converts from fasta to fastq 
rule C2Fast:
    input:  Organisms
    output: 'AllFastqs/'
    run:
        for i in range(len(Organisms)):

            # create output directory
            out = Workdir + str(output) +'O' + str(i) + '/'
            if not os.path.isdir(out):
                os.makedirs(out)

            # get the organism files
            directory = listdir(Organisms[i])

            for file in directory:
                # path to file
                filename = Organisms[i] + file

                if filename[-5:] == 'fasta':
                    # eliminates fasta extension and gives fastq 
                    newfilename = file[:-6] + ".fastq"
                    # creates the path to the new file
                    filepath = os.path.join(out,newfilename)

                    # the conversion 
                    with open(filename, "r") as fasta, open(filepath, "w") as fastq:
                        for record in SeqIO.parse(fasta, "fasta"):
                            record.letter_annotations["phred_quality"] = [40] * len(record)
                            SeqIO.write(record, fastq, "fastq")

                elif filename[-5:] == 'fastq':
                    # Moves the fastq to the all fastq folder while ensuring not do destroy original data
                    os.system("cp "+filename+" "+out)
                else:
                    sys.exit("Only FASTA or FASTQ files are allowed for the PIPIELINE!")

# If we want to convert from fastq to fasta we use:
# with open(filename, "rU") as input_handle:
#     with open(filepath, "w+") as output_handle:
#         sequences = SeqIO.parse(input_handle, "fasta")
#         count = SeqIO.write(sequences, output_handle, "fastq")


# This rule uses Trimmomatic to cut adapters which are not needed
rule trim:
    input:  'AllFastqs/'
    output: 'Trimmed/'
    run:
       
        for i in range(len(Organisms)):

            # create output directory
            out = Workdir + str(output) +'O' + str(i) + '/'
            if not os.path.isdir(out):
                os.makedirs(out)

            # get the organism files
            directory = listdir(Organisms[i])
        
            for file in directory:

                # get the file
                ip = Workdir + str(input) + 'O' + str(i) + '/' + file

                # reset output & add filename
                out = Workdir + str(output) +'O' + str(i) + '/'+file

                # execute command
                os.system("java -jar "+TrimPath+"trimmomatic-0.36.jar SE -threads 8 "+ip + " " +out+" ILLUMINACLIP:"+TrimPath+"adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
             

        # remove unneeded directory to clear space
        os.system("rm -r "+ Workdir + str(input))


# This rule normalizes the data
rule Diginorm:
    input: 'Trimmed/'
    output: 'Normalized/'
    run:
        for i in range(len(Organisms)):

            # create output directory
            out = Workdir + str(output) +'O' + str(i) + '/'
            if not os.path.isdir(out):
                os.makedirs(out)

            # get the organism files
            directory = listdir(Organisms[i])

            for file in directory:

                # reset output & add filename
                out = Workdir + str(output) +'O' + str(i) + '/'

                # get the file
                ip = Workdir + str(input) + 'O' + str(i) + '/' + file

                # countgraph needed for filter-abund
                sv = out + "khmer"

                #  the output of normalize-by-median
                op1 = out+ file[:-6] + 'norm.fastq'

                # final output
                op3 = out + file
                
                # add -M memory to use a lot
                os.system("normalize-by-median.py -k 17 -M "+MaxMemory+" -s "+sv + " "+ip + " -o " + op1)
                os.system("filter-abund.py "+sv + " "+op1 + " -o "+ op3 )
                # remove countgraph and first unneeded output
                os.system("rm "+sv)
                os.system("rm "+op1)

        # remove unneeded directory to clear space
        os.system("rm -r "+ Workdir + str(input))

# This rule decides how much CPU is needed for Trinity and creates a file with all the files together
rule GBdecider:
    input:  'Normalized/'
    output: 'Trinityinfo/'
    run:

        for i in range(len(Organisms)):

            # create output directory
            out = Workdir + str(output) +'O' + str(i) + '/'
            if not os.path.isdir(out):
                os.makedirs(out)

            # file with all path to files
            af = out + 'allfiles.txt'
            # file with the GB needed
            filepath = out + 'GB.txt'


            # get the organism files
            directory = listdir(Organisms[i])

            # list to contain all files with path
            allfiles = []

            for file in directory:
                # Get te file
                filename = Workdir + str(input) + 'O' + str(i) + '/' + file
                # add to file list
                allfiles.append(filename)

                fun= GBneeded(filename)

                # creo que not needed
                # filename = os.path.basename(filename)

                gb = fun /(1000000 * 76)
                # write gb to .txt to file
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
            # write all the files in .txt
            for e in range(len(allfiles)-1):
                a = open(af,'a')
                a.write(allfiles[e]+",")
                a.close()
            a = open(af,'a')
            a.write(allfiles[len(allfiles)-1])
            a.close()

            # Round number of GB and add G
            f  = open(filepath, 'r+')
            lines = f.read()
            line = float(lines) 

            if(line > MaxMemory):
                sys.exit("To process data we need "+str(line)+"GB. Look for a computer with more memory. Error on Organism: "+str(i))
                
            line = math.ceil(line)
            line = str(line) + "G"
            f.close()
            f  = open(filepath, 'w+')
            f.write(line)
            f.close()
        
        


# This rule runs the Trinity software explained in technical report
rule Trinity:
    input: 'Trinityinfo/'
    output: 'trinity/'
    run:

        for i in range(len(Organisms)):
            out = Workdir + str(output) 
            if not os.path.isdir(out):
                os.makedirs(out)

            al = Workdir + str(input) + 'O' + str(i) + '/allfiles.txt'
            gb = Workdir + str(input) + 'O' + str(i) + '/GB.txt'

            os.system("Trinity --seqType fq --single `cat "+al+"` --SS_lib_type F --max_memory `cat "+gb+"` --output "+out+" --full_cleanup ")

            out2 =Workdir +'Trinity_Organism' + str(i) + '.fasta'
            output_name = out[:-1]

            os.system("mv "+output_name+".Trinity.fasta"+ " "+out2)
            # os.system("rm -r "+out)
        
        if not os.path.isdir(out):
            os.makedirs(out)

        for i in range(len(Organisms)):
            os.system("mv "+Workdir+'Trinity_Organism'+str(i)+'.fasta '+out)

           
# This rule compares given organisms
rule sourmash:
    input: 'trinity/'
    output: 'Comparison/'
    run:
        for e in ComparisonList:
            try:
                compa = e.split()
            except:
                sys.exit("Nothing to compare in comparion list position: "+e)

            # Get Trinity file for Organisms
            file1 = Workdir + 'Trinity/Trinity_Organism' + str(compa[0]) + '.fasta'
            file2 = Workdir + 'Trinity/Trinity_Organism' + str(compa[1]) + '.fasta'
            out =  Workdir + str(output)
            if not os.path.isdir(out):
                os.makedirs(out)
            # Compare
            compare_sequences(file1,file2,out)

