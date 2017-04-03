SAMPLES = ["A", "B"]


rule all:
    input:
        "report.html"


rule concatenate:
    # input:
    #
    # output:
    #
    # run:
        # change to python or run by shell:
        # fasta files
        # for file in ~/Files/*.fasta
        # do
        #   echo $file >> fa.txt | sed '$s/ $/\n/'
        # done
        #
        # cat fa.txt | tr "\n" "," > allFA.txt
        # cat allFA.txt | head -c -1 > allFA.txt
        # rm fa.txt

rule GBdecider:
    input:

    output:

    run:
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





rule converttofasta:
    input:

    output:

    run:
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

rule Diginorm:
    input:

    output:

    run:
        # change to python or run shell:
        #
        # for filename in ./Files/* do
        #   normalize-by-median.py -N 4 -x 15e9
        # done
        #
        # for filename in ./Files/* do
        #   filter-abund.py -V normC20k20.kh *.keep
        # done

rule Trimmomatic:
    input:

    output:

    run:
        #  change to python or run shell:
        # for filename in ./Files/* do
        #   java -jar trimmomatic-0.35.jar SE -phred33 $filename $filename.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        # done
rule Trinity:
    input:

    output:

    run:
        #  change to python or run shell:
        # Trinity --seqType fa --single ` cat ~/allFA.txt`  --SS_lib_type F --CPU `cat ~/Hello/GBFA.txt` --output ~/Hello

rule Sourmash:
    # Courtesy of Angel Sanquiche
    input:

    output:

    run:
        import time

        import sys

        try:

        	import screed
        	import sourmash_lib

        except ImportError:
        	print("screed and sourmash need to be installed and in the path")
        	sys.exit(0)



        def usage():
        	print("Usage:\npython %s <file A> <file B>  "%(sys.argv[0]))
        	sys.exit(0)

        def compare_sequences(s1=None , s2=None):
        	#if imported
        	if s1 and s2:
        		seqfile1 = screed.open(s2)
        		seqfile2 = screed.open(s1)
        	#if ran standalone
        	elif (sys.argv[1] and sys.argv[2]):
        		seqfile1 = screed.open(sys.argv[1])
        		seqfile2 = screed.open(sys.argv[2])
        	#in case of some sort of input error
        	else:
        		usage()

        	#lists to hold all the sequences in the files
        	f1 , f2 = [] , []


        	t1 = time.time()

        	#load the file into the list
        	for read in seqfile1:
        		E = sourmash_lib.Estimators(n=20, ksize=3)
        		E.add_sequence(read.sequence)
        		f1.append((read.name , read.sequence , E))

        	t2 = time.time()

        	print("Loaded file %s in %s seconds"%(sys.argv[1] , t2 - t1))


        	t1 = time.time()

        	#load the other file into the list
        	for read in seqfile2:
        		E = sourmash_lib.Estimators(n=20, ksize=3)
        		E.add_sequence(read.sequence)
        		f2.append((read.name , read.sequence , E))

        	t2 = time.time()

        	print("Loaded file %s in %s seconds"%(sys.argv[2] , t2 - t1))

        	#output file
        	out = open("%s-in-%s.txt"%(sys.argv[1] , sys.argv[2]) , "w")

        	t1 = time.time()

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


        	t2 = time.time()


        	print("Finished comparing all the sequences in %s seconds"%(t2 - t1))
        	print("output:%s-in-%.txts"%(sys.argv[1] , sys.argv[2]))

        if (__name__ == "__main__"):
        	compare_sequences()


# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"
#
#
# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"
#
#
# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"
#
#
# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#     output:
#         "calls/all.vcf"
#     shell:
#         "samtools mpileup -g -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"
#
#
# rule report:
#     input:
#         "calls/all.vcf"
#     output:
#         "report.html"
#     run:
#         from snakemake.utils import report
#         with open(input[0]) as vcf:
#             n_calls = sum(1 for l in vcf if not l.startswith("#"))
#
#         report("""
#         An example variant calling workflow
#         ===================================
#
#         Reads were mapped to the Yeast
#         reference genome and variants were called jointly with
#         SAMtools/BCFtools.
#
#         This resulted in {n_calls} variants (see Table T1_).
#         """, output[0], T1=input[0])
