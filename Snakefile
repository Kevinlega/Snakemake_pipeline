SAMPLE = config['Sample']
CONVERT = config['converttofast']

rule all:

rule C2Fast:
    input:   conv = expand('Files/{sample}', sample=CONVERT))
    output: 'AllFastas/'
    run:
        import os
        import sys
        import argparse
        from Bio import SeqIO
        from Bio.SeqIO.QualityIO import FastqGeneralIterator

        filename = str(conv)
        print(filename)
        if filename[-5:] == 'fastq':
            newfilename = os.path.basename(filename)
            newfilename = newfilename[:-6]
            newfilename += ".fasta"

            print(newfilename)
            x = str(output)
            filepath = os.path.join(x,newfilename)

            print(filepath)


            with open(filename, "rU") as input_handle:
                with open(filepath, "w+") as output_handle:
                    sequences = SeqIO.parse(input_handle, "fastq")
                    count = SeqIO.write(sequences, output_handle, "fasta")
        else:
            os.system("mv {filename} AllFastas/")

rule trim:
    input: expand('AllFastas/{sample}.fasta', sample=SAMPLE)

    output: 'Trimmed/{sample}.trim.fasta'

    shell: """ java -jar trimmomatic-0.35.jar SE -phred33
    {input} {output} ILLUMINACLIP:TruSeq3-SE:2:30:10
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 """

rule Diginorm:
    input: expand('Trimmed/{sample}.trim.fasta', sample=SAMPLE)
    output: 'Normalized/{sample}.diginorm.fasta'
    shell:
          "normalize-by-median.py -N 4 -x {input} {output}"
          "filter-abund.py -V normC20k20.kh *.keep {input} {output}"



rule GBdecider:
    input: expand('Normalized/{sample}.diginorm.fasta', sample=SAMPLE)
    output: 'GB.txt'
    run:

        import os
        import sys
        import argparse
        from Bio import SeqIO

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

        filename = str(input)

        filepath = str(output)

        i = get_contig_dict(filenme)

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

rule concatenate:
    input: expand('Normalized/{sample}.diginorm.fasta', sample=SAMPLE)
    output: 'concatenate.txt'
    shell:      # need to change it to add file one at a time
    """
        for file in ~/Files/*.fasta
        do
          echo $file >> fa.txt | sed '$s/ $/\n/'
        done

        cat fa.txt | tr "\n" "," > allFA.txt
        cat allFA.txt | head -c -1 > allFA.txt
        rm fa.txt
        """



rule Trinity:
    input: needed='GB.txt', allfiles='concatenate.txt'
    output: 'Trinity/Trinity.fasta'
    shell: """  Trinity --seqType fa --single  {inputself.allfiles}
    --SS_lib_type F --CPU {input.needed} --output {output} """


#
# after all comparison has done we do sourmash but all above must be in a for loop
# but must run at least twice to have 2 Trinity.fasta to compare
#
#
# rule sourmash:
#     input: file1='Trinity.fasta', file2='Trinity.fasta'
#     output: #compare sequences???
#     run:
#         import time
#         import sys
#
#         try:
#
#         	import screed
#         	import sourmash_lib
#
#         except ImportError:
#         	print("screed and sourmash need to be installed and in the path")
#         	sys.exit(0)
#
#
#
#         def usage():
#         	print("Usage:\npython %s <file A> <file B>  "%(sys.argv[0]))
#         	sys.exit(0)
#
#         def compare_sequences(s1=None , s2=None):
#         	#if imported
#         	if s1 and s2:
#         		seqfile1 = screed.open(s2)
#         		seqfile2 = screed.open(s1)
#         	#if ran standalone
#         	elif (sys.argv[1] and sys.argv[2]):
#         		seqfile1 = screed.open(sys.argv[1])
#         		seqfile2 = screed.open(sys.argv[2])
#         	#in case of some sort of input error
#         	else:
#         		usage()
#
#         	#lists to hold all the sequences in the files
#         	f1 , f2 = [] , []
#
#
#         	t1 = time.time()
#
#         	#load the file into the list
#         	for read in seqfile1:
#         		E = sourmash_lib.Estimators(n=20, ksize=3)
#         		E.add_sequence(read.sequence)
#         		f1.append((read.name , read.sequence , E))
#
#         	t2 = time.time()
#
#         	print("Loaded file %s in %s seconds"%(sys.argv[1] , t2 - t1))
#
#
#         	t1 = time.time()
#
#         	#load the other file into the list
#         	for read in seqfile2:
#         		E = sourmash_lib.Estimators(n=20, ksize=3)
#         		E.add_sequence(read.sequence)
#         		f2.append((read.name , read.sequence , E))
#
#         	t2 = time.time()
#
#         	print("Loaded file %s in %s seconds"%(sys.argv[2] , t2 - t1))
#
#         	#output file
#         	out = open("%s-in-%s.txt"%(sys.argv[1] , sys.argv[2]) , "w")
#
#         	t1 = time.time()
#
#         	#for every sequence in file A
#         	for seq1 in f1:
#
#         		#set start asuming nothing matches
#         		max = (None,None,0)
#
#         		#for every sequence in file B
#         		for seq2 in f2:
#
#         			#get the jaccard index
#         			j = seq2[2].jaccard(seq1[2])
#
#         			#if the index is larger than the previous maximum
#         			if j > max[2]:
#
#         				#it becomes the new max secuence
#         				max = (seq1[0] , seq2[0] , j )
#
#         		#write out the maximum sequence in the output
#         		out.write("\n%s\n%s\ncorrespondence:%s\n"%(max[0] , max[1] , max[2]))
#
#
#         	t2 = time.time()
#
#
#         	print("Finished comparing all the sequences in %s seconds"%(t2 - t1))
#         	print("output:%s-in-%.txts"%(sys.argv[1] , sys.argv[2]))
#
#         if (__name__ == "__main__"):
#         	compare_sequences()
