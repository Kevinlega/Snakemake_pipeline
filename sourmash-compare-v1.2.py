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

