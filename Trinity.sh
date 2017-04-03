#!/bin/bash
# _file = "allFA.txt"
chmod +r allFA.txt
if [ -s "allFA.txt" ]
then
  Trinity --seqType fa --single ` cat ~/allFA.txt`  --SS_lib_type F --CPU `cat ~/Hello/GBFA.txt` --output ~/Hello
else
  Trinity --seqType fq --single ` cat ~/allFQ.txt`  --SS_lib_type F --CPU `cat ~/Hello/GBFQ.txt` --output ~/Hello
fi



# ./Trinity --seqType fa --single ~/Tnormalize.fasta --SS_lib_type F --max_memory 10G  --output ~/trinitytest
