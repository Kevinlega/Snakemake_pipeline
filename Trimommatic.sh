for filename in ./Files/* do
  java -jar trimmomatic-0.35.jar SE -phred33 $filename $filename.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
