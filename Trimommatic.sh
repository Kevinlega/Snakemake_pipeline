for filename in ./Files/* do
  java -jar trimmomatic-0.35.jar SE -phred33 $filename $filename.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


java -jar ~/ThisSemester/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 Files/ZHour.fastq ZHour_trim.fastq ILLUMINACLIP:ThisSemester/Programs/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:24 SLIDINGWINDOW:4:4 MINLEN:26
