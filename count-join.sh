#!/bin/bash

###Author: Niek Wit (University of Cambridge)###

#Unzip data
pigz -dkv ../raw-data/*fastq.gz

#Fastq quality control
fastqc --threads 30 -o ../fastqc ../raw-data/*fastq.gz
multiqc -o ../fastqc ../fastqc . 2> count.log

#Trims, aligns and counts reads
for file in ../raw-data/*.fastq
do 
	echo $file >> count.log
	file_name=${file##*/} #substring removal of file location
	file2=${file_name%.fastq} #substring removal of .fastq
	extension=".guidecounts.txt"
	output_file=$file2$extension
	cat $file | fastx_trimmer -l 20 -Q33 2>> count.log | bowtie --sam-nohead -5 0 -p40 -t -v0 ~/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUB_only-index - 2>> count.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > $output_file
done

#join count files to generate MAGeCK input file
cp *.guidecounts.txt ../join/
cd ../join
python3 join.py ../join/DUBlibrarysgRNAname.csv
cp counts-aggregated.tsv ../mageck
