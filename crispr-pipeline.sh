#!/bin/bash

###Author: Niek Wit (University of Cambridge)###

file_path=""
library=""
guides=""

usage() {                                    
  echo "Usage: $0 [ -p /path/to/data ] [ -l CRISPR library ]"
  exit 2
}

while getopts 'p:l:?h' c
do
  case $c in
    p) 
    	file_path=$OPTARG 
    	;;
    l) 
    	library=$OPTARG 
    	;;
    h|?) usage 
    	;;
  esac
done


if [ $library = "bassik" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-bowtie"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-guides-sorted.csv"
		echo "Bassik library selected"
elif [ $library = "moffat_tko1" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat_TKO1-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat-guideslist-sorted.csv"	
		echo "Moffat TKO1 library selected"
elif [ $library = "sabatini" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Sabatini/sabatini-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Sabatini/sabatini-guides-names.csv"	
		echo "Sabatini library selected"
elif [ $library = "dub-only" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUB_only-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUBlibrarysgRNAname.csv"
		echo "DUB only library selected"
fi


cd $file_path
mkdir -p {count,fastqc,mageck}

#Unzip data
echo "Decompressing .fastq.gz files"
pigz -dkv raw-data/*fastq.gz

#Fastq quality control
echo "Performing FastQC"
fastqc --threads 30 -o fastqc/ raw-data/*fastq.gz
echo "Performing MultiQC"
multiqc -o fastqc/ fastqc/ . 2> crispr.log

#Trims, aligns and counts reads
echo "Aligning reads to reference"
for file in raw-data/*.fastq
do 
	echo $file >> crispr.log
	file_name=${file##*/} #substring removal of file location
	file2=${file_name%.fastq} #substring removal of .fastq
	extension=".guidecounts.txt"
	output_file=$file2$extension
	cat $file | fastx_trimmer -l 20 -Q33 2>> crispr.log | bowtie --sam-nohead -5 0 -p40 -t -v0 $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
done

cp /home/niek/Documents/scripts/CRISPR-tools/mageck-join.py count/

#Creates MAGeCK inout file
python3 count/mageck-join.py $guides
mv counts-aggregated.tsv mageck/

#Removes uncompressed .fastq files
rm -r raw-data/*.fastq

#Performs MAGeCK



