#!/bin/bash

### Author: Niek Wit (University of Cambridge) 2020 ###

file_path=""
library=""
rename_config="NULL"
remove_fq=""
align_mm=0

usage() {                                    
  echo "Usage: $0 [ -p /path/to/data ] [ -l <CRISPR library> ] [ -n <rename.config> OPTIONAL] [-r Removes uncompressed fq files after analysis OPTIONAL][-m # of mismatches allowed for alignment (standard is zero) OPTIONAL]"
  echo -e "CRISPR library options:\nbassik (Morgens et al 2017 Nat. Comm.)\nmoffat_tko1 (Hart et al 2015 Cell)\nsabatini (Park et al 2016 Nat. Gen.)\ndub-only (Nathan lab, unpublished)"
  exit 2
}

while getopts 'rp:l:n:m:?h' c
do
  case $c in
    r)
    	remove_fq="TRUE"
    	;;
    p) 
    	file_path=$OPTARG 
    	;;
    l) 
    	library=$OPTARG 
    	;;
    n)	rename_config=$OPTARG 
    	;;
    m)  align_mm=$OPTARG
    	;;	
    h|?) usage 
    	;;
  esac
done

if [ $library = "bassik" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-bowtie"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-guides-sorted.csv"
		read_mod="clip"
		clip_seq="GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA"
		echo "Bassik library selected"
		
elif [ $library = "moffat_tko1" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat_TKO1-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat-guideslist-sorted.csv"
		read_mod="trim"	
		echo "Moffat TKO1 library selected"
elif [ $library = "sabatini" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Sabatini/sabatini-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Sabatini/sabatini-guides-names.csv"
		read_mod="trim"	
		echo "Sabatini library selected"
elif [ $library = "dub-only" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUB_only-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUBlibrarysgRNAname.csv"
		read_mod="trim"
		echo "DUB only library selected"
fi

start_time=$(date +%s)

cd $file_path
mkdir -p {count,fastqc,mageck}

#renames fastq.gz files
if [ $rename_config != "NULL" ];
	then
		input=$rename_config
		while IFS= read -r line
		do
		  original_file=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into original file name
		  new_file=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into new file name
		  mv "raw-data/${original_file}" "raw-data/${new_file}"
		done < "$input"
	else
		:
fi

#Unzip data
echo "Decompressing .fastq.gz files"
pigz -dkv raw-data/*fastq.gz

#Fastq quality control
echo "Performing FastQC"
fastqc --threads 30 -o fastqc/ raw-data/*fastq.gz
echo "Performing MultiQC"
multiqc -o fastqc/ fastqc/ . 2> crispr.log

#Trims, aligns and counts reads
echo "Aligning reads to reference ($align_mm mismatch(es) allowed)"
if [ $read_mod == "trim" ];
	then
		for file in raw-data/*.fastq
		do 
			echo $file >> crispr.log
			file_name=${file##*/} #substring removal of file location
			file2=${file_name%.fastq} #substring removal of .fastq
			extension=".guidecounts.txt"
			output_file=$file2$extension
			cat $file | fastx_trimmer -l 20 -Q33 2>> crispr.log | bowtie --sam-nohead -5 0 -p40 -t "-v$align_mm" $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
		done
elif [ $read_mod == "clip" ];
	then
		for file in raw-data/*.fastq
		do 
			echo $file >> crispr.log
			file_name=${file##*/} #substring removal of file location
			file2=${file_name%.fastq} #substring removal of .fastq
			extension=".guidecounts.txt"
			output_file=$file2$extension
			cat $file | fastx_clipper -Q33 -l 12 -a $clip_seq -v -n 2>> crispr.log | bowtie --sam-nohead -5 0 -p40 -t "-v$align_mm" $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
		done
fi

cp /home/niek/Documents/scripts/CRISPR-tools/mageck-join.py count/

#Creates MAGeCK input file
cd count/
python3 mageck-join.py $guides

#Normalises MAGeCK input file to total read count
working_dir=$(pwd)
Rscript /home/niek/Documents/scripts/CRISPR-tools/normalise.R $working_dir

cp counts-aggregated.tsv ../mageck/

#Optionally removes uncompressed .fastq files
if [[ $remove_fq == "TRUE" ]];
	then
		rm -r ../raw-data/*.fastq
	else
		:
fi

#Performs MAGeCK
##-c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample)
cd ../mageck
sed '1d' ../mageck.config > mageck.config #removes header from config file
input="mageck.config"
while IFS= read -r line
do
  test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
  control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
  mageck_output="${test_sample}_vs_${control_sample}"
  mkdir $mageck_output
  mageck test -k counts-aggregated.tsv -t $test_sample -c $control_sample --sort-criteria neg -n $mageck_output/$mageck_output 2>> ../crispr.log
done < "$input"

end_time=$(date +%s)
runtime=$((end_time-start_time))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

