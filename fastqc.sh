#!/bin/bash

threads=$1
file_extension=$2

if [[ ! -d fastqc ]];then
	echo "Performing FastQC/MultiQC"
	mkdir -p fastqc
	fastqc --threads $threads -o fastqc/ "raw-data/*$file_extension" 2>> fastqc.log
	multiqc -o "fastqc/" "fastqc/" . 2>> fastqc.log
elif [[ ! -n "$(ls -A fastqc 2>/dev/null)" ]];then
	echo "Performing FastQC/MultiQC"
	fastqc --threads $threads -o fastqc/ "raw-data/*$file_extension" 2>> fastqc.log
	multiqc -o "fastqc/" "fastqc/" . 2>> fastqc.log
elif [[ -n "$(ls -A fastqc 2>/dev/null)" ]];then
	echo "FastQC/MultiQC already performed"
fi
