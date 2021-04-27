#!/bin/bash

working_dir=$1
threads=$2
align_mm=$3
read_mod=$4
library=$5
sg_length=$6
index_path=$7
clip_seq=$8
fasta=$9
SCRIPT_DIR=${10}
input_extension=${11}

#Trims, aligns and counts reads
count_folder="count/"
if [[ ! -d  "$count_folder" ]];
	then
		echo "$library library selected"
		echo "Aligning reads to reference ($align_mm mismatch(es) allowed)"
		mkdir count
		if [ $read_mod == "trim" ];
			then
				count=0
				totalfiles=$(ls raw-data/* | wc -l)
				for file in raw-data/*$input_extension
				do
					count=$(($count + 1))
					echo $file >> crispr.log
					file_name=${file##*/} #substring removal of file location
					file2=${file_name%$input_extension} #substring removal of file extension
					extension=".guidecounts.txt"
					output_file=$file2$extension
					echo -ne "\rAligning file $count/$totalfiles"
					cutadapt -j $threads --quality-base 33 -l "$sg_length" -o - $file 2>> crispr.log | bowtie2 --no-hd -p $threads -t -N $align_mm -x $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
				done
				echo -e "\n"

		elif [ $read_mod == "clip" ];
			then
				count=0
				totalfiles=$(ls raw-data/*$input_extension | wc -l)
				for file in raw-data/*$input_extension
				do
					count=$(($count + 1))
					echo $file >> crispr.log
					file_name=${file##*/} #substring removal of file location
					file2=${file_name%$input_extension} #substring removal of file extension
					extension=".guidecounts.txt"
					output_file=$file2$extension
					echo -ne "\rAligning file $count/$totalfiles"
					cutadapt -j $threads --quality-base 33 -a $clip_seq -o - $file 2>> crispr.log | bowtie2 --trim5 1 --no-hd -p $threads -t -N $align_mm -x $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
				done
				echo -e "\n"
		fi

		#removes first line from guide count text files (bowtie2 artefact)
		for file in count/*.guidecounts.txt
		do
			sed '1d' $file > "$file.temp" && mv "$file.temp" $file
		done

		#generates csv file with only sorted guide names
		guides="${fasta%.fasta}-guide-names.csv"
		if [[ ! -e "$guides" ]];
		then
			cat "$fasta" | grep ">" | tr -d ">" | sort > "$guides"
		fi

		cp "${SCRIPT_DIR}/mageck-join.py" count/
		#Creates MAGeCK input file
		cd count/
		python3 mageck-join.py $guides
		#Normalises MAGeCK input file to total read count
		Rscript "${SCRIPT_DIR}/normalise.R" "$working_dir/count/"
		#Removes duplicate sgRNAs, needed for BAGEL2_DIR
		Rscript "${SCRIPT_DIR}/remove-dups.R" "$working_dir/count/"
else
	echo "Trimming, alignment and counting already performed"

fi
