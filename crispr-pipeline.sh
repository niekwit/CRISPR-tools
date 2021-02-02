#!/bin/bash

### Author: Niek Wit (University of Cambridge) 2020 ###

library=""
rename_config="NULL"
remove_fq=""
align_mm=0

SCRIPT_DIR="/home/niek/Documents/scripts/CRISPR-tools/"

usage() {                                    
  echo "Usage: $0 [ -l <CRISPR library> ] [ -r <rename.config> OPTIONAL] [-m <INT> number of mismatches allowed for alignment (standard is zero) OPTIONAL]"
  echo -e "CRISPR library options:\nbassik (Morgens et al 2017 Nat. Comm.)\nmoffat_tko1 (Hart et al 2015 Cell)\nmoffat_tko3 (Hart et al 2017 G3/Mair et al 2019 Cell Rep)\nsabatini (Park et al 2016 Nat. Gen.)\ndub-only (Nathan lab, unpublished)"
  exit 2
}

while getopts 'l:r:m:?h' c
do
  case $c in
    l) 
    	library=$OPTARG 
    	;;
    r)	
    	rename_config=$OPTARG 
    	;;
    m)  
    	align_mm=$OPTARG
    	;;	
    h|?) usage 
    	;;
  esac
done

if [[ $align_mm == 0 ]] || [[ $align_mm == 1 ]];
then
    :
else
    echo "ERROR: -m parameter must be 0 or 1"
    usage
    exit 1
fi

if [ $library = "bassik" ];
	then
		index_path="/home/niek/Documents/references/bowtie2-index/bassik/bassik"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-guides-sorted.csv"
		read_mod="clip"
		clip_seq="GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA"
		species="human"
		echo "Bassik library selected"
elif [ $library = "moffat_tko1" ];
	then
		index_path="/home/niek/Documents/references/bowtie2-index/moffat_tko1/moffat_tko1"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat-guideslist-sorted.csv"
		read_mod="trim"
		sg_length=20
		species="human"
		echo "Moffat TKO1 library selected"
elif [ $library = "moffat_tko3" ];
	then
		index_path="/home/niek/Documents/references/bowtie2-index/moffat_tko3/moffat_tko3"
		guides="/home/niek/Documents/references/bowtie2-index/moffat_tko3/Moffat_TKO3-guide-names.csv"
		read_mod="trim"
		sg_length=20
		species="human"	
		echo "Moffat TKO3 library selected"
elif [ $library = "sabatini" ];
	then
		index_path="/home/niek/Documents/references/bowtie2-index/sabatini/sabatini"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Sabatini/sabatini-guides-names.csv"
		read_mod="trim"
		sg_length=20
		species="human"	
		echo "Sabatini library selected"
elif [ $library = "dub-only" ];
	then
		index_path="/home/niek/Documents/references/bowtie2-index/dub-only/bowtie2-dub-only-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/DUB-only-lib/DUBlibrarysgRNAname.csv"
		read_mod="trim"
		sg_length=20
		species="human"
		echo "DUB only library selected"
fi

start_time=$(date +%s)

mkdir -p {count,fastqc,mageck}

#renames files
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

#checks if file extension is .fq or .fastq
input_format=$(ls raw-data/*.gz | head -1)

if [[ $input_format == *"fastq"* ]]
	then
  		input_extension=".fastq.gz"
elif [[ $input_format == *"fq"* ]]
	then
  		input_extension=".fq.gz"		
fi

#determines CPU thread count
max_threads=$(nproc --all)

#Fastq quality control
echo "Performing FastQC"
fastqc --threads $max_threads -o fastqc/ raw-data/*$input_extension
echo "Performing MultiQC"
multiqc -o fastqc/ fastqc/ . 2> crispr.log

#Trims, aligns and counts reads
echo "Aligning reads to reference ($align_mm mismatch(es) allowed)"
if [ $read_mod == "trim" ];
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
			echo -ne "\r$count/$totalfiles"
			cutadapt -j $max_threads --quality-base 33 -l 20 -o - $file 2>> crispr.log | bowtie2 --no-hd -p $max_threads -t -N $align_mm -x $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
		done
		echo -e \
elif [ $read_mod == "clip" ];
	then
		count=0
		totalfiles=$(ls raw-data/*$input_extension | wc -l)		
		for file in raw-data/$input_extension
		do 
			count=$(($count + 1))			
			echo $file >> crispr.log
			file_name=${file##*/} #substring removal of file location
			file2=${file_name%$input_extension} #substring removal of file extension
			extension=".guidecounts.txt"
			output_file=$file2$extension
			echo -ne "\r$count/$totalfiles"
			cutadapt -j $max_threads --quality-base 33 -a $clip_seq -o - $file 2>> crispr.log | bowtie2 --trim5 1 --no-hd -p $max_threads -t -N $align_mm -x $index_path - 2>> crispr.log | sed '/XS:/d' | cut -f3 | sort | uniq -c > count/$output_file
		done
		echo -e \
fi

#removes first line from guide count text files (bowtie2 artefact)
for file in count/*.guidecounts.txt
do
	sed '1d' $file > "$file.temp" && mv "$file.temp" $file
done

cp "${SCRIPT_DIR}mageck-join.py" count/

#Creates MAGeCK input file
cd count/
python3 mageck-join.py $guides

#Normalises MAGeCK input file to total read count
working_dir=$(pwd)
Rscript "${SCRIPT_DIR}normalise.R" $working_dir

cp counts-aggregated.tsv ../mageck/

#Performs CRISPR maxiprep sequencing analysis (only when samples are present in counts-aggregated.tsv)
##Name pre-amplication sample `pre` and post-amplification sample `post`
test_line=$(head -1 counts-aggregated.tsv)
if [[ "$test_line" == *"pre"* ]] && [[ "$test_line" == *"post"* ]]; 
	then
  		mkdir ../library-analysis
  		echo "Performing pre- and post-library amplification comparative analysis" 
  		python3 -W ignore "${SCRIPT_DIR}library-analysis.py"
fi

###Performs MAGeCK
#-c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample)
#If pipeline has been called with only pre and post CRISPR library amplification samples, then MAGeCK will not be executed 

#determines number of tabs in counts-aggregated.tsv (3 means only two samples in sheet)
sep="\t"
mageck_test=$(awk -F"${sep}" '{print NF-1}' <<< "${test_line}")

if [[ "$test_line" == *"pre"* ]] && [[ "$test_line" == *"post"* ]] && [[ $mageck_test == 3 ]]; 
	then
  		rm -r ../mageck
  		echo "No MAGeCK analysis performed"
	else
		cd ../mageck
		sed '1d' ../mageck.config > mageck.config #removes header from config file
		echo "Performing MAGeCK"
		input="mageck.config"
		while IFS= read -r line
		do
 		 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
  			control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
  			mageck_output="${test_sample}_vs_${control_sample}"
  			mkdir $mageck_output
  			mageck test -k counts-aggregated.tsv -t $test_sample -c $control_sample --sort-criteria neg -n $mageck_output/$mageck_output 2>> ../crispr.log
		done < "$input"	
fi

#plots results and marks top 10 hits, and performs GO analysis of signficant genes (standard FDR cut off is set at 0.25)
fdr=0.25
file_list=$(find . -name "*gene_summary.txt") #lists all paths to MAGeCK output files

for file in $file_list
		do 
			save_path=$(echo $file | sed 's|\(.*\)/.*|\1|') #removes file name from $file 
			save_path="${save_path}/"
			Rscript "${SCRIPT_DIR}plot-hits.R" $fdr $file $species $save_path 
		done

end_time=$(date +%s)
runtime=$((end_time-start_time))
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 )) 
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
