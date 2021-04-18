#!/bin/bash

fdr=$1
go=$2
go_test=$3
go_term=$4
SCRIPT_DIR=$5
email=$6
species=$7

file_list=$(find . -name "*gene_summary.txt") #lists all paths to MAGeCK output files

if [[ -n "$file_list" ]];
	then
		for file in $file_list
		do
			save_path=$(echo $file | sed 's|\(.*\)/.*|\1|') #removes file name from $file
			save_path="${save_path}/"
			Rscript "${SCRIPT_DIR}/plot-hits.R" $fdr $file $species $save_path
			#enter code for GO analysis depletion and enrichment:
			go_folder=$save_path"/GO-analysis-$go_test-$go_term"
			if [[ $go_analysis == "TRUE" ]];
				then
					if [[ ! -d $go_folder ]]
						then
							mkdir -p $go_folder
							Rscript "${SCRIPT_DIR}/go.R" $email $file $fdr $species $go_folder $go_test $go_term
					fi
			fi
		done
fi
