#!/bin/bash

working_dir=$1
SCRIPT_DIR=$2

lib_folder="library-analysis/"
if [[ ! -d  "$lib_folder" ]];
	then
		test_line=$(head -1 "$working_dir/count/counts-aggregated.tsv")
		if [[ "$test_line" == *"pre"* ]] && [[ "$test_line" == *"post"* ]];
			then
				mkdir -p "$working_dir/library-analysis"
				echo "Performing pre- and post-library amplification comparative analysis"
				cd "$working_dir/count"
				python3 -W ignore "${SCRIPT_DIR}/library-analysis.py"
				#Rscript "${SCRIPT_DIR}/gc-bias.R" "$working_dir/count/" "$fasta"
		fi
fi
