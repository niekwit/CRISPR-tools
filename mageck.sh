#!/bin/bash

fdr=$1
go_analysis=$2
go_test=$3
go_term=$4
SCRIPT_DIR=$5
email=$6
species=$7
working_dir=$8

mkdir -p "$working_dir/mageck"

###Performs MAGeCK
#-c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample)
#If pipeline has been called with only pre and post CRISPR library amplification samples, then MAGeCK will not be executed

#determines number of tabs in counts-aggregated.tsv (3 means only two samples in sheet)
sep="\t"
mageck_test=$(awk -F"${sep}" '{print NF-1}' <<< "${test_line}")

if [[ "$test_line" == *"pre"* ]] && [[ "$test_line" == *"post"* ]] && [[ $mageck_test == 3 ]];
	then
		rm -r "$working_dir/mageck"
		echo "No MAGeCK analysis performed: no non-library samples present"
elif [[ ! -e stats.config ]];
	then
		rm -r "$working_dir/mageck"
		echo "ERROR: No MAGeCK analysis performed: stats.config missing"
else
	cd "$working_dir/mageck"
	sed '1d' "$working_dir/stats.config" | awk '{if(NF>0) {print $0}}' > "$working_dir/mageck/mageck.config" #removes header from config file and removes empty line, if any
	
	echo "Performing MAGeCK"
	input="$working_dir/mageck/mageck.config"
	while IFS= read -r line
	do
	 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
		control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
		mageck_output="${test_sample}_vs_${control_sample}"
		if [[ ! -d $mageck_output ]];
			then
				mkdir -p $mageck_output
				mageck test -k "$working_dir/count/counts-aggregated.tsv" -t $test_sample -c $control_sample --sort-criteria neg -n $mageck_output/$mageck_output 2>> ../crispr.log
		fi
	done < "$input"
fi

#plot MAGeCK hits
file_list=$(find . -name "*gene_summary.txt") #lists all paths to MAGeCK output files

if [[ -n "$file_list" ]];
	then
		for file in $file_list
		do
			save_path=$(echo $file | sed 's|\(.*\)/.*|\1|') #removes file name from $file
			save_path="${save_path}/"
			Rscript "${SCRIPT_DIR}/plot-hits-mageck.R" $fdr $file $save_path
			go_folder=$save_path"/GO-analysis-$go_test-$go_term"
			if [[ $go_analysis == "True" ]];
				then
					if [[ ! -d $go_folder ]]
						then
							mkdir -p $go_folder
							Rscript "${SCRIPT_DIR}/go.R" $email $file $fdr $species $go_folder $go_test $go_term
					fi
			fi
		done
fi
