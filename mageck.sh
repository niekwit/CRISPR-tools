#!/bin/bash

working_dir=$1

mkdir -p "$working_dir/mageck"
cd "$working_dir"

###Performs MAGeCK
#-c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample)
#If pipeline has been called with only pre and post CRISPR library amplification samples, then MAGeCK will not be executed

#determines number of tabs in counts-aggregated.tsv (3 means only two samples in sheet)
sep="\t"
mageck_test=$(awk -F"${sep}" '{print NF-1}' <<< "${test_line}")

if [[ "$test_line" == *"pre"* ]] && [[ "$test_line" == *"post"* ]] && [[ $mageck_test == 3 ]];
	then
		rm -r "$working_dir/mageck"
		echo "No MAGeCK analysis performed"
elif [[ ! -e mageck.config ]];
	then
		rm -r "$working_dir/mageck"
		echo "No MAGeCK analysis performed (mageck.config missing)"
else
	cd "$working_dir/mageck"
	sed '1d' "$working_dir/mageck.config" > "$working_dir/mageck/mageck.config" #removes header from config file
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
