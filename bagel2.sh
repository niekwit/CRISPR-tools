#!/bin/bash
SCRIPT_DIR=$1
WORKING_DIR=$2
fasta=$3

mkdir -p bagel2
Rscript "$SCRIPT_DIR/convert4bagel.R" $WORKING_DIR $fasta


sed '1d' "$WORKING_DIR/bagel2.config" > "$WORKING_DIR/bagel2/bagel2.config" #removes header from config file
input="$WORKING_DIR/bagel2/bagel2.config"
while IFS= read -r line
do
 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	bagel2_output="${test_sample}_vs_${control_sample}"
	mkdir -p "bagel2/$bagel2_output"
	BAGEL.py fc -i "bagel2/counts-aggregated-bagel.txt" -o "bagel2/$bagel2_output/$bagel2_output" -c $control_sample
	
done < "$input"
