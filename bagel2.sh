#!/bin/bash
SCRIPT_DIR=$1
WORKING_DIR=$2
BAGEL2_DIR=$4
fasta=$3

mkdir -p bagel2

#converts count table to BAGEL2 format
Rscript "$SCRIPT_DIR/convert4bagel.R" $WORKING_DIR $fasta 

#performs BAGEL2
sed '1d' "$WORKING_DIR/mageck.config" > "$WORKING_DIR/bagel2/bagel2.config" #removes header from config file
input="$WORKING_DIR/bagel2/bagel2.config"
while IFS= read -r line
do
 	counter=1
 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	bagel2_output="${test_sample}_vs_${control_sample}"
	mkdir -p "bagel2/$bagel2_output"
	
	#to do: create counts-aggregated-bagel with only test and control sample
	head -1 bagel2/counts-aggregated-bagel.txt | tr '\t' '\012' | nl > bagel2/columns.txt #gets columns names and number
	sed -n "${counter}p" bagel2/columns.txt
	
	
	"$BAGEL2_DIR/BAGEL.py" fc -i "bagel2/counts-aggregated-bagel.txt" -o "bagel2/$bagel2_output/$bagel2_output" -c $control_sample
	"$BAGEL2_DIR/BAGEL.py" bf -i "bagel2/$bagel2_output/$bagel2_output.foldchange" -o "bagel2/$bagel2_output/$bagel2_output.bf" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt" -c 1,2,3
	"$BAGEL2_DIR/BAGEL.py" pr -i "bagel2/$bagel2_output/$bagel2_output.bf" -o "bagel2/$bagel2_output/$bagel2_output.pr" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt"
	
	counter=$($number+1 | bc) #counter
done < "$input"
