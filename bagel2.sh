#!/bin/bash
SCRIPT_DIR=$1
WORKING_DIR=$2
BAGEL2_DIR=$4
fasta=$3

mkdir -p bagel2

###converts count table to BAGEL2 format
Rscript "$SCRIPT_DIR/convert-count.R" $WORKING_DIR $fasta "bagel2"

###performs BAGEL2
head -1 bagel2/counts-aggregated-bagel2.txt | tr '\t' '\012' | nl | sed 's/\"//g' > bagel2/columns.txt #gets columns names and number
declare -A array #array to store column names and numbers
input="bagel2/columns.txt"
while IFS= read -r line
do
 	column_name=$(echo $line | awk '{print $(NF)}')
 	column_number=$(echo $line | awk '{print $(NF-1)}')
 	let column_number=column_number-2 #convert to sample number column
	array[$column_name]=$column_number

done < "$input"

sed '1d' "$WORKING_DIR/stats.config" | awk '{if(NF>0) {print $0}}' > "$WORKING_DIR/bagel2/bagel2.config" #removes header from config file and removes empty line, if any
input="$WORKING_DIR/bagel2/bagel2.config"
while IFS= read -r line
do
 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	control_sample_column=${array["$control_sample"]}
	bagel2_output="${test_sample}_vs_${control_sample}"
	mkdir -p "bagel2/$bagel2_output"
	echo "Generating fold change table $bagel2_output"
	"$BAGEL2_DIR/BAGEL.py" fc -i "bagel2/counts-aggregated-bagel2.txt" -o "bagel2/$bagel2_output/$bagel2_output" -c $control_sample_column
done < "$input"

head -1 "bagel2/$bagel2_output/$bagel2_output.foldchange" | tr '\t' '\012' | nl | sed 's/\"//g' > bagel2/columns-foldchange.txt #gets columns names and number
declare -A array2 #array to store column names and numbers
input="bagel2/columns-foldchange.txt"
while IFS= read -r line
do
 	column_name=$(echo $line | awk '{print $(NF)}')
 	column_number=$(echo $line | awk '{print $(NF-1)}')
 	let column_number=column_number-2 #convert to sample number column
	array2[$column_name]=$column_number
done < "$input"

input="$WORKING_DIR/bagel2/bagel2.config"
while IFS= read -r line
do
	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	test_sample_column=${array2["$test_sample"]}
	bagel2_output="${test_sample}_vs_${control_sample}"
	echo "Calculating Bayes Factors $bagel2_output"
	"$BAGEL2_DIR/BAGEL.py" bf -i "bagel2/$bagel2_output/$bagel2_output.foldchange" -o "bagel2/$bagel2_output/$bagel2_output.bf" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt" -c $test_sample_column
	echo "Calculating precision-recall $bagel2_output"
	"$BAGEL2_DIR/BAGEL.py" pr -i "bagel2/$bagel2_output/$bagel2_output.bf" -o "bagel2/$bagel2_output/$bagel2_output.pr" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt"
	Rscript "${SCRIPT_DIR}/plot-hits-bagel2.R" "bagel2/$bagel2_output/$bagel2_output.pr" "bagel2/$bagel2_output/" $bagel2_output
done < "$input"
