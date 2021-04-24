#!/bin/bash
SCRIPT_DIR=$1
WORKING_DIR=$2
BAGEL2_DIR=$4
fasta=$3

mkdir -p bagel2

###converts count table to BAGEL2 format
Rscript "$SCRIPT_DIR/convert4bagel.R" $WORKING_DIR $fasta 

###performs BAGEL2
#generate array with column names and their column numbers (needed for BAGEL.py bf)
#to do: create counts-aggregated-bagel with only test and control sample
head -1 bagel2/counts-aggregated-bagel.txt | tr '\t' '\012' | nl > bagel2/columns.txt #gets columns names and number
sed -n "${counter}p" bagel2/columns.txt #stores column number and column name in file
#split the rows of the files in column number and column name (store as variable) (awk)
declare -A array
input="bagel2/columns.txt"
while IFS= read -r line
do
 	column_name=$(echo $line | awk '{print $(NF)}')
 	column_number=$(echo $line | awk '{print $(NF-1)}')
    array[$column_name]=$column_number

done < "$input"

sed '1d' "$WORKING_DIR/mageck.config" > "$WORKING_DIR/bagel2/bagel2.config" #removes header from config file
input="$WORKING_DIR/bagel2/bagel2.config"
while IFS= read -r line
do
 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	control_sample_column=$(${array["$control_sample"]}-2 | bc)
	test_sample_column=$(${array["$test_sample"]}-2 | bc)
	bagel2_output="${test_sample}_vs_${control_sample}"
	mkdir -p "bagel2/$bagel2_output"
	"$BAGEL2_DIR/BAGEL.py" fc -i "bagel2/counts-aggregated-bagel.txt" -o "bagel2/$bagel2_output/$bagel2_output" -c $control_sample_column
	"$BAGEL2_DIR/BAGEL.py" bf -i "bagel2/$bagel2_output/$bagel2_output.foldchange" -o "bagel2/$bagel2_output/$bagel2_output.bf" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt" -c $test_sample_column
	"$BAGEL2_DIR/BAGEL.py" pr -i "bagel2/$bagel2_output/$bagel2_output.bf" -o "bagel2/$bagel2_output/$bagel2_output.pr" -e "$BAGEL2_DIR/CEGv2.txt" -n "$BAGEL2_DIR/NEGv1.txt"
	
done < "$input"
