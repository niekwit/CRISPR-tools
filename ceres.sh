#!/bin/bash

SCRIPT_DIR=$1
WORKING_DIR=$2
BAGEL2_DIR=$4
fasta=$3
ceres_cn_file=$5

mkdir -p ceres

#create cell line list to check if sample cell line is in CCLE list
if [[ ! -e "$SCRIPT_DIR/CERES/CCLE_cell-line_list.txt" ]]; then
	cat $ceres_cn_file | awk {'print($1)'} | uniq | awk -F "_" '{if ($1 != "CCLE") print($1)}' | sort > "$SCRIPT_DIR/CERES/CCLE_cell-line_list.txt"
fi

#check if sample cell line is in CCLE list

#input="ceres-replicate-list.tsv"
reference=$(cat "$SCRIPT_DIR/CERES/CCLE_cell-line_list.txt" | tr "\n" " ")
test_cell_line=$(cat "$SCRIPT_DIR/config.yml" | shyaml get-value CERES)
#test_cell_line=$(echo $line | awk '{if (NR!=1) {print $2}}' | awk -F "_" '{print($1)}')
if [[ "$reference" != *"$test_cell_line"* ]]; then
	echo "Cell line config.yml not in CCLE list"
	exit 1
fi




###converts count table to BAGEL2 format
Rscript "$SCRIPT_DIR/convert-count.R" $WORKING_DIR $fasta "ceres"

###performs BAGEL2 FC (required for CERES)
#generate array with column names and their column numbers 
head -1 ceres/counts-aggregated-ceres.txt | tr '\t' '\012' | nl | sed 's/\"//g' > ceres/columns.txt #gets columns names and number
declare -A array #array to store column names and numbers
input="ceres/columns.txt"
while IFS= read -r line
do
 	column_name=$(echo $line | awk '{print $(NF)}')
 	column_number=$(echo $line | awk '{print $(NF-1)}')
 	let column_number=column_number-2 #convert to sample number column
	array[$column_name]=$column_number

done < "$input"

sed '1d' "$WORKING_DIR/stats.config" | awk '{if(NF>0) {print $0}}' > "$WORKING_DIR/ceres/ceres.config" #removes header from config file and removes empty line, if any
input="$WORKING_DIR/ceres/ceres.config"
while IFS= read -r line
do
 	test_sample=$(echo "$line" | cut -d ";" -f 1) #splits line of config file into test sample name
	control_sample=$(echo "$line" | cut -d ";" -f 2) #splits line of config file into control sample name
	control_sample_column=${array["$control_sample"]}
	bagel2_output="${test_sample}_vs_${control_sample}"
	mkdir -p "ceres/$bagel2_output"
	echo "Generating fold change table $bagel2_output"
	"$BAGEL2_DIR/BAGEL.py" fc -i "ceres/counts-aggregated-ceres.txt" -o "ceres/$bagel2_output/$bagel2_output" -c $control_sample_column
	echo "Running CERES for $bagel2_output"
	Rscript "$SCRIPT_DIR/convert-fc.R" "ceres/$bagel2_output/$bagel2_output.foldchange" "ceres/$bagel2_output/" $bagel2_output
	#add sgRNA count and sample count to fold-change file
	sgrna_count=$(cat "ceres/$bagel2_output/$bagel2_output.ceres.foldchange" | wc -l)
	let sgrna_count=$sgrna_count-1
	sample_count=$(awk '{print NF}' "ceres/$bagel2_output/$bagel2_output.ceres.foldchange" | sort -nu | tail -n 1)
	let sample_count=$sample_count-2
	cat <(echo -e "$sgrna_count\t$sample_count") "ceres/$bagel2_output/$bagel2_output.ceres.foldchange" > "ceres/$bagel2_output/$bagel2_output.ceres-temp.foldchange"
	cat <(echo -e "#1.2") "ceres/$bagel2_output/$bagel2_output.ceres-temp.foldchange" > "ceres/$bagel2_output/$bagel2_output.ceres.gct"
	rm "ceres/$bagel2_output/$bagel2_output.ceres-temp.foldchange" "ceres/$bagel2_output/$bagel2_output.foldchange" "ceres/$bagel2_output/$bagel2_output.ceres.foldchange"
	#Rscript "$SCRIPT_DIR/ceres.R" "ceres/$bagel2_output/" "ceres/$bagel2_output/$bagel2_output.ceres.gct" $bagel2_output
done < "$input"
