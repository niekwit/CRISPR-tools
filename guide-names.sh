#!/bin/bash -

#Converts a CRISPR library fasta file to just guide names (CSV) for use in CRISPR pipeline

if (($# != 1));
	then 
		>&2 echo "Error. Usage: $ ./guide-names.sh <library.fasta>"
		exit
fi

file_out=$1
file_out="${file_out%.fasta}.csv"
cat $1 | grep ">" | tr -d ">" > $file_out
