# CRISPR screen bioinformatic tools


### Software Dependencies:
- python3 
- R
- FASTQC
- Bowtie 
- MultiQC
- FASTX-Toolkit
- pigz

### Instructions:

Usage: ./crispr-pipeline.sh [ -p /path/to/data ] [ -l <CRISPR library> ] [ -n <rename.config> OPTIONAL] [-r Removes uncompressed fq files after analysis OPTIONAL][-m # of mismatches allowed for alignment (standard is zero) OPTIONAL]

CRISPR libraries can be set in the script from line 38, for example:
```
if [ $library = "bassik" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-bowtie"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Bassik/bowtie-lib/bassik-guides-sorted.csv"
		read_mod="clip"
		clip_seq="GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA"
		echo "Bassik library selected"
```
or
```
elif [ $library = "moffat_tko1" ];
	then
		index_path="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat_TKO1-index"
		guides="/home/niek/Documents/references/bowtie-index/CRISPR/Moffat_TKO1/moffat-guideslist-sorted.csv"
		read_mod="trim"
		sg_length=20	
		echo "Moffat TKO1 library selected"
```




`index_path` is the file path to the Bowtie index. `guides` is the path to a sorted csv file that contains each guide name (e.g. A1BG_sgA1BG_1) on a new line. `read_mod` ("clip" or "trim") sets the method or removing the vector sequence. Use "clip" for libraries with variable guide length (vector sequence to be removed is `clip_seq`) and use "trim" for fixed guide lengths (set with `sg_length`).

