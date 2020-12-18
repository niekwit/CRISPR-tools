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
