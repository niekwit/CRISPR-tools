# CRISPR screen bioinformatic pipeline


This bioinformatic pipeline will automate analysis of NGS data from CRISPR-Cas9 screen experiments. It uses MAGeCK for statistical analysis.

### Software dependencies:
- Python 3
	- Pandas
	- Numpy
	- Matplotlib
	- Seaborn
- R
- FASTQC
- MultiQC
- Bowtie2 
- Cutadapt
- pigz (if not installed gunzip will be used, but will be slower)
- MAGeCK

### Instructions:

Installation (command line): `https://github.com/niekwit/CRISPR-tools.git`

Create a main folder (can be any name) for the analysis that contains the subfolder `raw-data`, which contains the fastq.gz files.

Navigate to the analysis folder in the command line and type: `path/to/crispr-pipeline.sh [ -l <CRISPR library> ] [ -n <rename.config> OPTIONAL] [-m <INT> number of mismatches allowed for alignment (standard is zero) OPTIONAL]`

CRISPR libraries can be set in the script from line 42, for example:
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



For each library, a Bowtie2 index has to be generated beforehand using `bowtie2 build` with the library fasta file as input. Additionally, a sorted csv file is required that contains each guide name (e.g. A1BG_sgA1BG_1) on a new line.
`index_path` is the file path to the Bowtie2 index. `guides` is the path to the sorted csv file. `read_mod` ("clip" or "trim") sets the method or removing the vector sequence. Use "clip" for libraries with variable guide length (vector sequence to be removed is `clip_seq`) and use "trim" for fixed guide lengths (set guide length with `sg_length`).

Files can be optionally be renamed using the `-n` flag and a config file (`rename.config`), as file names are used to generate sample names for MAGeCK. `-m INT` sets the number of mismatches that are allowed during the alignement step (if not called, zero mismatches are set).
The mageck.config file contains the comparisons that MAGeCK will perform: -c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample).

If you have sequencing samples that represent pre and post CRISPR library amplifications, then these can be named `pre` and `post`, respectively. This will initiate a comparative analysis of these two samples to test for a skew in guide numbers (GINI). When only these two samples are present, no MAGeCK analysis will be performed.
