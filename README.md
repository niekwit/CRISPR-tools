# CRISPR screen bioinformatic pipeline


This bioinformatic pipeline will automate analysis of NGS data from CRISPR-Cas9 screen experiments. It uses MAGeCK for statistical analysis.

### Software dependencies:
- Python 3
	- Pandas
	- Numpy
	- Matplotlib
	- Seaborn
	- shyaml
- R
- FASTQC
- MultiQC
- Bowtie2 
- Cutadapt
- pigz (if not installed gunzip will be used, but will be slower)
- MAGeCK

### Installation:

Installation from the command line:
> `git clone https://github.com/niekwit/CRISPR-tools.git`

Dependencies can be installed by running:
> `./setup.sh`

The `crispr-pipeline.sh` can be permamently added to $PATH by adding the following line to `~/.bashrc`:
> `export PATH=/home/path/to/CRISPR-tools:$PATH`



### Usage:

Create a main folder (can be any name) for the analysis that contains the subfolder `raw-data`, which contains the fastq.gz files.

Navigate to the analysis folder in the command line and type: `path/to/crispr-pipeline.sh [ -l <CRISPR library> ] [ -n <rename.config> OPTIONAL] [-m <INT> number of mismatches allowed for alignment (standard is zero) OPTIONAL]`

CRISPR libraries can be configured in the `config.yaml` file, as follows:
```
bassik:
  fasta: ""
  index_path: "/home/niek/Documents/references/bowtie2-index/bassik/bassik"
  read_mod: "clip"
  clip_seq: "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAA"
  sg_length: ""
  species: "human"
moffat_tko1:
  fasta: ""
  index_path: "/home/niek/Documents/references/bowtie2-index/moffat_tko1/moffat_tko1"
  read_mod: "trim"
  clip_seq: ""
  sg_length: 20
  species: "human"
```
For each library, a Bowtie2 index has to be generated beforehand using `bowtie2 build` with the library fasta file as input. 

Fasta files for a variety of CRISPR libraries can be found in the `Addgene_CRISPR_libraries_FASTA` folder.
`index_path` is the file path to the Bowtie2 index.`read_mod` ("clip" or "trim") sets the method or removing the vector sequence. Use "clip" for libraries with variable guide length (vector sequence to be removed is `clip_seq`) and use "trim" for fixed guide lengths (set guide length with `sg_length`). If one variable is not used (for example, when `read_mod: "trim"`, then `sg_length: ""`) the non-relevant variable should be left empty. 

Files can be optionally be renamed using the `-n` flag and a config file (`rename.config`), as file names are used to generate sample names for MAGeCK. `-m INT` sets the number of mismatches that are allowed during the alignement step (if not called, zero mismatches are set).
The mageck.config file contains the comparisons that MAGeCK will perform: -c reference sample, -t test sample: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample).

If you have sequencing samples that represent pre and post CRISPR library amplifications, then these can be named `pre` and `post`, respectively. This will initiate a comparative analysis of these two samples to test for a skew in guide numbers (GINI) and calculates the G/C bias in missing/overreppresented sgRNAs. When only these two samples are present, no MAGeCK analysis will be performed.
