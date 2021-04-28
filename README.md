# CRISPR screen bioinformatic pipeline

This bioinformatic pipeline will automate analysis of NGS data from CRISPR-Cas9 screen experiments. It uses MAGeCK for statistical analysis.


## Index

1. [Software dependencies](https://github.com/niekwit/CRISPR-tools#software-dependencies)
2. [Installation](https://github.com/niekwit/CRISPR-tools#installation)
3. [Configuration](https://github.com/niekwit/CRISPR-tools#configuration)
4. [Usage](https://github.com/niekwit/CRISPR-tools#usage)
5. [Output](https://github.com/niekwit/CRISPR-tools#output)

## Software dependencies:

- [Python 3](https://www.python.org/)
	- [Pandas](https://pandas.pydata.org/) 
	- [Numpy](https://numpy.org/)
	- [Matplotlib](https://matplotlib.org/stable/index.html)
	- [seaborn](https://seaborn.pydata.org/index.html)
	- [PyYAML](https://pyyaml.org/)
- [R](https://www.r-project.org/)
	- [Tidyverse](https://www.tidyverse.org/)
	- [ggrepel](https://www.rdocumentation.org/packages/ggrepel/versions/0.9.1)
	- [rJava](https://cran.r-project.org/web/packages/rJava/index.html)
	- [RDAVIDWebService](https://www.bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html)
	- [bioMart](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
	- [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8)
	- [GOplot](https://wencke.github.io/)
	- [stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0)
	- [gridExtra](https://cran.r-project.org/web/packages/gridExtra/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [pigz](https://zlib.net/pigz/) (if not installed gunzip will be used, but will be slower)
- [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/)

The Python3 and R dependencies can be installed with the `setup.sh` script. Other dependencies have to be installed manually. Non-Python3/R dependencies should also be set in your `$PATH`,
for example, to add the Bowtie2 binary to your `$PATH`, add the following line to `~/.bashrc`:
> `export PATH=/path/to/bowtie2-2.4.2-linux-x86_64:$PATH`


## Installation:

Installation from the command line:
> `git clone https://github.com/niekwit/CRISPR-tools.git`

Dependencies can be installed by running (can also be done manually):
> `./setup.sh`

The `crispr-pipeline.sh` can be permamently added to $PATH by adding the following line to `~/.bashrc`:
> `export PATH=/home/path/to/CRISPR-tools:$PATH`

OPTIONAL: to enable auto-completion of the command line options for the CRISPR library (-l) and statistical analyses (-a):
Add this line to your `~/.bashrc` file:
> `source /path/to/CRISPR-tools/auto-complete.sh`

## Configuration:

CRISPR libraries can be configured in the `library.yaml` file (located in the `CRISPR-tools` folder), as follows:
```
bassik:
  fasta: "/home/niek/Documents/references/fasta/Human/Bassik-library/bassik_lib.fasta"
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
Explanation of `library.yaml`:
* For each library a Bowtie2 index has to be generated using `bowtie2 build` with the library fasta file as input. The path to this index has to be set as `index_path` in `library.yaml`. 
* The path to the fasta file must be set with `fasta` (fasta files for a variety of CRISPR libraries can be found in the `Addgene_CRISPR_libraries_FASTA` folder).
* If a CRISPR library has a fixed sgRNA length, then the length of the sgRNA must be set with the `sg_length` variable. Additionaly, set `read_mod` as "trim".
* If a CRISPR library has variable sgRNA lengths, then `read_mod` should be set and "clip" and `clip_seq` should contain the vector sequence downstream of the sgRNA sequence.

Important: when a variable is not used (e.g. `clip_seq` for a fixed sgRNA length CRISPR library), it should be left empty, see example.

For each experiment, a settings file (`config.yml`) has to be provided in the main experimental folder. This file contains the following:
```
mageck-fdr: 0.25
GO:
  email: "your@email.com" #create an account first on https://david.ncifcrf.gov/webservice/register.htm
  species: "human"
  test: "Benjamini" #choose Bonferroni or Benjamini
  term: "BP" #BP, CC or MF
BAGEL2dir: "/home/niek/Documents/scripts/bagel"
```
The `mageck-fdr` variable sets the statistical cut off for plotting MAGeCK hits. 
## Usage:

1. Create a main folder (can be any name) for the analysis that contains the subfolder `raw-data`, which contains the fastq.gz files.

2. If you want to rename your sequencing files (the files names will be used as sample names for the MAGeCK analsysis so it is recommended to abbreviate them), then this can be set with the `rename.config` file that should be located in the main analysis folder.
On each line put the existing file name and the desired new file name, seperated by a semi-colon (do not include any white space), for example:
```
S25_S2_L001_R1_001.fastq.gz;S25.fastq.gz
L8_S1_L001_R1_001.fastq.gz;L8.fastq.gz
S8_S3_L001_R1_001.fastq.gz;S8.fastq.gz
S15_S4_L001_R1_001.fastq.gz;S15.fastq.gz
```
The `rename.config` file should be placed in the main analysis folder.

If your samples contain sequencing data from amplifications of the CRISPR library itself, then these can be named `pre` and `post`, with `pre` being the pre-amplification DNA (i.e. what was delivered from Addgene), and  with `post` being the post-amplification DNA (i.e. your own maxiprep of the CRISPR library). Renaming these samples in this way, will trigger a comparative analysis of these two samples that will show any skew in sgRNA numbers (depicted by the [GINI index](https://en.wikipedia.org/wiki/Gini_coefficient)). A good library amplification will maintain the same sgRNA number skew as the original prep.

3. If you want to perform a comparative analysis between samples using MAGeCK or BAGEL2, then a `stats.config` file has to be created, for example:
```
t;c
S8;L8
S15;L8
S25;L8
```
c: reference sample, t: test sample. In the MAGeCK output file: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample).

The `stats.config` file should be placed in the main analysis folder.

4. To get an overview of all the options for the CRISPR analysis, type `path/to/crispr.py -h, --help` in the command line: 
```
path/to/crispr.py [-h] -l {CRISPR library} [-t N] [-r] [-m N] [-a {mageck,bagel2}] [-g]

optional arguments:
  -h, --help                                      show this help message and exit
  -l {CRISPR library}, --library {CRISPR library} CRISPR library
  -t N, --threads N                               Number (N) of CPU threads to use (default is 1). Use max to apply all
                                                  available CPU threads
  -r, --rename                                    Rename fq files according to rename.config
  -m N, --mismatch N                              Number of mismatches (0 or 1) allowed during alignment
  -a {mageck,bagel2}, --analysis {mageck,bagel2}  Statistical analysis with MAGeCK or BAGEL2. Default is MAGeCK
  -g, --go                                        GO analysis with DAVID
```
To start an analysis, with for example the Bassik library, navigate to main analysis folder in the command line and run:

> `path/to/crispr.py -l bassik -r -t max`

This initiates a run that will rename your samples according to `rename.config`, allows no mismatches during alignment, will use all available CPU threads for the analysis, and uses MAGeCK for statistical analysis
If you also want to use BAGEL2 for statistical analysis simply run afterwards:
> `path/to/crispr.py -l bassik -a bagel2`
This will only run BAGEL2 and skip all steps that are common with MAGeCK.

## Output:

Several folder/files will be generated:

* `fastqc`: contains FastQC and MultiQC analyses on the fastq.gz files.
* `count`: contains the sgRNA counts in individual files, and counts of all files collated in one file (`counts-aggregated.tsv`, MAGeCK input file). It also contains a normalised version of `counts-aggregated.tsv`.
* `library analysis`: contains the analyses of the pre and post library amplification samples.
* `mageck`: contains the MAGeCK output files. It will also contain plots of the results with the top 10 genes marked.
* `bagel2`: contains the BAGEL2 output files. It will also contain plots of the results with the top 10 genes marked.

