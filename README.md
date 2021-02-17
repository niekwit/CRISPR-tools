# CRISPR screen bioinformatic pipeline

This bioinformatic pipeline will automate analysis of NGS data from CRISPR-Cas9 screen experiments. It uses MAGeCK for statistical analysis.


### Software dependencies:

- Python 3
	- [Pandas](https://pandas.pydata.org/) 
	- [Numpy](https://numpy.org/)
	- [Matplotlib](https://matplotlib.org/stable/index.html)
	- [seaborn](https://seaborn.pydata.org/index.html)
	- [shyaml](https://pypi.org/project/shyaml/)
- [R](https://www.r-project.org/)
	- [Tidyverse](https://www.tidyverse.org/)
	- [ggrepel](https://www.rdocumentation.org/packages/ggrepel/versions/0.9.1)
	- [GO.db](https://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [pigz](https://zlib.net/pigz/) (if not installed gunzip will be used, but will be slower)
- [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/)


### Installation:

Installation from the command line:
> `git clone https://github.com/niekwit/CRISPR-tools.git`

Dependencies can be installed by running:
> `./setup.sh`

The `crispr-pipeline.sh` can be permamently added to $PATH by adding the following line to `~/.bashrc`:
> `export PATH=/home/path/to/CRISPR-tools:$PATH`


### Configuration:

CRISPR libraries can be configured in the `config.yml` file (located in the `CRISPR-tools` folder), as follows:
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
Explanation of `config.yml`:
* For each library a Bowtie2 index has to be generated using `bowtie2 build` with the library fasta file as input. The path to this index has to be set as `index_path` in `config.yaml`. 
* The path to the fasta file must be set with `fasta` (fasta files for a variety of CRISPR libraries can be found in the `Addgene_CRISPR_libraries_FASTA` folder).
* If a CRISPR library has a fixed sgRNA length, then the length of the sgRNA must be set with the `sg_length` variable. Additionaly, set `read_mod` as "trim".
* If a CRISPR library has variable sgRNA lengths, then `read_mod` should be set and "clip" and `clip_seq` should contain the vector sequence downstream of the sgRNA sequence.

Important: when a variable is not used (e.g. `clip_seq` for a fixed sgRNA length CRISPR library), it should be left empty, see example.


### Usage:

1. Create a main folder (can be any name) for the analysis that contains the subfolder `raw-data`, which contains the fastq.gz files.

2. If you want to rename your sequencing files (the files names will be used as sample names for the MAGeCK analsysis so it is recommended to abbreviate them), then this can be set with the `rename.config` file.
On each line put the existing file name and the desired new file name, seperated by a semi-colon (do not include any white space), for example:
```
S25_S2_L001_R1_001.fastq.gz;S25.fastq.gz
L8_S1_L001_R1_001.fastq.gz;L8.fastq.gz
S8_S3_L001_R1_001.fastq.gz;S8.fastq.gz
S15_S4_L001_R1_001.fastq.gz;S15.fastq.gz
```
The `rename.config` file should be placed in the main analysis folder.

If your samples contain sequencing data from amplifications of the CRISPR library itself, then these can be named `pre` and `post`, with `pre` being the pre-amplification DNA (i.e. what was delivered from Addgene), and  with `post` being the post-amplification DNA (i.e. your own maxiprep of the CRISPR library). Renaming these samples in this way, will trigger a comparative analysis of these two samples that will show any skew in sgRNA numbers (depicted by the [GINI index](https://en.wikipedia.org/wiki/Gini_coefficient)). A good library amplification will maintain the same sgRNA number skew as the original prep.

3. If you want to perform a comparative analysis between samples using MAGeCK, then a `mageck.config` file has to be created, for example:
```
t;c
S8;L8
S15;L8
S25;L8
```
c: reference sample, t: test sample. In the MAGeCK output file: neg rank(genes that drop out in test sample)/pos rank(genes that are overrepresented in test sample).

The `mageck.config` file should be placed in the main analysis folder.

4. To start the analysis, navigate to the analysis folder in the command line and type: 
> `path/to/crispr-pipeline.sh -l bassik -r`

This will run the pipeline with the Bassik sgRNA library and will rename the NGS files according to `rename.config`. If `-m` and `-t` are not called, it will allow zero mismatches when aligning the sample to the index, and it will use all the CPU threads that are available on your system when applicable, respectively. 

To allow one mismatch during alignment use the `-m 1` flag.
To set a fixed number of CPU threads use the `-t <INT>` flag.


### Output:

Several folder/files will be generated:

* `fastqc`: contains FastQC and MultiQC analyses on the fastq.gz files.
* `count`: contains the sgRNA counts in individual files, and counts of all files collated in one file (`counts-aggregated.tsv`, MAGeCK input file). It also contains a normalised version of `counts-aggregated.tsv`.
* `library analysis`: contains the analyses of the pre and post library amplification samples.
* `mageck`: contains the MAGeCK output files. It will also contain plots of the results with the top 10 genes marked.

