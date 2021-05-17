library(ceres)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)









### Setup

# Edit this line to point to data directory
data_dir <- "/home/niek/Downloads/example_data"

cn_seg_file <- file.path(data_dir, "/CCLE_copynumber_2013-12-03.seg.txt")
gene_annot_file <- file.path(data_dir, "/CCDS.current.txt")

# Set bowtie index directory. Not needed if $BOWTIE_INDEXES environmental variable is set and includes hg19 index.
Sys.setenv(BOWTIE_INDEXES = file.path(data_dir, "bowtie_indexes"))


gecko_dep_file <- file.path(data_dir, "Gecko.gct")
gecko_rep_map <- file.path(data_dir, "Gecko_replicate_map.tsv")

wang_dep_file <- file.path(data_dir, "Wang2017.gct")
wang_rep_map <- file.path(data_dir, "Wang2017_replicate_map.tsv")



### Run CERES on Gecko data

gecko_inputs_dir <- file.path("./data/gecko_ceres_inputs", Sys.Date())

prepare_ceres_inputs(inputs_dir=gecko_inputs_dir,
                     dep_file=gecko_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=gecko_rep_map,
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad",
                     bowtie_exe = "/home/niek/Downloads/example_data/bowtie-1.3.0-linux-x86_64/bowtie",
                     samtools_exe = "/home/niek/Downloads/example_data/samtools-1.12/bin/samtools")

gecko_ceres <-
  wrap_ceres(sg_path=file.path(gecko_inputs_dir, "guide_sample_dep.Rds"),
             cn_path=file.path(gecko_inputs_dir, "locus_sample_cn.Rds"),
             guide_locus_path=file.path(gecko_inputs_dir, "guide_locus.Rds"),
             locus_gene_path=file.path(gecko_inputs_dir, "locus_gene.Rds"),
             replicate_map_path=file.path(gecko_inputs_dir, "replicate_map.Rds"),
             run_id="Gecko",
             params=list(lambda_g=0.68129207))

gecko_ceres_scaled <-
  scale_to_essentials(gecko_ceres$gene_essentiality_results$ge_fit)
