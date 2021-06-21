#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import multiprocessing
import yaml
import glob
import timeit
import time
import pickle

def main():
    ###command line argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-l", "--library", required=True, choices=library_list,
       help="CRISPR library")
    ap.add_argument("-t", "--threads", required=False, default=1, metavar="<int>",
       help="Number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    ap.add_argument("-r", "--rename", required=False, action='store_true',
       help="Rename fastq files according to rename.config")
    ap.add_argument("-m", "--mismatch", required=False,choices=("0","1"), metavar="N",
       help="Number of mismatches (0 or 1) allowed during alignment", default=0)
    ap.add_argument("-a", "--analysis", required=False, default="mageck",
                    choices=["mageck","bagel2"],
                    help="Statistical analysis with MAGeCK or BAGEL2. Default is MAGeCK")
    ap.add_argument("-f","--fdr", required=False, metavar="<FDR value>", default=0.25,
       help="Set FDR cut off for MAGeCK hits (default is 0.25)")
    ap.add_argument("-c","--cnv", required=False, metavar="<CCLE cell line>",default=None,
       help="Activate CNV correction for MAGeCK/BAGEL2 with given cell line")
    ap.add_argument("--go", required=False, action='store_true', default=None,
       help="Gene set enrichment analysis with enrichR")
    ap.add_argument("--skipfastqc", required=False, action='store_true', default=False,
          help="Skip FastQC/MultiQC analysis")

    args = vars(ap.parse_args())

    ###check if software requirements are met
    try:
        exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
        dep_list=("fastqc","bowtie2","mageck","bagel2")
        for i in dep_list:
            if not i in exe_dict:
                sys.exit("ERROR: %s directory not found\nRun setup.py again\n" % i)
    except FileNotFoundError:
        sys.exit("setup.py not run before first analysis")

    ###set thread count for processing
    threads=utils.set_threads(args)

    ###run modules based on parsed arguments:
    ##rename files
    rename=args["rename"]
    if rename == True:
        utils.rename(work_dir)

    #determine file extension raw data
    file_extension=utils.get_extension(work_dir)

    ##Run FastQC/MultiQC
    skip_fastqc=args["skipfastqc"]
    if not skip_fastqc:
        utils.fastqc(work_dir,threads,file_extension,exe_dict)
    else:
        print("Skipping FastQC/MultiQC analysis")

    ##count reads
    #get parsed arguments
    crispr_library=args["library"]
    #check if bowtie2 index is build for CRISPR library
    utils.check_index(library,crispr_library,script_dir,exe_dict,work_dir)
    #check if file with just guide names exists
    utils.guide_names(library,crispr_library)
    #count sgRNAs
    mismatch=args["mismatch"]
    utils.count(library,crispr_library,mismatch,threads,script_dir,work_dir,exe_dict)
    #plot alignment rates
    utils.plot_alignment_rate(work_dir)
    #join count files
    if not utils.file_exists(os.path.join(work_dir,"count",'counts-aggregated.tsv')):
        utils.join_counts(work_dir,library,crispr_library)
    #normalise read count table
    if not utils.file_exists(os.path.join(work_dir,"count","counts-aggregated-normalised.csv")):
        utils.normalise(work_dir)
    ##run library analysis
    utils.lib_analysis(work_dir,library,crispr_library,script_dir)

    ##run stats on counts
    analysis=args["analysis"]
    go=args["go"]

    if analysis == "mageck":
        #set FDR
        fdr=float(args["fdr"])
        if fdr > 0.25:
            print("WARNING: MAGeCK FDR cut off set higher than default 0.25")

        cnv=args["cnv"]
        if not cnv:
            utils.mageck(work_dir,script_dir,cnv,fdr)

    elif analysis == "bagel2":
        print("Running BAGEL2")
        utils.remove_duplicates(work_dir)
        utils.convert4bagel(work_dir,library,crispr_library)
        utils.bagel2(work_dir,script_dir,exe_dict)


    if go == True:
        utils.go(work_dir,script_dir,analysis,fdr)

if __name__ == "__main__":
    #start run timer
    start = timeit.default_timer()

    script_dir=os.path.abspath(os.path.dirname(__file__))
    work_dir=os.getcwd()

    #adds script directory to runtime for importing modules
    sys.path.append(script_dir)
    import crispr_utils as utils

    ###loads available CRISPR libraries from library.yaml
    with open(os.path.join(script_dir,"library.yaml")) as file:
        library=yaml.full_load(file)
    library_list=list(library.keys())

    main()

    #print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print('Total run time: ', res)
