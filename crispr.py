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
    ap.add_argument("-m", "--mismatch", required=False,choices=[0,1], metavar="N",
       help="Number of mismatches (0 or 1) allowed during alignment", default=0)
    ap.add_argument("-a", "--analysis", required=False, default="mageck",
                    choices=["mageck","bagel2"],
                    help="Statistical analysis with MAGeCK or BAGEL2. Default is MAGeCK")
    ap.add_argument("-c", "--cnv", required=False, action='store_true',
       help="Activate CNV correction for MAGeCK/BAGEL2")
    ap.add_argument("-g", "--go", required=False, action='store_true',
       help="GO analysis with DAVID")

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
    utils.fastqc(work_dir,threads,file_extension,exe_dict)

    ##count reads
    #get parsed arguments
    crispr_library=args["library"]
    #check if bowtie2 index is build for CRISPR library
    utils.check_index(library,crispr_library,script_dir,exe_dict,work_dir)
    #check if file with just guide names exists
    utils.guide_names(library,crispr_library)
    #count sgRNAs
    mismatch=args["mismatch"]
    utils.count(library,crispr_library,mismatch,threads,script_dir,work_dir)
    #plot alignment rates
    utils.plot_alignment_rate(work_dir)
    #join count files
    utils.join_counts(work_dir,library,crispr_library)
    #normalise read count table
    utils.normalise(work_dir)
    ##run library analysis
    utils.lib_analysis(work_dir)

    ##run stats on counts
    analysis=args["analysis"]
    go=args["go"]

    if analysis == "mageck":
        print("Statistical analysis with MAGeCK selected")
        cnv=args["cnv"]
        utils.mageck(work_dir,script_dir,cnv)
        #GO analysis
        if go == True:
            utils.go(work_dir,script_dir)
    elif analysis == "bagel2":
        print("Statistical analysis with BAGEL2 selected")
        utils.remove_duplicates(work_dir)
        utils.convert4bagel(work_dir,library,crispr_library)
        utils.bagel2(work_dir,script_dir,exe_dict)


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
