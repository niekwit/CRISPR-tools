#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import multiprocessing
import yaml
import glob
import urllib.request
import timeit
import time

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
                    choices=["mageck","bagel2", "ceres"],
                    help="Statistical analysis with MAGeCK or BAGEL2. Default is MAGeCK")
    ap.add_argument("-g", "--go", required=False, action='store_true',
       help="GO analysis with DAVID")

    args = vars(ap.parse_args())

    ####set thread count for processing
    threads=utils.set_threads(args)

    ###run modules based on parsed arguments:
    ##rename files
    rename=args["rename"]
    if rename == True:
        utils.rename(work_dir)

    #determine file extension raw data
    file_extension=utils.get_extension(work_dir)

    ##Run FastQC/MultiQC
    utils.fastqc(work_dir,threads,file_extension)

    ##count reads
    #get parsed arguments
    crispr_library=args["library"]
    #check if bowtie2 index is build for CRISPR library
    utils.check_index(library,crispr_library,script_dir)
    #check if file with just guide names exists
    utils.guide_names(library,crispr_library)
    #count sgRNAs
    mismatch=args["mismatch"]
    utils.count(library,crispr_library,mismatch,threads,script_dir,work_dir,file_extension)
    #join count files
    utils.join_counts(work_dir,library,crispr_library)
    #normalise read count table
    utils.normalise(work_dir)
    ##run library analysis
    utils.lib_analysis(work_dir)

    ##run stats on counts
    analysis=args["analysis"]


    if analysis == "mageck":
        print("Statistical analysis with MAGeCK selected")
        utils.mageck(work_dir,script_dir)
    elif analysis == "bagel2":
        print("Statistical analysis with BAGEL2 selected")
        bagel2_dir=config.get("BAGEL2dir")

        utils.convert4bagel(work_dir,library)

        bagel2_script=script_dir+"/bagel2.sh"
        subprocess.run([bagel2_script,script_dir,work_dir,fasta,bagel2_dir])
    elif analysis == "ceres":
        print("Statistical analysis with CERES selected")
        if os.path.isdir(script_dir+'/CERES') == False:
            print("ERROR: no CCLE copy number file present")
            url="https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_2013-12-03.seg.txt"
            download=input("Download CCLE copy numer file from https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/ y/n?")
            if download in ["yes","Yes","y","Y"]:
                os.mkdir(script_dir+"/CERES")
                urllib.request.urlretrieve(url, script_dir+"/CERES/CCLE_copynumber_2013-12-03.seg.txt")
                print("Download finished")
            else:
                sys.exit()
        ceres_cn=glob.glob(script_dir+"/CERES/*.txt")
        if len(ceres_cn) > 1:
            print("ERROR: more than one CCLE copy number file present. Keep only one.")
            sys.exit()
        elif len(ceres_cn) == 0:
            print("ERROR: no CCLE copy number file present.")
            url="https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_2013-12-03.seg.txt"
            download=input("Download CCLE copy numer file from https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/ y/n?")
            if download in ["yes","Yes","y","Y"]:
                urllib.request.urlretrieve(url, script_dir+"/CERES/CCLE_copynumber_2013-12-03.seg.txt")
                print("Download finished")
            else:
                sys.exit()
        print("Running CERES")
        ceres_cn_file=script_dir+"/CERES/CCLE_copynumber_2013-12-03.seg.txt"
        ceres_script=script_dir+"/ceres.sh"
        bagel2_dir=config.get("BAGEL2dir")
        subprocess.run([ceres_script,script_dir,work_dir,fasta,bagel2_dir,ceres_cn_file])

    #GO analysis
    go=args["go"]

if __name__ == "__main__":
    #start run timer
    start = timeit.default_timer()

    script_dir=os.path.abspath(os.path.dirname(__file__))
    work_dir=os.getcwd()

    #adds script directory to runtime for importing modules
    sys.path.append(script_dir)
    import utils

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
