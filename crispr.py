#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:54:01 2021

@author: Niek Wit (University of Cambridge)
"""

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

start = timeit.default_timer()#initiate timing of run

script_dir=os.path.abspath(os.path.dirname(__file__))
work_dir=os.getcwd()

###loads available CRISPR libraries from library.yaml
with open(script_dir+"/library.yaml") as file: library=yaml.full_load(file)

library_list=list(library.keys())

###command line argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-l", "--library", required=True, choices=library_list,
   help="CRISPR library")
ap.add_argument("-t", "--threads", required=False, default=1, metavar="<int>", type=int,
   help="Number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files according to rename.config")
ap.add_argument("-m", "--mismatch", required=False,choices=[0,1], metavar="N",
   help="Number of mismatches (0 or 1) allowed during alignment")
ap.add_argument("-a", "--analysis", required=False, default="mageck", 
                choices=["mageck","bagel2", "ceres"], 
                help="Statistical analysis with MAGeCK or BAGEL2. Default is MAGeCK")
ap.add_argument("-g", "--go", required=False, action='store_true',
   help="GO analysis with DAVID")

args = vars(ap.parse_args())

####set thread count for processing
max_threads=str(multiprocessing.cpu_count())
threads=args["threads"]
if threads == "max":
    threads=max_threads

###run modules based on parsed arguments:
##rename files
rename=args["rename"]
rename_script=script_dir+"/rename.sh"
if rename == True:
    subprocess.run([rename_script,script_dir])

##run FastQC/MultiQC
fastqc_script=script_dir+"/fastqc.sh"
subprocess.run([fastqc_script,str(threads),work_dir])

##count reads
#first determine file extension
file_list=glob.glob("raw-data/*")
test_file=file_list[0]
extension_index=test_file.index(".",0)
file_extension=test_file[extension_index:]

#load settings
crispr_library=args["library"]
mismatch=args["mismatch"]
read_mod=library.get(crispr_library, {}).get('read_mod')
sg_length=library.get(crispr_library, {}).get('sg_length')
index_path=library.get(crispr_library, {}).get('index_path')
clip_seq=library.get(crispr_library, {}).get('clip_seq')
fasta=library.get(crispr_library, {}).get('fasta')

#run count script
count_script=script_dir+"/count.sh"
subprocess.run([count_script,work_dir,str(threads),str(mismatch),
                read_mod,crispr_library,str(sg_length),index_path,
                clip_seq,fasta,script_dir,file_extension])

##run library analysis
library_script=script_dir+"/lib-analysis.sh"
subprocess.run([library_script,work_dir,script_dir])

##run MAGeCK
#get settings
if os.path.exists("config.yml") == True:
        with open("config.yml") as file: config=yaml.full_load(file)
else:
        print("ERROR: config.yml not found. Please provide this file for further analysis.")
        sys.exit()


analysis=args["analysis"]
go=args["go"]
fdr=config.get("mageck-fdr")
go_test=config.get("GO", {}).get('test')
go_term=config.get("GO", {}).get('term')
email=config.get("GO", {}).get('email')
species=config.get("GO", {}).get('species')
if analysis == "mageck":
    print("Statistical analysis with MAGeCK selected")
    mageck_script=script_dir+"/mageck.sh"
    subprocess.run([mageck_script,str(fdr),str(go),go_test,
                go_term,script_dir,email,species,work_dir])
elif analysis == "bagel2":
    print("Statistical analysis with BAGEL2 selected")
    bagel2_dir=config.get("BAGEL2dir")
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


###print total run time
stop = timeit.default_timer()
total_time = stop - start
ty_res = time.gmtime(total_time)
res = time.strftime("%H:%M:%S",ty_res)
print('Total run time: ', res)