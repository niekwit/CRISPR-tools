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
import timeit
import time

start = timeit.default_timer()#initiate timing of run

script_dir=os.path.abspath(os.path.dirname(__file__))
work_dir=os.getcwd()

#loads available CRISPR libraries from library.yaml
with open(script_dir+"/library.yaml") as file: library=yaml.full_load(file)

library_list=list(library.keys())

#command line argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-l", "--library", required=True, choices=library_list,
   help="CRISPR library")
ap.add_argument("-t", "--threads", required=False, default=1,
   help="<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files according to rename.config")
ap.add_argument("-m", "--mismatch", required=False,choices=[0,1],
   help="<INT> number of mismatches allowed during alignment")
ap.add_argument("-g", "--go", required=False, action='store_true',
   help="GO analysis with DAVID")
args = vars(ap.parse_args())

#set thread count for processing
max_threads=str(multiprocessing.cpu_count())
threads=args["threads"]
if threads == "max":
    threads=max_threads

#check if file extension of raw data is fq or fastq
for fname in os.listdir(work_dir+"/raw-data"):
    if fname.endswith('.fq.gz'):
        extension=".fq.gz"
        break
    elif fname.endswith('.fastq.gz'):
        extension=".fastq.gz"

###run modules based on parsed arguments:
#rename files
rename=args["rename"]
rename_script=script_dir+"/rename.sh"
if rename == True:
    subprocess.run([rename_script,script_dir])

#run FastQC/MultiQC
fastqc_script=script_dir+"/fastqc.sh"
subprocess.run([fastqc_script,str(threads),work_dir])

#count reads
crispr_library=args["library"]
mismatch=args["mismatch"]
read_mod=library.get(crispr_library, {}).get('read_mod')
sg_length=library.get(crispr_library, {}).get('sg_length')
index_path=library.get(crispr_library, {}).get('index_path')
clip_seq=library.get(crispr_library, {}).get('clip_seq')
fasta=library.get(crispr_library, {}).get('fasta')

count_script=script_dir+"/count.sh"
subprocess.run([count_script,work_dir,str(threads),str(mismatch),
                read_mod,crispr_library,str(sg_length),index_path,
                clip_seq,fasta,script_dir])

#run library analysis
library_script=script_dir+"/lib-analysis.sh"
subprocess.run([library_script,work_dir,script_dir])

#run MAGeCK
mageck_script=script_dir+"/mageck.sh"
subprocess.run([mageck_script,work_dir])

#plot top 10 MAGeCK hits and perform GO analysis
with open(script_dir+"/config.yaml") as file: config=yaml.full_load(file)
go=args["go"]
fdr=config.get("mageck-fdr")
plot_script=script_dir+"/plot-hits.sh"
go_test=config.get("GO", {}).get('test')
go_term=config.get("GO", {}).get('term')
email=config.get("GO", {}).get('email')
species=config.get("GO", {}).get('species')

subprocess.run([plot_script,fdr,go,go_test,
                go_term,script_dir,email,species])



#print total run time
stop = timeit.default_timer()
total_time = stop - start
ty_res = time.gmtime(total_time)
res = time.strftime("%H:%M:%S",ty_res)
print('Total run time: ', res)