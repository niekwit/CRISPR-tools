#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:54:01 2021

@author: Niek Wit (University of Cambridge)
"""

#import os
import sys
import argparse
import subprocess
import multiprocessing
import yaml
import timeit
import time

start = timeit.default_timer()

#Setting variablea for parsing to bash scripts
script_dir=os.path.abspath(os.path.dirname(__file__))


library_list=[]

ap = argparse.ArgumentParser()

ap.add_argument("-l", "--library", required=True, choices=library_list,
   help="CRISPR library")
ap.add_argument("-t", "--threads", required=False,
   help="<INT> number of CPU threads to use")
ap.add_argument("-r", "--rename", required=False, action='store_true',
   help="Rename fq files")
ap.add_argument("-m", "--mismatch", required=False,
   help="<INT> number of mismatches allowed during alignment")
ap.add_argument("-g", "--go", required=False, action='store_true',
   help="GO analysis with DAVID")
args = vars(ap.parse_args())


#set thread count for processing
max_threads=str(multiprocessing.cpu_count())
threads=args["threads"]
if threads == "max":
    threads=max_threads






#print total run time
stop = timeit.default_timer()
total_time = stop - start
ty_res = time.gmtime(total_time)
res = time.strftime("%H:%M:%S",ty_res)
print('Total run time: ', res)