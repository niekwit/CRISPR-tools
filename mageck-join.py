#!/usr/bin/env python3
# -*- coding: utf-8 -*-
docstring="""
2019
@author: Niek Wit (MRC-LMB)

This script will join all your individual Bowtie alignment files (.txt) into one file that is suitable for MAGeCK analysis

Instructions:
1. Create a new folder
2. Copy in your sorted library file that only contains gene/sgRNA info (like this: A1BG_sgA1BG_1) in one column
3. Copy in your Bowtie output .txt files (can be any name) that need to be joined for MAGeCK analysis
    - IMPORTANT: do not include any other .txt files
4. In terminal type: python3 join.py <sgRNA_name_list.csv>
5. The output file (counts-aggregated.tsv) can be used for MAGeCK analysis
"""

import glob
import pandas as pd
import sys
import csv

if len(sys.argv) != 2:
    sys.exit('\nJust one argument required (see instructions)%s' %(docstring))
    
sgrnas_list00 = list(csv.reader(open(sys.argv[1])))

sgrnas_list0 = []

for x in sgrnas_list00: #Flattens the list
    for y in x:
        sgrnas_list0.append(y)

#Generates sgRNA and gene columns for final output
sgRNA_output = []
gene_output = []

for n in sgrnas_list0:
    s,g = n.split("_", 1)
    sgRNA_output.append(g)
    gene_output.append(s)

#Generates reference Pandas data frame from sgRNA list library file
d0 = {'sgRNA':pd.Series(sgRNA_output),'gene':pd.Series(gene_output),'sgRNA2':pd.Series(sgrnas_list0)}
dfjoin1 = pd.DataFrame(d0) #sgRNA/gene column required for MAGeCK, sgRNA2 is needed for join operation (deleted later)

#Generates a list of all count .txt files
file_list = glob.glob('*.txt')
file_list.sort()
file_list2 = [w.replace('.txt','') for w in file_list] #this list will generate the column headers for the output file (removes .txt)

#Counts number of .txt files in script folder
txtnumber = len(file_list)

#Generates list of lists for join function output
cycle = 1
master_count_list0 = []
while cycle <= txtnumber:
    master_count_list0.append("count_list"+ str(cycle))
    cycle +=1
master_count_list1 = []
for i in master_count_list0:
    master_count_list1.append([i])
    
cycle = 1
master_sgrna_list0 = []
while cycle <= txtnumber:
    master_sgrna_list0.append("sgrna_list"+ str(cycle))
    cycle +=1
master_sgrna_list1 = []
for i in master_sgrna_list0:
    master_sgrna_list1.append([i])

#Generates Pandas data frame and adds each of the count files in the folder to it after joining
counter = 0
while counter < txtnumber:
    #Opens count files and extract counts and sgRNA names
    file = list(csv.reader(open(file_list [counter])))
    
    for x in file:
        a = str(x)
        if a.count(' ') > 1:
            z,b,c = a.split()
            bint = int(b)
            cmod = c.replace("']","")
            master_count_list1 [counter].append(bint)
            master_sgrna_list1 [counter].append(cmod)
        else:
            b,c = a.split()
            bint = b.replace("['","")
            bint = int(bint)
            cmod = c.replace("']","")
            master_count_list1 [counter].append(bint)
            master_sgrna_list1 [counter].append(cmod)
    
    #Generates Pandas data frame for the data
    d1 = {'sgRNA2':pd.Series(master_sgrna_list1 [counter]),
        file_list2 [counter]:pd.Series(master_count_list1 [counter])}
    df1 = pd.DataFrame(d1)

    #Performs left join to merge Pandas data frames sets:
    dfjoin1 = pd.merge(dfjoin1, df1, on='sgRNA2', how='left')
    dfjoin1 = dfjoin1.fillna(0) #Replaces nan with zero
       
    counter +=1

#Deletes sgRNA2 column from dataframe (only needed for joining, not for MAGeCK)
dfjoin2 = dfjoin1.drop(columns='sgRNA2')
 
#Writes all data to a single .tsv file, ready for MAGeCK
dfjoin2.to_csv('counts-aggregated.tsv', sep='\t',index=False)

