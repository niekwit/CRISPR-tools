#!/usr/bin/env python3

import pkg_resources
import os
import subprocess
import yaml
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set(style="whitegrid")

def install_packages(): #check for required python packages; installs if absent
    required = {"pyyaml","pandas","numpy","matplotlib","seaborn"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        subprocess.check_call([python, '-m', 'pip3', 'install', *missing], stdout=subprocess.DEVNULL)

def set_threads():
    max_threads=str(multiprocessing.cpu_count())
    threads=args["threads"]
    if threads == "max":
        threads=max_threads
    return threads

def rename(work_dir):
    file=open(os.path.join(work_dir,"rename.config"), "r")
    lines=file.readlines()
    count=0
    for line in lines: #removes newline characters
        lines[count]=line.replace("\n","")
        count+=1

    for line in lines:#rename files
        old_name,new_name=line.split(";")
        os.rename(os.path.join(work_dir,"raw-data",old_name),os.path.join(work_dir,"raw-data",new_name))

def get_extension(work_dir):
    file_list=glob.glob(os.path.join(work_dir,"raw-data","*"))
    test_file=file_list[0]
    extension_index=test_file.index(".",0)
    file_extension=test_file[extension_index:]
    return file_extension

def file_exists(file):
    if os.path.exists(file):
        print("Skipping "+file+" (already exists/analysed)")
        return(True)
    else:
        return(False)

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep=" ",file=file)

def fastqc(work_dir,threads):
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(work_dir+"/fastqc",exist_ok=True)
        fastqc_command="fastqc --threads "+str(threads)+" --quiet -o fastqc/ raw-data/*.fastq.gz"
        multiqc_command=["multiqc","-o","fastqc/","fastqc/"]
        #log commands
        with open(os.path.join(work_dir,"commands.log"),"w") as file:
            file.write("FastQC: ")
            print(fastqc_command, file=file)
            file.write("MultiQC: ")
            print(*multiqc_command, sep=" ", file=file)
        print("Running FastQC on raw data")
        subprocess.run(fastqc_command, shell=True)
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")

def count():
    pass

def normalise():
    df=pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv"))
    column_range=range(2,len(df.columns))
    for i in column_range:
        column_sum=df.iloc[:,i].sum()
        df.iloc[:,i]=df.iloc[:,i] / column_sum * 1E8
        df.iloc[:,i]=df.iloc[:,i].astype(int)
    df.to_csv(os.path.join(work_dir,"count","counts-aggregated-normalised.csv"),index=False,header=True)

def remove_duplicates(work_dir):
    df=pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv"))
    df=pd.DataFrame.drop_duplicates(df)
    df.to_csv(os.path.join(work_dir,"count","counts-aggregated.tsv"),index=False,header=True,sep="\t")

def mageck():
    pass

def bagel2():
    pass

def ceres():
    pass

def lib_analysis():
    import warnings
    warnings.filterwarnings("ignore")

    df = pd.read_csv('counts-aggregated.tsv',sep='\t')

    #determines total read count per column
    pre_lib_sum = df['pre'].sum()
    pre_lib_sum = int(pre_lib_sum)
    post_lib_sum = df['post'].sum()
    post_lib_sum = int(post_lib_sum)

    #normalises guide counts to total guide count
    #X: pre-amplification library, Y: post-amplification library
    datax = df['pre']
    datax = datax.sort_values(ascending=True)
    X = datax.to_numpy()
    X = X / pre_lib_sum

    datay = df['post']
    datay = datay.sort_values(ascending=True)
    Y = datay.to_numpy()
    Y = Y / pre_lib_sum

    index_len = len(df.index)

    #plots data:
    ax = sns.lineplot(x=range(index_len),y=X, color='navy',label='Pre-amplification library')
    ax = sns.lineplot(x=range(index_len),y=Y, color='green',label='Post-amplification library')
    ax.set_yscale('log')
    ax.legend(loc='lower right')
    ax.set(ylabel='Normalised sgRNA count', xlabel='sgRNA')
    plt.savefig('../library-analysis/normalised-guides-frequency.pdf')
    plt.close()

    #
    datax2 = df['pre']
    X2 = datax2.to_numpy()
    X2 = X2 / pre_lib_sum

    datay2 = df['post']
    Y2 = datay2.to_numpy()
    Y2 = Y2 / pre_lib_sum

    data2 = X2 / Y2
    data2 = np.sort(data2)

    ax = sns.lineplot(x=range(index_len),y=data2, color='navy')
    ax.set(ylabel='Normalised \n pre-amplification/post-amplification', xlabel='sgRNA')
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('../library-analysis/normalised-pre-amplification-post-amplification.pdf')
    plt.close()

    #Calculates Gini index of data sets
    ##Code taken and adapted from https://zhiyzuo.github.io/Plot-Lorenz/

    #function to calculate Gini index
    def gini(arr):
        ## first sort
        sorted_arr = arr.copy()
        sorted_arr.sort()
        n = arr.size
        coef_ = 2. / n
        const_ = (n + 1.) / n
        weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
        return coef_*weighted_sum/(sorted_arr.sum()) - const_

        pre_gini_index = gini(X)
        pre_gini_index = round(pre_gini_index, 3)
        post_gini_index = gini(Y)
        post_gini_index = round(post_gini_index, 3)

        X_lorenz = X.cumsum() / X.sum()
        X_lorenz = np.insert(X_lorenz, 0, 0)
        X_lorenz[0], X_lorenz[-1]

        Y_lorenz = Y.cumsum() / Y.sum()
        Y_lorenz = np.insert(Y_lorenz, 0, 0)
        Y_lorenz[0], Y_lorenz[-1]

    #plots Lorenz curve
    fig, ax = plt.subplots(figsize=[6,6])
    ax.plot(np.arange(X_lorenz.size)/(X_lorenz.size-1), X_lorenz, label='Library pre-amplification',
            color='green')
    ax.plot(np.arange(Y_lorenz.size)/(Y_lorenz.size-1), Y_lorenz, label='Library post-amplification',
            color='red')
    ax.plot([0,1], [0,1], color='k', label='Ideal library')#line plot of equality
    ax.set(ylabel='Cumulative fraction of reads represented', xlabel='sgRNAs ranked by abundance')
    plt.text(0.075, 0.9, 'pre-amplification Gini index = '+str(pre_gini_index))
    plt.text(0.075, 0.85, 'post-amplification Gini index = '+str(post_gini_index))
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('../library-analysis/lorenz-curve.pdf')
    plt.close()

def go():
    pass
