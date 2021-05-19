#!/usr/bin/env python3

import warnings
import pkg_resources
import os
import subprocess
import multiprocessing
import yaml
import sys
import csv
import glob
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

def set_threads(args):
    max_threads=str(multiprocessing.cpu_count())
    threads=args["threads"]
    if threads == "max":
        threads=max_threads
    threads=str(threads)
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
        print("\tSkipping "+file+" (already exists/analysed)")
        return(True)
    else:
        return(False)

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep="",file=file)

def fastqc(work_dir,threads,file_extension):
    if not os.path.isdir(os.path.join(work_dir,"fastqc")) or len(os.listdir(os.path.join(work_dir,"fastqc"))) == 0:
        os.makedirs(work_dir+"/fastqc",exist_ok=True)
        fastqc_command="fastqc --threads "+str(threads)+" --quiet -o fastqc/ raw-data/*"+file_extension
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

def check_index(library,crispr_library,script_dir):
    index_path=library[crispr_library]["index_path"]
    fasta=library[crispr_library]["fasta"]

    print(crispr_library+" library selected")

    if index_path == "":
        print("No index file found for "+crispr_library)
        if fasta == "":
            sys.exit("ERROR:No fasta file found for "+crispr_library)
        else:
            index_base=os.path.join(script_dir,"index",crispr_library,crispr_library+"-index")
            bowtie2_build_command=["bowtie2-build",fasta,index_base]
            write2log(work_dir,bowtie2_build_command,"Bowtie2-build: ")

            #Write bowtie2 index file location to library.yaml
            with open(os.path.join(script_dir,"library.yaml")) as f:
                doc=yaml.safe_load(f)
            doc[crispr_library]["index_path"]=index_base
            with open(os.path.join(script_dir,"library.yaml"), "w") as f:
                yaml.dump(doc,f)

            print("Building Bowtie2 index for "+crispr_library+"library")
            subprocess.run(bowtie2-build-command) #build index

def guide_names(library,crispr_library):
    fasta=library[crispr_library]["fasta"]
    output_name=fasta.replace(".fasta","-guide_names.csv")

    if not os.path.exists(output_name):
        library = pd.read_csv(fasta, names = ['guide'])
        #creates new dataframe with only guide names:
        library = library[library['guide'].str.contains('>')]
        #removes '<' character from each row:
        library = library['guide'].str.strip('>')
        #saves guide names to a .csv file:
        library.to_csv(output_name, index=False, header=False)

def count(library,crispr_library,mismatch,threads,script_dir,work_dir,file_extension):
    os.makedirs(os.path.join(work_dir,"count"),exist_ok=True)

    read_mod=library[crispr_library]["read_mod"]
    sg_length=library[crispr_library]["sg_length"]
    sg_length=str(sg_length)
    index_path=library[crispr_library]["index_path"]
    clip_seq=library[crispr_library]["clip_seq"]
    fasta=library[crispr_library]["fasta"]
    mismatch=str(mismatch)

    print("Aligning reads to reference (mismatches allowed: "+mismatch+")")

    #bowtie2 and bash commands (common to both trim and clip)
    bowtie2="bowtie2 --no-hd -p"+threads+" -t -N "+mismatch+" -x "+index_path+" - 2>> crispr.log | "
    bash="sed '/XS:/d' | cut -f3 | sort | uniq -c > "

    #trim, align and count
    if read_mod == "trim":
        file_list=glob.glob(os.path.join(work_dir,"raw-data","*"+file_extension))
        for file in file_list:
            base_file=os.path.basename(file)
            out_file=os.path.join(work_dir,"count",base_file.replace(file_extension,".guidecounts.txt"))
            if not file_exists(out_file):
                print("Aligning "+base_file)
                print(base_file+":", file=open("crispr.log", "a"))
                cutadapt="cutadapt -j "+threads+" --quality-base 33 -l "+sg_length+" -o - "+file+" 2>> crispr.log | "
                cutadapt=str(cutadapt)
                bowtie2=str(bowtie2)
                count_command=cutadapt+bowtie2+bash+out_file
                write2log(work_dir,count_command,"Count: ")
                subprocess.run(count_command,shell=True)
    elif read_mod == "clip":
        file_list=glob.glob(os.path.join(work_dir,"raw-data","*"+file_extension))
        for file in file_list:
            base_file=os.path.basename(file)
            out_file=os.path.join(work_dir,"count",file.replace(base_file,".guidecounts.txt"))

            if not file_exists(out_file):
                print("Aligning "+base_file)
                print(base_file, file=open("crispr.log", "a"))
                cutadapt="cutadapt -j "+threads+" --quality-base 33 -a "+clip_seq+" -o - "+file+" 2>> crispr.log | "
                cutadapt=str(cutadapt)
                bowtie2=str(bowtie2)
                count_command=cutadapt+bowtie2+bash+out_file
                write2log(work_dir,count_command,"Count: ")
                subprocess.run(count_command, shell=True)

    #remove first line from guide count text files (bowtie2 artefact)
    count_list=glob.glob(os.path.join(work_dir,"count","*guidecounts.txt"))
    for file in count_list:
        command="sed '1d' "+file+" > "+file+".temp "+"&& mv "+file+".temp "+file
        subprocess.run(command, shell=True)

def normalise(work_dir):
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

def join_counts(work_dir,library,crispr_library):
    #load sgRNA names, used for merging data2
    fasta=library[crispr_library]["fasta"]
    guide_name_file=fasta.replace(".fasta","-guide_names.csv")

    sgrnas_list00 = list(csv.reader(open(guide_name_file)))

    sgrnas_list0 = []

    for x in sgrnas_list00: #Flattens the list
        for y in x:
            sgrnas_list0.append(y)

    #Generates sgRNA and gene columns for final output
    sgRNA_output = []
    gene_output = []

    for n in sgrnas_list0:
        #print(n)
        s,g = n.split("_", 1)
        sgRNA_output.append(g)
        gene_output.append(s)

    #Generates reference Pandas data frame from sgRNA list library file
    d0 = {'sgRNA':pd.Series(sgRNA_output),'gene':pd.Series(gene_output),'sgRNA2':pd.Series(sgrnas_list0)}
    dfjoin1 = pd.DataFrame(d0) #sgRNA/gene column required for MAGeCK, sgRNA2 is needed for join operation (deleted later)

    #Generates a list of all count .txt files
    file_list = glob.glob(os.path.join(work_dir,"count",'*.guidecounts.txt'))
    file_list.sort()
    file_list2 = [w.replace('.guidecounts.txt','') for w in file_list] #this list will generate the column headers for the output file (removes .txt)

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

    #only keep base name as column names
    column_number=len(dfjoin2.columns)
    column_range=range(2,column_number)
    for i in column_range:
        old_name=dfjoin2.columns[i]
        new_name=os.path.basename(old_name)
        dfjoin2.rename(columns={list(dfjoin2)[i]:new_name},inplace=True)

    #Writes all data to a single .tsv file, ready for MAGeCK
    dfjoin2.to_csv(os.path.join(work_dir,"count",'counts-aggregated.tsv'), sep='\t',index=False)

def mageck(work_dir,script_dir):
    #check for stats.config
    stats_config=os.path.join(work_dir,"stats.config")
    if not os.path.exists(stats_config):
        print("ERROR: stats.config not found (MAGeCK comparisons)")
        return(None)

    #check for config.yml
    config_yml=os.path.join(work_dir,"config.yml")
    if not os.path.exists(config_yml):
        print("ERROR: config.yml not found (MAGeCK settings)")
        return(None)

    #determine number of samples in count table
    header=subprocess.check_output(["head", "-1",os.path.join(work_dir,"count","counts-aggregated.tsv")])
    header=header.decode("utf-8")
    sample_count=header.count("\t") - 1
    if sample_count ==3:
        if "pre" and "post" in header:
            print("Skipping MAGeCK (only CRISPR library samples present)")
            return(None)

    #open config.yml
    with open(os.path.join(work_dir,"config.yml")) as file:
        config=yaml.full_load(file)

    fdr=config["mageck-fdr"]

    #create MAGeCK dir
    os.makedirs(os.path.join(work_dir,"mageck"),exist_ok=True)

    #load MAGeCK comparisons and run MAGeCK
    df=pd.read_csv(os.path.join(work_dir,"stats.config"),sep=";")
    sample_number=len(df)
    sample_range=range(sample_number)

    for i in sample_range:
        test_sample=df.loc[i]["t"]
        control_sample=df.loc[i]["c"]
        mageck_output=test_sample+"_vs_"+control_sample
        print(mageck_output)

        if not file_exists(os.path.join(work_dir,"mageck",mageck_output)):
            os.makedirs(os.path.join(work_dir,"mageck",mageck_output),exist_ok=True)

        prefix=os.path.join(work_dir,"mageck",mageck_output,mageck_output)
        input=os.path.join(work_dir,"count","counts-aggregated.tsv")
        log=" 2>> "+os.path.join(work_dir,"crispr.log")
        mageck_command="mageck test -k "+input+" -t "+test_sample+" -c "+control_sample+" -n "+prefix+log
        write2log(work_dir,mageck_command,"MAGeCK: ")
        subprocess.run(mageck_command, shell=True)

    #plot MAGeCK hits
    file_list=glob.glob(os.path.join(work_dir,"mageck","*","*gene_summary.txt"))

    for file in file_list:
        save_path=os.path.dirname(file)
        plot_command="Rscript "+os.path.join(script_dir,"plot-hits.R ")+work_dir+" "+file+" mageck "+save_path
        write2log(work_dir,plot_command,"Plot hits MAGeCK: ")
        subprocess.run(plot_command, shell=True)


def convert4bagel(work_dir,library):
    #obtain sequences of each guide
    fasta=library[crispr_library]["fasta"]
    df_fasta=pd.read_csv(fasta, header=None)

    df_name=df_fasta[df_fasta[0].str.contains(">")]
    names=df_name.squeeze()#convert to series
    names=names.reset_index(drop=True)#reset index
    names=names.str.replace(">","")#remove >
    names.name="sgRNA"
    df_seq=df_fasta[~df_fasta[0].str.contains(">")]
    seq=df_seq.squeeze()#convert to series
    seq=seq.reset_index(drop=True)#reset index
    seq.name="SEQUENCE"
    df_join=pd.concat([names,seq],axis=1)#create df with names and sequences
    #reformat sgRNA to just gene name

    #open MAGeCK formatted count table
    count_file=os.path.join(work_dir,"count","counts-aggregated.tsv")
    df_master=pd.read_csv(count_file, sep="\t")


def bagel2():
    pass

def ceres():
    pass

def lib_analysis(work_dir):
    #determine whether count file contains library samples pre and post
    header=subprocess.check_output(["head", "-1",os.path.join(work_dir,"count","counts-aggregated.tsv")])
    if "pre" and "post" in str(header):
        #check if analysis has been performed already
        out_file=os.path.join(work_dir,"library-analysis","normalised-guides-frequency.pdf")
        if not os.path.exists(out_file):
            print("Analysing CRISPR library quality")
        else:
            print("CRISPR library quality already analysed")
            return(None)

        os.makedirs(os.path.join(work_dir,"library-analysis"),exist_ok=True)
        warnings.filterwarnings("ignore")

        df = pd.read_csv(os.path.join(work_dir,"count",'counts-aggregated.tsv'),sep='\t')

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
        plt.savefig(os.path.join(work_dir,"library-analysis","normalised-guides-frequency.pdf"))
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
        plt.savefig(os.path.join(work_dir,"library-analysis","normalised-pre-amplification-post-amplification.pdf"))
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
        plt.savefig(os.path.join(work_dir,"library-analysis","lorenz-curve.pdf"))
        plt.close()
    else:
        return None

def go():
    pass
