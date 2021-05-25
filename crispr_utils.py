#!/usr/bin/env python3

import warnings
import pkg_resources
import os
import subprocess
import multiprocessing
from zipfile import ZipFile
import tarfile
import pickle
import yaml
import sys
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set(style="whitegrid")

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep="",file=file)

def install_fastqc(script_dir):
    url="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
    download_file=os.path.join(script_dir,"fastqc_v0.11.9.zip")
    download_command="wget "+url+" --output-document="+download_file
    download=input("Download and install FastQC? y/n")

    if download in ["yes","Yes","y","Y"]:
        try:
            write2log(work_dir,download_command,"Download FastQC: ")
            subprocess.run(download_command, shell=True)
            #unzip FastQC file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)
        except:
            print("FastQC installation failed, check log")
    else:
        sys.exit("FastQC installation aborted")
    fastqc_folder=os.path.join(script_dir,"fastqc_v0.11.9",FastQC)
    return(fastqc_folder)

def install_mageck(script_dir):
    url="https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.4.tar.gz/download"
    download=input("Download and install MAGeCK? y/n")
    download_file=os.path.join(script_dir,"mageck-0.5.9.4.tar.gz")
    download_command="wget "+url+" --output-document="+download_file

    if download in ["yes","Yes","y","Y"]:
        try:
            write2log(work_dir,download_command,"Download MAGeCK: ")
            subprocess.run(download_command, shell=True)
            #unpack MAGeCK file
            tar = tarfile.open(download_file, "r:gz")
            tar.extractall()
            tar.close()
        except:
            print("MAGeCK installation failed, check log")
    else:
        sys.exit("MAGeCK installation aborted")
    mageck_folder=os.path.join(script_dir,"mageck-0.5.9.4/bin")
    return(mageck_folder)

def install_bowtie2(script_dir):
    if sys.platform in ["linux","linux2"]:
        url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-linux-x86_64.zip/download"
        download_file=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64.zip")
    elif sys.platform == "darwin":
        url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-macos-x86_64.zip/download"
        download_file=os.path.join(script_dir,"bowtie2-2.4.3-macos-x86_64.zip")
    download_command="wget "+url+" --output-document="+download_file
    download=input("Download and install Bowtie2? y/n")
    if download in ["yes","Yes","y","Y"]:
        try:
            write2log(work_dir,download_command,"Download Bowtie2: ")
            subprocess.run(download_command, shell=True)
            #unzip Bowtie2 file
            with ZipFile(download_file, 'r') as zip_ref:
                zip_ref.extractall(script_dir)
        except:
            print("Bowtie2 installation failed, check log")
    else:
        sys.exit("Bowtie2 installation aborted")
    bowtie2_folder=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64")
    return(bowtie2_folder)


def install_bagel2(script_dir):
    print("Installing BAGEL2 to "+script_dir)
    bagel2_git="https://github.com/hart-lab/bagel.git"
    clone_command="git "+"clone "+"https://github.com/hart-lab/bagel.git "+os.path.join(script_dir,"bagel2")
    write2log(work_dir,clone_command,"Clone BAGEL2 git: ")
    subprocess.run(clone_command, shell=True)

def check_env(script_dir,work_dir): #check if required binary directories are set in $PATH
    env=dict(os.environ)
    required={"FastQC","mageck-0.5.9.4","bowtie2"}

    exe_dict=dict() #to store FastQC and MAGeCK binary locations

    for i in required:
        if not i in env["PATH"]:
            find_command="find "+"$HOME "+"-name "+i
            i_dir=subprocess.check_output(find_command, shell=True)
            i_dir=i_dir.decode("utf-8")
            dir_count=i_dir.count("\n")
            if dir_count == 0:
                print("Warning: "+i+" directory not found")
                if i == "fastqc":
                    #installs FastQC and writes dir to exe_dict
                    fastqc_dir=install_fastqc(script_dir)
                    exe_dict.update({"fastqc":fastqc_dir})
                elif i == "mageck":
                    #installs MAGeCK and writes dir to exe_dict
                    mageck_dir=install_mageck(script_dir)
                    exe_dict.update({"mageck":mageck_dir})
                elif i == "bowtie2":
                    #installs Bowtie2 and writes dir to exe_dict
                    bowtie2_dir=install_bowtie2(script_dir)
                    exe_dict.update({"bowtie2":bowtie2_dir})
            elif dir_count == 1:
                if i == "mageck":
                    i_dir=i_dir.replace("\n","")
                    exe_dict.update({"mageck":i_dir})
                elif i == "bowtie2":
                    i_dir=i_dir.replace("\n","")
                    exe_dict.update({"bowtie2":i_dir})
            elif dir_count == 2:
                if i == "fastqc":
                    i_dir=list(i_dir.split("\n"))[0]
                else:
                    sys.exit("ERROR: multiple "+i+" directories found (keep only one)")
            elif dir_count > 2:
                sys.exit("ERROR: multiple "+i+" directories found (keep only one)")

    #write exe_dict to file in script_dir for future reference
    try:
        pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
    except pickle.PicklingError:
        print("Storing of dictionary with dependency locations failed")

    #open with
    #exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))

def install_python_packages(): #check for required python packages; installs if absent
    required = {"pyyaml","pandas","numpy","matplotlib","seaborn","multiqc","cutadapt"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        try:
            install_command=[python, '-m', 'pip3', 'install', *missing]
            write2log(work_dir,install_command,"Missing package installation: ")
            subprocess.check_call(install_command, stdout=subprocess.DEVNULL)
        except:
            sys.exit("ERROR: package installation failed, check log")

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
        print("Skipping "+file+" (already exists/analysed)")
        return(True)
    else:
        return(False)

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
        try:
            print("Running FastQC on raw data")
            subprocess.run(fastqc_command, shell=True)
        except:
            sys.exit("ERROR: FastQC failed, check logs")
        print("Running MultiQC")
        subprocess.run(multiqc_command)
    else:
        print("Skipping FastQC/MultiQC (already performed)")

def check_index(library,crispr_library,script_dir):
    try:
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
                try:
                    subprocess.run(bowtie2-build-command) #build index
                except:
                    sys.exit("ERROR: bpwtie2-build failed, check logs")
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")

def guide_names(library,crispr_library):
    try:
        fasta=library[crispr_library]["fasta"]
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")
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
    try:
        read_mod=library[crispr_library]["read_mod"]
        sg_length=library[crispr_library]["sg_length"]
        sg_length=str(sg_length)
        index_path=library[crispr_library]["index_path"]
        clip_seq=library[crispr_library]["clip_seq"]
        fasta=library[crispr_library]["fasta"]
        mismatch=str(mismatch)
    except KeyError:
        sys.exit("ERROR: CRISPR library not specified in command line")

    print("Aligning reads to reference (mismatches allowed: "+mismatch+")")

    #bowtie2 and bash commands (common to both trim and clip)
    bowtie2="bowtie2 --no-hd -p "+threads+" -t -N "+mismatch+" -x "+index_path+" - 2>> crispr.log | "
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
                try:
                    subprocess.run(count_command,shell=True)
                except:
                    sys.exit("ERROR: read count failed, check logs")
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
                try:
                    subprocess.run(count_command, shell=True)
                except:
                    sys.exit("ERROR: read count failed, check logs")

    #remove first line from guide count text files (bowtie2 artefact)
    count_list=glob.glob(os.path.join(work_dir,"count","*guidecounts.txt"))
    for file in count_list:
        command="sed '1d' "+file+" > "+file+".temp "+"&& mv "+file+".temp "+file
        try:
            subprocess.run(command, shell=True)
        except:
            sys.exit("ERROR: removal of first line of count file failed")

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

def mageck(work_dir,script_dir,cnv):
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

    #open config.yml
    with open(os.path.join(work_dir,"config.yml")) as file:
        config=yaml.full_load(file)

    def cnv_com(script_dir,config): #generate MAGeCK command for CNV correction
        #check if specified cell line is in CCLE data list
        cell_line_list=subprocess.check_output(["head","-1",os.path.join(script_dir,"CCLE","CCLE_copynumber_byGene_2013-12-03.txt")])
        cell_line=config["CNV-cell-line"]
        cnv_command=" --cnv-norm "+ccle_ref+" --cell-line "+cell_line
        return(cnv_command)

    #check if CNV correction is requested and perform checks
    if cnv == True:
        ccle_ref=os.path.join(script_dir,"CCLE","CCLE_copynumber_byGene_2013-12-03.txt")
        if not os.path.exists(ccle_ref):
            print("WARNING: no CCLE copy number file found")
            url=" https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_copynumber_byGene_2013-12-03.txt"
            download_command="wget --directory-prefix="+os.path.join(script_dir,"CCLE")+url
            download=input("Download CCLE copy numer file from "+url+" ? y/n")
            if download in ["yes","Yes","y","Y"]:
                write2log(work_dir,download_command,"Download CCLE file: ")
                try:
                    subprocess.run(download_command, shell=True)
                except:
                    sys.exit("ERROR: download failed, check log and url")
                cnv_command=cnv_com(script_dir,config)
                if not cell_line in cell_line_list:
                    print("ERROR: specified cell line not found in CCLE reference file")
                    print("Skipping CNV correction for MAGeCK")
                    cnv=False
            else:
                    print("Skipping CNV correction for MAGeCK")
                    cnv=False
        else:
            cnv_command=cnv_com(script_dir,config)


    #determine number of samples in count table
    header=subprocess.check_output(["head", "-1",os.path.join(work_dir,"count","counts-aggregated.tsv")])
    header=header.decode("utf-8")
    sample_count=header.count("\t") - 1
    if sample_count ==3:
        if "pre" and "post" in header:
            print("Skipping MAGeCK (only CRISPR library samples present)")
            return(None)

    #load FDR cut off from config
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

        if not file_exists(os.path.join(work_dir,"mageck",mageck_output)):
            os.makedirs(os.path.join(work_dir,"mageck",mageck_output),exist_ok=True)

        prefix=os.path.join(work_dir,"mageck",mageck_output,mageck_output)
        input=os.path.join(work_dir,"count","counts-aggregated.tsv")
        log=" 2>> "+os.path.join(work_dir,"crispr.log")
        mageck_command="mageck test -k "+input+" -t "+test_sample+" -c "+control_sample+" -n "+prefix+log
        if cnv == True:
            prefix=os.path.join(work_dir,"mageck-cnv",mageck_output,mageck_output)
            mageck_command="mageck test -k "+input+" -t "+test_sample+" -c "+control_sample+" -n "+prefix+log
            mageck_command=mageck_command+cnv_command
        write2log(work_dir,mageck_command,"MAGeCK: ")
        subprocess.run(mageck_command, shell=True)

    #plot MAGeCK hits
    file_list=glob.glob(os.path.join(work_dir,"mageck","*","*gene_summary.txt"))

    for file in file_list:
        save_path=os.path.dirname(file)
        plot_command="Rscript "+os.path.join(script_dir,"R","plot-hits.R ")+work_dir+" "+file+" mageck "+save_path+" "+mageck_output
        write2log(work_dir,plot_command,"Plot hits MAGeCK: ")
        try:
            subprocess.run(plot_command, shell=True)
        except:
            sys.exit("ERROR: plotting hits failed, check log")


def convert4bagel(work_dir,library,crispr_library): #convert MAGeCK formatted count table to BAGEL2 format
    count_table_bagel2=os.path.join(work_dir,"bagel",'counts-aggregated-bagel2.tsv')

    if not file_exists(count_table_bagel2):
        #obtain sequences of each guide
        try:
            fasta=library[crispr_library]["fasta"]
        except KeyError:
            sys.exit("ERROR: CRISPR library not specified in command line")
        df_fasta=pd.read_csv(fasta, header=None)

        #place sgRNA name and sequence is separate columns
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

        #create gene column
        df_join["gene"]=df_join["sgRNA"].str.split(pat="_").str[0]
        #df_join.rename(columns = {"sgRNA":"gene"}, inplace = True)

        #open MAGeCK formatted count table
        count_file=os.path.join(work_dir,"count","counts-aggregated.tsv")
        df_master=pd.read_csv(count_file, sep="\t")

        #format sgRNA column df_join to same format as df_master
        df_join["sgRNA"]=df_join["sgRNA"].str.split(pat="_",n=1).str[1]

        #sort data frames
        df_join=df_join.sort_values(by="sgRNA")
        df_master=df_master.sort_values(by="sgRNA")

        #merge data frames and format data frame for BAGEL2
        df_join=df_join.drop(["gene"], axis=1)#remove gene column
        df_merge=pd.merge(df_master,df_join, on="sgRNA", how="left")
        df_merge=df_merge.drop(["sgRNA"], axis=1)#remove gene column
        df_merge.rename(columns={"gene":"GENE"}, inplace=True)
        cols=list(df_merge)
        cols.insert(0,cols.pop(cols.index("SEQUENCE")))
        df_merge=df_merge.loc[:,cols]

        #save df to file
        os.makedirs(os.path.join(work_dir,"bagel"),exist_ok=True)
        df_merge.to_csv(count_table_bagel2, sep='\t',index=False)

def bagel2(work_dir,script_dir):
    #find BAGEL2 dir if not in .exe_dict_obj
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if "bagel2" not in exe_dict:
        find_command="find "+"$HOME "+"-type "+"d "+"-name "+"bagel2"
        bagel2_dir=subprocess.check_output(find_command, shell=True)
        bagel2_dir=bagel2_dir.decode("utf-8")
        dir_count=bagel2_dir.count("\n")

        if dir_count == 0:
            print("ERROR: BAGEL2 directory not found")
            download=input("Download and install BAGEL2? y/n")
            if download in ["yes","Yes","y","Y"]:
                install_bagel2(script_dir)
                bagel2_exe=os.path.join(script_dir,"bagel2","BAGEL.py")
                exe_dict.update({"bagel2":bagel2_exe})
            else:
                sys.exit("BAGEL2 install aborted")
        elif dir_count > 1:
            sys.exit("ERROR: multiple BAGEL2 directories found (keep only one)")
        elif dir_count == 1:
            bagel2_dir=bagel2_dir.replace("\n","")
            bagel2_exe=os.path.join(bagel2_dir,"BAGEL.py")
            exe_dict.update({"bagel2":bagel2_exe})
    else:
        bagel2_exe=exe_dict["bagel2"]

    #get sample names from BAGEL2 count table
    header=subprocess.check_output(["head", "-1",os.path.join(work_dir,"bagel","counts-aggregated-bagel2.tsv")])
    header=header.decode("utf-8")
    header=header.replace("\n","")
    header=list(header.split("\t"))# convert string into list

    #create dictionary that holds column name (key) and column index
    column_dict={key: i for i, key in enumerate(header)}
    column_dict={key: column_dict[key] - 1 for key in column_dict} #first sample column should have value 1

    count_table=os.path.join(work_dir,"bagel",'counts-aggregated-bagel2.tsv')

    #reference genes files for Bayes Factor calculation
    essential_genes=os.path.join(bagel2_dir,"CEGv2.txt")
    nonessential_genes=os.path.join(bagel2_dir,"NEGv1.txt")

    #load stats.config for sample comparisons
    df=pd.read_csv(os.path.join(work_dir,"stats.config"),sep=";")
    sample_number=len(df)
    sample_range=range(sample_number)

    #run BAGEL2 for each comparison in stats.config
    for i in sample_range:
        test_sample=df.loc[i]["t"]
        control_sample=df.loc[i]["c"]
        bagel2_output=test_sample+"_vs_"+control_sample
        bagel2_output_base=os.path.join(work_dir,"bagel",bagel2_output,bagel2_output)
        test_sample_column=column_dict[test_sample]
        control_sample_column=column_dict[control_sample]

        if not file_exists(os.path.join(work_dir,"bagel",bagel2_output)):
            os.makedirs(os.path.join(work_dir,"bagel",bagel2_output),exist_ok=True)

        print("Generatig fold change table for "+bagel2_output)
        fc_file=os.path.join(bagel2_output_base+".foldchange")
        if not file_exists(fc_file):
            bagel2fc_command="python3 "+bagel2_exe+" fc"+" -i "+count_table+" -o "+bagel2_output_base+" -c "+str(control_sample_column)
            write2log(work_dir,bagel2fc_command,"BAGEL2 fc: ")
            try:
                subprocess.run(bagel2fc_command, shell=True)
            except:
                sys.exit("ERROR: generation of BAGEL2 fc file failed, check log")

        print("Calculating Bayes Factors for "+bagel2_output)
        bf_file=os.path.join(bagel2_output_base+".bf")
        if not file_exists(bf_file):
            #get sample names from BAGEL2 foldchange table
            header2=subprocess.check_output(["head", "-1",os.path.join(bagel2_output_base+".foldchange")])
            header2=header2.decode("utf-8")
            header2=header2.replace("\n","")
            header2=list(header2.split("\t"))# convert string into list

            #create dictionary that holds column name (key) and column index
            column_dict2={key: i for i, key in enumerate(header2)}
            column_dict2={key: column_dict2[key] - 1 for key in column_dict2} #first sample column should have value 1
            test_sample_column2=column_dict2[test_sample]

            bagel2bf_command="python3 "+bagel2_exe+" bf"+" -i "+fc_file+" -o "+bf_file+" -e "+essential_genes+" -n "+nonessential_genes+" -c "+str(test_sample_column2)
            write2log(work_dir,bagel2bf_command,"BAGEL2 bf: ")
            try:
                subprocess.run(bagel2bf_command, shell=True)
            except:
                sys.exit("ERROR: Calculation of Bayes Factors failed, check log")

        print("Calculating precision-recall for "+bagel2_output)
        pr_file=os.path.join(bagel2_output_base+".pr")
        if not file_exists(pr_file):
            bagel2pr_command="python3 "+bagel2_exe+" pr"+" -i "+bf_file+" -o "+pr_file+" -e "+essential_genes+" -n "+nonessential_genes
            write2log(work_dir,bagel2pr_command,"BAGEL2 pr: ")
            try:
                subprocess.run(bagel2pr_command, shell=True)
            except:
                sys.exit("ERROR: Calculation of precision-recall failed, check log")

        print("Plotting BAGEL2 results for "+bagel2_output)
        plot_file=os.path.join(work_dir,"bagel",bagel2_output,"PR-"+bagel2_output+".pdf")
        if not file_exists(plot_file):
            plot_script=os.path.join(script_dir,"R","plot-hits.R")
            plot_command="Rscript "+plot_script+" "+work_dir+" "+pr_file+" bagel2 "+os.path.join(work_dir,"bagel",bagel2_output)+" "+bagel2_output
            write2log(work_dir,plot_command,"BAGEL2 plot: ")
            try:
                subprocess.run(plot_command, shell=True)
            except:
                sys.exit("ERROR: Calculation of precision-recall failed, check log")


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

def go(work_dir,script_dir):
    #load GO analysis settings
    with open(os.path.join(work_dir,"config.yml")) as file:
        config=yaml.full_load(file)

    email=config["GO"]["email"]
    fdr=config["mageck-fdr"]
    species=config["GO"]["species"]
    go_test=config["GO"]["test"]
    go_term=config["GO"]["term"]

    #get list og MAGeCK gene summary file_exists
    find_command="find "+work_dir+" -name "+"*gene_summary.txt"
    file_list=subprocess.check_output(find_command, shell=True)
    file_list=file_list.decode("utf-8")
    file_list=list(file_list.split("\n"))

    for file in file_list:
        save_path=os.path.dirname(file)
        go_command=os.path.join(script_dir,"R","go.R ")+email+" "+file+" "+fdr+" "+species+" "+save_path+" "+go_test+" "+go_term
        write2log(work_dir,go_command,"MAGeCK GO: ")
        try:
            subprocess.run(go_command, shell=True)
        except:
            sys.exit("ERROR: GO analysis failed, check log")