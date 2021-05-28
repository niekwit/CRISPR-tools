#!/usr/bin/env python3

import glob
import subprocess
import pkg_resources
import os
import stat
import sys
import pickle
from zipfile import ZipFile
import urllib.request
import shutil

def install_python_packages(): #check for required python packages; installs if absent
    required = {"shyaml","pyyaml","pandas","numpy","matplotlib","seaborn","multiqc","cutadapt","scipy","scikit-learn"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        try:
            install_command=[python, '-m', 'pip', 'install', *missing]
            subprocess.check_call(install_command, stdout=subprocess.DEVNULL)
        except:
            sys.exit("ERROR: package installation failed, check log")
    else:
        print("All required Python3 packages already installed")

def install_fastqc(script_dir,fastqc_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "fastqc" in exe_dict:
        url="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
        download_file=os.path.join(script_dir,"fastqc_v0.11.9.zip")
        print("Installing FastQC to "+script_dir)
        urllib.request.urlretrieve(url,download_file)
        #unzip FastQC file
        with ZipFile(download_file, 'r') as zip_ref:
            zip_ref.extractall(script_dir)
        #add FastQC directory to exe_dict
        exe_dict["fastqc"]=fastqc_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of FastQC dir to dictionary with dependency locations failed")
        #chmod +x fastqc
        fastqc_file=os.path.join(fastqc_dir,"fastqc")
        st = os.stat(fastqc_file)
        os.chmod(fastqc_file, st.st_mode | stat.S_IEXEC)
        #remove download file
        os.remove(download_file)
    else:
        print("FastQC already installed")
        return(None)

def install_mageck(script_dir,mageck_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "mageck" in exe_dict:
        url="https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.4.tar.gz/download"
        download_file=os.path.join(script_dir,"mageck-0.5.9.4.tar.gz")
        print("Installing MAGeCK to "+script_dir)
        urllib.request.urlretrieve(url,download_file)
        #unpack MAGeCK file
        tar_command="tar -xzf "+ download_file
        subprocess.run(tar_command,shell=True)
        curr_dir=os.getcwd()
        os.chdir(os.path.join(script_dir,"mageck-0.5.9.4"))
        mageck_install_command="python3 "+os.path.join(script_dir,"mageck-0.5.9.4","setup.py"+" install --user")
        subprocess.run(mageck_install_command,shell=True)
        os.chdir(curr_dir)
        #add MAGeCK directory to exe_dict
        exe_dict["mageck"]=mageck_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of MAGeCK dir to dictionary with dependency locations failed")
        #remove download file and folder
        os.remove(download_file)
        shutil.rmtree(os.path.join(script_dir,"mageck-0.5.9.4"))
    else:
        print("MAGeCK already installed")
        return(None)

def install_bowtie2(script_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "bowtie2" in exe_dict:
        if sys.platform in ["linux","linux2"]:
            url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-linux-x86_64.zip/download"
            download_file=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64.zip")
            bowtie2_dir=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64")
        elif sys.platform == "darwin":
            url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-macos-x86_64.zip/download"
            download_file=os.path.join(script_dir,"bowtie2-2.4.3-macos-x86_64.zip")
            bowtie2_dir=os.path.join(script_dir,"bowtie2-2.4.3-macos-x86_64")
        print("Installing Bowtie2 to"+script_dir)
        urllib.request.urlretrieve(url,download_file)
        #unzip Bowtie2 file
        unzip_command="unzip -qq "+download_file+" -d "+script_dir
        subprocess.run(unzip_command,shell=True)

        #add Bowtie directory to exe_dict
        exe_dict["bowtie2"]=bowtie2_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of Bowtie2 dir to dictionary with dependency locations failed")

        #remove download file
        os.remove(download_file)
    else:
        print("Bowtie2 already installed")
        return(None)

def install_bagel2(script_dir,bagel2_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "bagel2" in exe_dict:
        print("Installing BAGEL2 to "+script_dir)
        bagel2_git="https://github.com/hart-lab/bagel.git"
        clone_command="git "+"clone --quiet "+"https://github.com/hart-lab/bagel.git "+os.path.join(script_dir,"bagel2")
        subprocess.run(clone_command, shell=True)
        #add Bowtie directory to exe_dict
        exe_dict["bagel2"]=bagel2_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of BAGEL2 dir to dictionary with dependency locations failed")
    else:
        print("BAGEL2 already installed")
        return(None)

def install_R_packages():
    install_command="Rscript "+os.path.join("R","install_R_packages.R")
    subprocess.run(install_command,shell=True)

def check_env(script_dir,work_dir):
    fastqc_dir=os.path.join(script_dir,"FastQC")
    mageck_dir=os.path.join(script_dir,"mageck-0.5.9.4","bin")
    bagel2_dir=os.path.join(script_dir,"bagel2")

    #create dictionary for storing paths and save to file
    exe_dict=dict()
    exe_dict_file=os.path.join(script_dir,".exe_dict.obj")
    try:
        pickle.dump(exe_dict, file=open(exe_dict_file,"wb"))
    except pickle.PicklingError:
        print("Storing of dictionary with dependency locations failed")

    #download dependencies and store path to exe_dict
    install_fastqc(script_dir,fastqc_dir)
    install_mageck(script_dir,mageck_dir)
    install_bowtie2(script_dir)
    install_bagel2(script_dir,bagel2_dir)
    print("Installation finished")

if __name__ == "__main__":
    work_dir=os.getcwd()
    script_dir=os.path.abspath(os.path.dirname(__file__))

    install_python_packages()
    install_R_packages()
    check_env(script_dir,work_dir)
