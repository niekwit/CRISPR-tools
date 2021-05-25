#!/usr/bin/env python3

import subprocess
import pkg_resources
import os
import sys
import pickle

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep="",file=file)

def install_python_packages(): #check for required python packages; installs if absent
    required = {"shyaml","pyyaml","pandas","numpy","matplotlib","seaborn","multiqc","cutadapt"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        try:
            install_command=[python, '-m', 'pip', 'install', *missing]
            #write2log(work_dir,install_command,"Missing package installation: ")
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
        download_command="wget "+url+" --output-document="+download_file
        write2log(work_dir,download_command,"Download FastQC: ")
        subprocess.run(download_command, shell=True)
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

def install_mageck(script_dir,mageck_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "mageck" in exe_dict:
        url="https://sourceforge.net/projects/mageck/files/0.5/mageck-0.5.9.4.tar.gz/download"
        download_file=os.path.join(script_dir,"mageck-0.5.9.4.tar.gz")
        download_command="wget "+url+" --output-document="+download_file
        subprocess.run(download_command, shell=True)
        #unpack MAGeCK file
        tar = tarfile.open(download_file, "r:gz")
        tar.extractall()
        tar.close()
        #add MAGeCK directory to exe_dict
        exe_dict["mageck"]=mageck_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of MAGeCK dir to dictionary with dependency locations failed")

def install_bowtie2(script_dir,bowtie2_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "bowtie2" in exe_dict:
        if sys.platform in ["linux","linux2"]:
            url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-linux-x86_64.zip/download"
            download_file=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64.zip")
        elif sys.platform == "darwin":
            url="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-macos-x86_64.zip/download"
            download_file=os.path.join(script_dir,"bowtie2-2.4.3-macos-x86_64.zip")
        download_command="wget "+url+" --output-document="+download_file
        subprocess.run(download_command, shell=True)
        #unzip Bowtie2 file
        with ZipFile(download_file, 'r') as zip_ref:
            zip_ref.extractall(script_dir)
        #add Bowtie directory to exe_dict
        exe_dict["bowtie2"]=bowtie2_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of Bowtie2 dir to dictionary with dependency locations failed")

def install_bagel2(script_dir,bagel2_dir):
    exe_dict=pickle.load(open(os.path.join(script_dir,".exe_dict.obj"),"rb"))
    if not "bagel2" in exe_dict:
        print("Installing BAGEL2 to "+script_dir)
        bagel2_git="https://github.com/hart-lab/bagel.git"
        clone_command="git "+"clone "+"https://github.com/hart-lab/bagel.git "+os.path.join(script_dir,"bagel2")
        subprocess.run(clone_command, shell=True)
        #add Bowtie directory to exe_dict
        exe_dict["bagel2"]=bowtie2_dir
        #save exe_dict to file
        try:
            pickle.dump(exe_dict, file=open(os.path.join(script_dir,".exe_dict.obj"),"wb"))
        except pickle.PicklingError:
            print("Storing of BAGEL2 dir to dictionary with dependency locations failed")

def check_env(script_dir,work_dir):
    fastqc_dir=os.path.join(script_dir,"fastqc_v0.11.9","FastQC")
    mageck_dir=os.path.join(script_dir,"mageck-0.5.9.4","bin")
    bagel2_dir=os.path.join(script_dir,"bagel2")
    if sys.platform in ["linux","linux2"]:
        bowtie2_dir=os.path.join(script_dir,"bowtie2-2.4.3-linux-x86_64")
    elif sys.platform == "darwin":
        bowtie2_dir=os.path.join(script_dir,"bowtie2-2.4.3-macos-x86_64")

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
    install_bowtie2(script_dir,bowtie2_dir)
    install_bagel2(script_dir,bagel2_dir)

if __name__ == "__main__":
    work_dir=os.getcwd()
    script_dir=os.path.abspath(os.path.dirname(__file__))

    install_python_packages()
    check_env(script_dir,work_dir)
