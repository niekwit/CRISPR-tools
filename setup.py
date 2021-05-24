#!/usr/bin/env python3

import subprocess
import pkg_resources
import os
import sys

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

if __name__ == "__main__":
    work_dir=os.getcwd()
    script_dir=os.path.abspath(os.path.dirname(__file__))

    install_python_packages()
    check_env(script_dir,work_dir)
