#!/usr/bin/env python3

import subprocess
import pkg_resources
import os

def write2log(work_dir,command,name):
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep="",file=file)

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

if __name__ == "__main__":
    work_dir=os.getcwd()
    install_python_packages()
