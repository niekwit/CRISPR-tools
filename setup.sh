#!/bin/bash

#Installs Python3 dependencies
pip3_path=$(which pip3)

if [ -z ${pip3_path} ];
	then
		echo "pip3 not found, please install with 'sudo apt install python3-pip'"
		exit 1
fi

python_deps="pandas numpy matplotlib seaborn"

for dep in $python_deps
	do
		pip3 install $dep
	done

#Installs R dependencies
r_cran_deps="ggplot2 ggrepel plyr dplyr BiocManager"

for dep in $r_cran_deps
	do
		echo "install.packages('${dep}')" | R --no-save
	done

r_biocman_deps="org.Hs.eg.db org.Mm.eg.db"

for dep in $r_biocman_deps
	do
		echo "BiocManager::install('${dep}')" | R --no-save
	done
