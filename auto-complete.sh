#!/bin/bash

#This script enables autocompletion for the CRISPR library options (-l) in crispr-pipeline.sh
#Add the following line to your .bashrc file: source /path/to/CRISPR-tools/auto-complete.sh

#finds CRISPR libary names
SCRIPT_DIR=$(find $HOME -type d -name "CRISPR-tools")
lib_list=$(cat "$SCRIPT_DIR/config.yml" | shyaml keys | tr "\n" " ")
lib_list=$(echo ${lib_list/mageck-fdr /})

#enables autocompletion of `-l` flag
function libs()
{
case $3 in
    -l) COMPREPLY=($(compgen -W "$lib_list" "${COMP_WORDS[$COMP_CWORD]}"));;
  esac
}

complete -F libs crispr-pipeline.sh
