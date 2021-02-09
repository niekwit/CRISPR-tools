#!/bin/bash

#This script enables autocompletion for the CRISPR library options (-l) in crispr-pipeline.sh
#Add the following line to your .bashrc file: source /path/to/CRISPR-tools/auto-complete.sh

function libs()
{
  # $1 is the name of the command 
  # $2 is the word being completed
  # $3 is the word preceding the word being completed

case $3 in
    -l) COMPREPLY=($(compgen -W "bassik moffat_tko1 moffat_tko3 sabatini dub-only slc-mito-2ogdd" "${COMP_WORDS[$COMP_CWORD]}"));;
  esac
}

complete -F libs crispr-pipeline.sh
