#!/bin/bash

#This script enables autocompletion for the CRISPR library options (-l) in crispr-pipeline.sh

function libs()
{
  # $1 is the name of the command 
  # $2 is the word being completed
  # $3 is the word preceding the word being completed

case $3 in
    -l) COMPREPLY=($(compgen -W "bassik moffat_tko1 moffat_tko3 sabatini dub-only" "${COMP_WORDS[$COMP_CWORD]}"));;
  esac
}

complete -F libs crispr-pipeline #refers to symbolic link in /usr/bin
