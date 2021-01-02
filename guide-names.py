#!/usr/bin/env python3
# -*- coding: utf-8 -*-
docstring="""
@author: Niek Wit (University of Cambridge)
This function returns a .csv file of guides names extracted from the FASTA CRISPR library file
Run as follows from the command line: $ python3 /path/to/guide-names.py library.fasta"""

import pandas as pd
import sys

if len(sys.argv) != 2:
    sys.exit('\nJust one argument required (see instructions)%s' %(docstring))

library = pd.read_csv(sys.argv[1], names = ['guide'])
#creates new dataframe with only guide names:
library = library[library['guide'].str.contains('>')]
#removes '<' character from each row:
library = library['guide'].str.strip('>')
#saves guide names to a .csv file:
library_name=sys.argv[1]
library_name=library_name.replace('.fasta','')
library.to_csv(library_name+'_guides.csv', index=False)    
