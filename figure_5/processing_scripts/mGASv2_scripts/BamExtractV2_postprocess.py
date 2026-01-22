# Import modules
import gzip
import subprocess
import os
import sys
import csv
import re
import itertools
from collections import Counter, OrderedDict
from operator import itemgetter
from datetime import datetime
from optparse import OptionParser, OptionGroup

import numpy as np
import edlib

import pandas as pd

# Additional functions
def complement(seq):
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement_map.get(base, 'N') for base in seq])

def reverse_complement(seq):
    return complement(seq[::-1])

# Main processing
# Read cell_target_nRead_table from the BAM Extract file

c = 0
cell_target_nRead_table = {}

with open('mGASv2_Lane3_CellByTape_10X_bamExtractV2_t3.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        if c == 0:
            c += 1
            continue
        if c > 0:
            CellBC = row<source_id data="0" title="mEBv2_BamExtractV2_postprocess-Copy1.ipynb" />
            TargetBC = row[1]
            nRead = row[8]
            cell_target = CellBC + ',' + TargetBC
            try:
                if cell_target_nRead_table[cell_target] < nRead:
                    cell_target_nRead_table[cell_target] = nRead
            except KeyError:
                cell_target_nRead_table[cell_target] = nRead

# Construct cell_TAPE information
c = 0
cell_TAPE = {}

with open('mGASv2_Lane3_CellByTape_10X_bamExtractV2_t3.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        if c == 0:
            c += 1
            continue
        if c > 0:
            CellBC = row<source_id data="0" title="mEBv2_BamExtractV2_postprocess-Copy1.ipynb" />
            TargetBC = row[1]
            nRead = row[8]
            cell_target = CellBC + ',' + TargetBC
            # TAPE assignment logic (from your notebook code)
            if row[2][3:] == 'GGA':
                if row[3][3:] == 'GGA':
                    if row[4][3:] == 'GGA':
                        if row[5][3:] == 'GGA':
                            if row[6][3:] == 'GGA':
                                if row[7][3:] == 'GGA':
                                    TAPE = ','.join(row[2:8])
                                else:
                                    TAPE = ','.join(row[2:7]) + ',None'
                            else:
                                TAPE = ','.join(row[2:6]) + ',None,None'
                        else:
                            TAPE = ','.join(row[2:5]) + ',None,None,None'
                    else:
                        TAPE = ','.join(row[2:4]) + ',None,None,None,None'
                else:
                    TAPE = ','.join(row[2:3]) + ',None,None,None,None,None'
            else:
                TAPE = 'None,None,None,None,None,None'
            if nRead == cell_target_nRead_table[cell_target]:
                cell_TAPE[cell_target] = TAPE

# Write output CSV
with open('mGASv2_Lane3_CellByTape_10X_bamExtractV2_t3_collapse.csv', 'w', newline='') as f0:
    f0.write('Cell,TargetBC,Site1,Site2,Site3,Site4,Site5,Site6\n')
    for key, value in cell_TAPE.items():
        f0.write(key + ',' + value + '\n')

# Print a preview (first 100 lines)
c = 0
for k, v in cell_TAPE.items():
    if c < 100:
        print(k, v)
        c += 1