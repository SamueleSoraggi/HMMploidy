#!usr/bin/python3

import pybam
import sys

bam_file = sys.argv[1]

# bam_data = pybam.bgunzip(bam_file)

for alignment in pybam.read(bam_file):
    print(alignment.sam_seq + ' ' + alignment.sam_qual)
