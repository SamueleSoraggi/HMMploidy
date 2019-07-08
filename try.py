import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("input",help="file containing the list of basenames for gzipped     mpileup files for use in analysis to be used")
#parser.add_argument("-ft","--fileType",help="file type of the input file, mpileup o    r bam")
args = parser.parse_args()
#fileType = args.fileType
#
#out = "asd" + "." + fileType
#print(out)
#
#
#print(type(fileType))
#fileTypes = ["mpileup.gz", "bam"]
#if fileType not in fileTypes:
#    sys.exit(fileType + " is not supported")

input = args.input

with open(input,'rb') as f:
    for line in f:
        line = line.decode().strip('\n')
        if line.endswith(".bam"):
            fileType = "bam"
        elif line.endswith(".mpileup.gz"):
            fileType = "mpileup.gz"
        elif line.endswith(".mpileup"):
            fileType = "mpileup"
        else:
            fileType = line
            sys.exit(line + "... Try using .mpileup, .mpileup.gz or .bam files instead.")

print(fileType)
