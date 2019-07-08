import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-ft","--fileType",help="file type of the input file, mpileup o    r bam")
args = parser.parse_args()
fileType = args.fileType

out = "asd" + "." + fileType
print(out)

fileTypes = ["mpileup.gz", "bam"]
if fileType not in fileTypes:
    sys.exit(fileType + " is not supported")
