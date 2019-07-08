import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-ft","--fileType",help="file type of the input file, mpileup o    r bam")
args = parser.parse_args()
fileType = args.fileType

if fileType == "bam":
    print("y")
else:
    print("n")
