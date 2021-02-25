#! usr/bin/python3
import os
import sys
import gzip
import generics
import numpy as np
import math
import scipy.stats
import random
from statistics import mode
import argparse

class Site:
    def __init__(self,chrom,position,reference):
        self.chrom = str(chrom)
        self.position = int(position)
        self.reference = str(reference)

class Reads:
    def __init__(self,base,base_quality):
        self.base = str(base)
        self.base_quality = str(base_quality)

alleles = ['A','C','G','T']

parser = argparse.ArgumentParser()
parser.add_argument("input",help="File containing the list of file names for gzipped mpileup files, mpileup files or bam files to be used in analysis")
parser.add_argument("-p", "--ploidyFile",help="File containing the list of ploidy levels to be used in analysis")
parser.add_argument("-o","--outFolder",help="output folder",default=0)
parser.add_argument("-i","--Inbreeding",help="Inbreeding coefficients for samples e.g 0.1x3,0.2 = 0.1,0.1,0.1,0.2 ")
parser.add_argument("-d","--downsampling",help="Fraction of data to be used in the calculations",default=1)
parser.add_argument("-m","--min_non_major_freq",type=float,help="Set the minimum frequency of non major alleles for bases to be included in the calculations",default=0.2)
parser.add_argument("-q","--min_quality_score",type=int,help="Set the minimum quality score of a read to be included in the calculation",default=1)
parser.add_argument("-dp","--min_global_depth",type=float,help="Set the minimum global depth of a base to be included in calculations",default=0)
parser.add_argument("-dpInd","--min_ind_depth",type=float,help="Set the minimum individual depth of a base to be included in calculations",default=0)
parser.add_argument("-M2","--max_minor2_freq",type=float,help="Set the maximum frequency of third most prolific alleles for bases to be included in the calculations",default=0.1)
parser.add_argument("-M3","--max_minor3_freq",type=float,help="Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations",default=0.1)
parser.add_argument("-s","--random_seed",type=int,help="Set the random seed to be included in the calculations")
parser.add_argument("-r","--ref_fasta",help="File name of the reference fasta file. Use this if bam files are being used. e.g TAIR10.fa")
args = parser.parse_args()

if args.random_seed:
    seed = args.random_seed
    random.seed(seed) # set the seed for calculations using random function
    print("Seed for calculations using random function is: " + str(seed))
else:
    print("Seed is not set.")

if args.ploidyFile:
    ploidyFile = args.ploidyFile
    list_of_ploidy=[]
    with open(ploidyFile,'rb') as ploidyList:
        for line in ploidyList:
            line = int(line.decode().strip('\n')) # convert bytes into integers
            list_of_ploidy.append(line)
        ploidy = list_of_ploidy
        print("Ploidy levels to be tested in analysis are: " + str(ploidy))
else:
    ploidy = [1,2,3,4,5,6] # default ploidy levels
    print("Default ploidy levels to be tested in analysis are: " + str(ploidy))

inputs = args.input # input file in form of mpileup, gzipped mpileup or bam
list_of_inputs=[]
fileType = 0 # initial file type
fileTypes = {
        1:"bam",
        2:"mpileup",
        3:"mpileup.gz"
        }

with open(inputs,'rb') as filelist:
    for line in filelist:
        line = line.decode().strip('\n') # convert bytes into strings
        if line.endswith(".bam"):
            if fileType != 1:
                if fileType != 0: # if it is no the first line
                    sys.exit("Error in file " + line + ". Input file should contain files from the same file type.")
            fileType = 1 # bam
            list_of_inputs.append(line)
        elif line.endswith(".mpileup"):
            if fileType != 2:
                if fileType != 0:
                    sys.exit("Error in file " + line + ". Input file should contain files from the same file type.")
            fileType = 2 # mpileup
            list_of_inputs.append(line)
        elif line.endswith(".mpileup.gz"):
            if fileType != 3:
                if fileType != 0:
                    sys.exit("Error in file " + line + ". Input file should contain files from the same file type.")
            fileType = 3 # mpileup.gz
            list_of_inputs.append(line)
        else:
            sys.exit(line + " file is not supported. Supported file types are '.mpileup', '.mpileup.gz' and '.bam'.")
    Nfiles=len(list_of_inputs)
    if Nfiles == 1:
        print('1 file found')
    else:
        print('%d files found' %Nfiles)
    print(list_of_inputs)
print("File type detected: " + fileTypes[fileType])

extensionLen = {
        1:4, #.bam
        2:8, #.mpileup
        3:11 #.mpileup.gz
        }

exl = extensionLen[fileType] # extension length of fileType including dots
outFolder = args.outFolder
    
if fileType == 1: # bam file
    if args.ref_fasta:
        refFile = args.ref_fasta
    else:
        sys.exit("Reference fasta for bam file is not found. Please add your reference file in fasta format by using '-r' parameter and try again.")
    for g1 in list_of_inputs: # for every filename 
        directory = '/'.join(g1.split('/')[:-1])
        if len(directory) == 0:
            g = "./" + g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = "./" + g2 + ".genolikes.gz"
                mpFile = "./" + g2 + ".mpileup"
            else:
                output = outFolder+'/'.join(g2.split('/')[-1])+".genolikes.gz"
                mpFile = outFolder+'/'.join(g2.split('/')[-1])+".mpileup"
        else:
            g = g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = g2 + ".genolikes.gz"
                mpFile = g2 + ".mpileup"
            else:
                output = outFolder+"/"+g2+".genolikes.gz"
                mpFile = outFolder+"/"+g2+".mpileup"
        print("Output file is: " + output)
        bashCommand = "samtools mpileup -f " + refFile + " "  + g1 + " > " + mpFile
        try:
            os.system(bashCommand)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
            # handle file not found error
                sys.exit("'samtools' not found. Please make sure that you have samtools installed.")
            else:
            # Something else went wrong while trying to run `samtools`
                sys.exit("Something went wrong with 'samtools'.")
        with open(mpFile) as mp:
            first_line = mp.readline()
            Data = first_line.strip('\n')
            l = Data.split('\t')
            NSAMS=int((len(l)-3)/3)
            if args.Inbreeding:
                inbreed = args.Inbreeding
            else:
                inbreed = "0x{}".format(str(NSAMS))
            # parse inbreeding coeffients
            F=[]
            temp=inbreed.split(',')
            for t in temp:
                vals=t.split('x')
                if len(vals)==2:
                    F+=list(np.repeat(float(vals[0]),int(vals[1])))
                else:
                    F+=[vals[0]]
            print(F)
            downsampling=float(args.downsampling) # fraction of data to be used (0-1]. 
            #Original_sample_number=NSAMS
            win=50 # window size for calculating ploidy
            phredscale=33
            NUMSITES=np.zeros(NSAMS+1,int)
            NUMSITES_HWE=np.zeros(NSAMS+1,int)
            ExpectedPloidy=[[] for i in range(NSAMS+1)]
            Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float) # (NSAMS+1)xploidies array for probabilities of each ploidy for each sample and overall ploidy probabilities
            Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
            delta_prob=np.zeros((NSAMS+1,len(ploidy)),float)
            counts=np.zeros((NSAMS+1,len(ploidy)),float)+1 # counts of bases for each ploidy being most likely
            base_number=0 # count for bases
            list_of_window=[]
            list_of_window2=[]
            no_bases=0
            total_bases=0
            
            for line in mp:
                Data=line.strip('\n')# convert bytes into string
                l = Data.split('\t')
                mySite = Site(str(l[0]),int(l[1]),str(l[2]))
                myReads = Reads("","")
                # pooled reads for first level filtering (global depth) and estimation of minor/major alleles and allele frequencies
                individualDepth = np.zeros(NSAMS,float)
                for n in range(NSAMS):
                    n=n+1
                    subReads= Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                    individualDepth[n-1] = len(subReads.base)
                    myReads.base=str(myReads.base+subReads.base)
                    myReads.base_quality=str(myReads.base_quality+subReads.base_quality)

                # convert to bases
                [bases,indexDelN] = generics.convertSyms(myReads,mySite)
                myReads=Reads(bases,myReads.base_quality)
                if len(myReads.base)!=len(myReads.base_quality):
                    sys.exit("Conversion not succesful")
                # filter by quality
                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                myReads=Reads(bases,qualities)

                # find all indexes of occurances to be filtered out
                index_of_X=[]
                index=-1
                while True:
                        index=myReads.base.find('X',index+1)
                        if index == -1:
                            break  # all occurrences have been found
                        index_of_X.append(index)
                myReads.base=myReads.base.replace('X','')
                # remove all corresponding base qualities
                count=0
                for i in index_of_X:
                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                    count+=1
                globalDepth = len(myReads.base)
                if ((globalDepth > args.min_global_depth) & (min(individualDepth)>args.min_ind_depth)):
                    total_bases+=globalDepth
                    no_bases+=1
                    # counts of non-major bases
                    nonMajorCount = generics.calcNonMajorCounts(myReads)
                    nonMajorProp = nonMajorCount/len(myReads.base)
                    # filter the site based on global depth
                    Set_min_prop = args.min_non_major_freq # minimum proportion of nonMajorCount
                    if nonMajorProp>Set_min_prop: # remove bases where more that 1-Set_min_prop are major allele i.e monomorphic bases
                        prob_of_ancestral_allelle_maj=1-nonMajorProp

                        haploid = generics.calcGenoLogLike1(myReads,mySite)
                        tri_ref = haploid[4] # retrieve reference value for if the base is not triallilic
                        haploid=haploid[:4] # remove reference value
                        # keep reference allele as one possible allele so always assume KeepRef=0
                        [major,minor,minor2,minor3] = [haploid.index(sorted(haploid,reverse=True)[0]),haploid.index(sorted(haploid,reverse=True)[1]),haploid.index(sorted(haploid,reverse=True)[2]),haploid.index(sorted(haploid,reverse=True)[3])]
                        # remove sites with >0.1 frequency of minor 2 or minor 3 allele to remove non biallilic sites (0.1 error built in for sequencing error)
                        minor2_prop=generics.calcAlleleFreq(minor2,myReads)/len(myReads.base) # Calculate allele frequencies of minor2&3 alleles
                        minor3_prop=generics.calcAlleleFreq(minor3,myReads)/len(myReads.base)

                        if(minor2_prop<args.max_minor2_freq and minor3_prop<args.max_minor3_freq):
                            P_bar=0
                            Q_bar=0
                            for read in range(len(myReads.base)):
                                if myReads.base[read]==alleles[major]:
                                    P_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))
                                elif myReads.base[read]==alleles[minor]:
                                    Q_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))

                            P = P_bar/(P_bar+Q_bar) # proportion of major allele weigted by read quality
                            Q = Q_bar/(P_bar+Q_bar) # proportion of minor allele weigted by read quality

                            base_number+=1 # count number of SNPs included in data
                            for n in range(NSAMS):

                                # retrieve bases for this particular sample
                                n=n+1
                                myReads = Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                                (bases, indexDelN) = generics.convertSyms(myReads,mySite)
                                myReads = Reads(bases, myReads.base_quality)

                                # filter by quality
                                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                                myReads=Reads(bases,qualities)
                                # find all indexes of occurances to be filtered out
                                index_of_X=[]
                                index=-1
                                while True:
                                        index=myReads.base.find('X',index+1)
                                        if index == -1:
                                            break  # all occurrences have been found
                                        index_of_X.append(index)
                                myReads.base=myReads.base.replace('X','')
                                # remove all corresponding base qualities
                                count=0
                                for i in index_of_X:
                                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                                    count+=1
                                major_count=generics.calcAlleleFreq(major,myReads)
                                minor_count=generics.calcAlleleFreq(minor,myReads)
                                # take a sample of the bases so that the proportion of data used is as required
                                if downsampling<1:
                                    data_prop = math.ceil(len(myReads.base)*downsampling) # calculate how many bases to include for proportion of sample
                                    rand_samp = random.sample(range(0,len(myReads.base)),data_prop)
                                    base = ""
                                    qualities = ""
                                    for r in rand_samp:
                                        base+=myReads.base[r]
                                        qualities+=myReads.base_quality[r]
                                    myReads=Reads(base,qualities)
                                # find sample depth of filtered data
                                sampleDepth = len(myReads.base)
                                NUMSITES[0]+=sampleDepth # count the number of reads for each sample
                                NUMSITES[n]+=sampleDepth
                                sep="\t"
                                content=""
                                content=(mySite.chrom,str(mySite.position),str(n),mySite.reference,str(sampleDepth),alleles[major],alleles[minor],str(major_count),str(minor_count))
                                content=sep.join(content)
                                content+="\t"

                                for ip in ploidy:
                                    Nploid = generics.calcGenoLogLikeN_MajorMinor(ip, myReads, mySite, major, minor)
                                    content2 ="\t".join(map(str,Nploid))
                                    content += content2
                                    content += "\t"
                                content += "\n"

                                # write file of genotype likelihoods
                                with gzip.open(output,'at+') as f:
                                    f.write(content)


                                # end likelihood calc
                            # end for sample
                        # end for if max minor
                    # end if not filtered for global depth
                # end for line

        try:
            os.remove(mpFile) # remove unnecessary mpileup file
        except:
            pass
    
elif fileType == 2: # mpilup file
    for g1 in list_of_inputs: # for every filename 
        directory = '/'.join(g1.split('/')[:-1])
        if len(directory) == 0:
            g = "./" + g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = "./" + g2 + ".genolikes.gz"
            else:
                output = outFolder+'/'.join(g2.split('/')[-1])+".genolikes.gz"
        else:
            g = g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = g2 + ".genolikes.gz"
            else:
                output = outFolder+"/"+g2+".genolikes.gz"
        print("Output file is: " + output)
        with open(g) as mp:
            first_line = mp.readline()
            Data = first_line.strip('\n')
            l = Data.split('\t')
            NSAMS=int((len(l)-3)/3)
            if args.Inbreeding:
                inbreed = args.Inbreeding
            else:
                inbreed = "0x{}".format(str(NSAMS))
            # parse inbreeding coeffients
            F=[]
            temp=inbreed.split(',')
            for t in temp:
                vals=t.split('x')
                if len(vals)==2:
                    F+=list(np.repeat(float(vals[0]),int(vals[1])))
                else:
                    F+=[vals[0]]
            print(F)
            downsampling=float(args.downsampling) # fraction of data to be used (0-1]. 
            #Original_sample_number=NSAMS
            win=50 # window size for calculating ploidy
            phredscale=33
            NUMSITES=np.zeros(NSAMS+1,int)
            NUMSITES_HWE=np.zeros(NSAMS+1,int)
            ExpectedPloidy=[[] for i in range(NSAMS+1)]
            Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float) # (NSAMS+1)xploidies array for probabilities of each ploidy for each sample and overall ploidy probabilities
            Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
            delta_prob=np.zeros((NSAMS+1,len(ploidy)),float)
            counts=np.zeros((NSAMS+1,len(ploidy)),float)+1 # counts of bases for each ploidy being most likely
            base_number=0 # count for bases
            list_of_window=[]
            list_of_window2=[]
            no_bases=0
            total_bases=0
            
            for line in mp:
                Data=line.strip('\n')# convert bytes into string
                l = Data.split('\t')
                mySite = Site(str(l[0]),int(l[1]),str(l[2]))
                myReads = Reads("","")
                # pooled reads for first level filtering (global depth) and estimation of minor/major alleles and allele frequencies
                individualDepth = np.zeros(NSAMS,float)
                for n in range(NSAMS):
                    n=n+1
                    subReads= Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                    individualDepth[n-1] = len(subReads.base)
                    myReads.base=str(myReads.base+subReads.base)
                    myReads.base_quality=str(myReads.base_quality+subReads.base_quality)

                # convert to bases
                [bases,indexDelN] = generics.convertSyms(myReads,mySite)
                myReads=Reads(bases,myReads.base_quality)
                if len(myReads.base)!=len(myReads.base_quality):
                    sys.exit("Conversion not succesful")
                # filter by quality
                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                myReads=Reads(bases,qualities)

                # find all indexes of occurances to be filtered out
                index_of_X=[]
                index=-1
                while True:
                        index=myReads.base.find('X',index+1)
                        if index == -1:
                            break  # all occurrences have been found
                        index_of_X.append(index)
                myReads.base=myReads.base.replace('X','')
                # remove all corresponding base qualities
                count=0
                for i in index_of_X:
                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                    count+=1
                globalDepth = len(myReads.base)
                if ((globalDepth > args.min_global_depth) & (min(individualDepth)>args.min_ind_depth)):
                    total_bases+=globalDepth
                    no_bases+=1
                    # counts of non-major bases
                    nonMajorCount = generics.calcNonMajorCounts(myReads)
                    nonMajorProp = nonMajorCount/len(myReads.base)
                    # filter the site based on global depth
                    Set_min_prop = args.min_non_major_freq # minimum proportion of nonMajorCount
                    if nonMajorProp>Set_min_prop: # remove bases where more that 1-Set_min_prop are major allele i.e monomorphic bases
                        prob_of_ancestral_allelle_maj=1-nonMajorProp

                        haploid = generics.calcGenoLogLike1(myReads,mySite)
                        tri_ref = haploid[4] # retrieve reference value for if the base is not triallilic
                        haploid=haploid[:4] # remove reference value
                        # keep reference allele as one possible allele so always assume KeepRef=0
                        [major,minor,minor2,minor3] = [haploid.index(sorted(haploid,reverse=True)[0]),haploid.index(sorted(haploid,reverse=True)[1]),haploid.index(sorted(haploid,reverse=True)[2]),haploid.index(sorted(haploid,reverse=True)[3])]
                        # remove sites with >0.1 frequency of minor 2 or minor 3 allele to remove non biallilic sites (0.1 error built in for sequencing error)
                        minor2_prop=generics.calcAlleleFreq(minor2,myReads)/len(myReads.base) # Calculate allele frequencies of minor2&3 alleles
                        minor3_prop=generics.calcAlleleFreq(minor3,myReads)/len(myReads.base)

                        if(minor2_prop<args.max_minor2_freq and minor3_prop<args.max_minor3_freq):
                            P_bar=0
                            Q_bar=0
                            for read in range(len(myReads.base)):
                                if myReads.base[read]==alleles[major]:
                                    P_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))
                                elif myReads.base[read]==alleles[minor]:
                                    Q_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))

                            P = P_bar/(P_bar+Q_bar) # proportion of major allele weigted by read quality
                            Q = Q_bar/(P_bar+Q_bar) # proportion of minor allele weigted by read quality

                            base_number+=1 # count number of SNPs included in data
                            for n in range(NSAMS):

                                # retrieve bases for this particular sample
                                n=n+1
                                myReads = Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                                (bases, indexDelN) = generics.convertSyms(myReads,mySite)
                                myReads = Reads(bases, myReads.base_quality)

                                # filter by quality
                                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                                myReads=Reads(bases,qualities)
                                # find all indexes of occurances to be filtered out
                                index_of_X=[]
                                index=-1
                                while True:
                                        index=myReads.base.find('X',index+1)
                                        if index == -1:
                                            break  # all occurrences have been found
                                        index_of_X.append(index)
                                myReads.base=myReads.base.replace('X','')
                                # remove all corresponding base qualities
                                count=0
                                for i in index_of_X:
                                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                                    count+=1
                                major_count=generics.calcAlleleFreq(major,myReads)
                                minor_count=generics.calcAlleleFreq(minor,myReads)
                                # take a sample of the bases so that the proportion of data used is as required
                                if downsampling<1:
                                    data_prop = math.ceil(len(myReads.base)*downsampling) # calculate how many bases to include for proportion of sample
                                    rand_samp = random.sample(range(0,len(myReads.base)),data_prop)
                                    base = ""
                                    qualities = ""
                                    for r in rand_samp:
                                        base+=myReads.base[r]
                                        qualities+=myReads.base_quality[r]
                                    myReads=Reads(base,qualities)
                                # find sample depth of filtered data
                                sampleDepth = len(myReads.base)
                                NUMSITES[0]+=sampleDepth # count the number of reads for each sample
                                NUMSITES[n]+=sampleDepth
                                sep="\t"
                                content=""
                                content=(mySite.chrom,str(mySite.position),str(n),mySite.reference,str(sampleDepth),alleles[major],alleles[minor],str(major_count),str(minor_count))
                                content=sep.join(content)
                                content+="\t"

                                for ip in ploidy:
                                    Nploid = generics.calcGenoLogLikeN_MajorMinor(ip, myReads, mySite, major, minor)
                                    content2 ="\t".join(map(str,Nploid))
                                    content += content2
                                    content += "\t"
                                content += "\n"

                                # write file of genotype likelihoods
                                with gzip.open(output,'at+') as f:
                                    f.write(content)

                                # end likelihood calc
                            # end for sample
                        # end for if max minor
                    # end if not filtered for global depth
                # end for line

elif fileType == 3: # mpileup.gz file
    for g1 in list_of_inputs: # for every filename 
        directory = '/'.join(g1.split('/')[:-1])
        if len(directory) == 0:
            g = "./" + g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = "./" + g2 + ".genolikes.gz"
            else:
                output = outFolder+'/'.join(g2.split('/')[-1])+".genolikes.gz"
        else:
            g = g1
            g2 = g1[:-exl]
            if outFolder==0:
                output = g2 + ".genolikes.gz"
            else:
                output = outFolder+"/"+g2+".genolikes.gz"
        print("Output file is: " + output)
        with gzip.open(g,'rb') as gz:# opens the mpilup.gz file. Use mpileup.read() to display content
            first_line = gz.readline()
            Data=first_line.decode().strip('\n') # Convert bytes into string
            l = Data.split('\t')
            NSAMS=int((len(l)-3)/3)
            if args.Inbreeding:
                inbreed = args.Inbreeding
            else:
                inbreed = "0x{}".format(str(NSAMS))
            # parse inbreeding coeffients
            F=[]
            temp=inbreed.split(',')
            for t in temp:
                vals=t.split('x')
                if len(vals)==2:
                    F+=list(np.repeat(float(vals[0]),int(vals[1])))
                else:
                    F+=[vals[0]]
            print(F)
            downsampling=float(args.downsampling) # fraction of data to be used (0-1]. 
            #Original_sample_number=NSAMS
            win=50 # window size for calculating ploidy
            phredscale=33
            NUMSITES=np.zeros(NSAMS+1,int)
            NUMSITES_HWE=np.zeros(NSAMS+1,int)
            ExpectedPloidy=[[] for i in range(NSAMS+1)]
            Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float) # (NSAMS+1)xploidies array for probabilities of each ploidy for each sample and overall ploidy probabilities
            Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
            delta_prob=np.zeros((NSAMS+1,len(ploidy)),float)
            counts=np.zeros((NSAMS+1,len(ploidy)),float)+1 # counts of bases for each ploidy being most likely
            base_number=0 # count for bases
            list_of_window=[]
            list_of_window2=[]
            no_bases=0
            total_bases=0
            for line in gz:
                Data=line.decode().strip('\n')# convert bytes into string
                l = Data.split('\t')
                mySite = Site(str(l[0]),int(l[1]),str(l[2]))
                myReads = Reads("","")

                # pooled reads for first level filtering (global depth) and estimation of minor/major alleles and allele frequencies
                individualDepth = np.zeros(NSAMS,float)
                for n in range(NSAMS):
                    n=n+1
                    subReads= Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                    individualDepth[n-1] = len(subReads.base)
                    myReads.base=str(myReads.base+subReads.base)
                    myReads.base_quality=str(myReads.base_quality+subReads.base_quality)

                # convert to bases
                [bases,indexDelN] = generics.convertSyms(myReads,mySite)
                myReads=Reads(bases,myReads.base_quality)
                if len(myReads.base)!=len(myReads.base_quality):
                    sys.exit("Conversion not succesful")
                # filter by quality
                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                myReads=Reads(bases,qualities)

                # find all indexes of occurances to be filtered out
                index_of_X=[]
                index=-1
                while True:
                        index=myReads.base.find('X',index+1)
                        if index == -1:
                            break  # all occurrences have been found
                        index_of_X.append(index)
                myReads.base=myReads.base.replace('X','')
                # remove all corresponding base qualities
                count=0
                for i in index_of_X:
                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                    count+=1

                globalDepth = len(myReads.base)
                if ((globalDepth > args.min_global_depth) & (min(individualDepth)>args.min_ind_depth)):
                    total_bases+=globalDepth
                    no_bases+=1
                    # counts of non-major bases
                    nonMajorCount = generics.calcNonMajorCounts(myReads)
                    nonMajorProp = nonMajorCount/len(myReads.base)
                    # filter the site based on global depth
                    Set_min_prop = args.min_non_major_freq # minimum proportion of nonMajorCount
                    if nonMajorProp>Set_min_prop: # remove bases where more that 1-Set_min_prop are major allele i.e monomorphic bases
                        prob_of_ancestral_allelle_maj=1-nonMajorProp
                        haploid = generics.calcGenoLogLike1(myReads,mySite)
                        tri_ref = haploid[4] # retrieve reference value for if the base is not triallilic
                        haploid=haploid[:4] # remove reference value
                        # keep reference allele as one possible allele so always assume KeepRef=0
                        [major,minor,minor2,minor3] = [haploid.index(sorted(haploid,reverse=True)[0]),haploid.index(sorted(haploid,reverse=True)[1]),haploid.index(sorted(haploid,reverse=True)[2]),haploid.index(sorted(haploid,reverse=True)[3])]
                        # remove sites with >0.1 frequency of minor 2 or minor 3 allele to remove non biallilic sites (0.1 error built in for sequencing error)
                        minor2_prop=generics.calcAlleleFreq(minor2,myReads)/len(myReads.base) # Calculate allele frequencies of minor2&3 alleles
                        minor3_prop=generics.calcAlleleFreq(minor3,myReads)/len(myReads.base)

                        if(minor2_prop<args.max_minor2_freq and minor3_prop<args.max_minor3_freq):
                            P_bar=0
                            Q_bar=0
                            for read in range(len(myReads.base)):
                                if myReads.base[read]==alleles[major]:
                                    P_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))
                                elif myReads.base[read]==alleles[minor]:
                                    Q_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))

                            P = P_bar/(P_bar+Q_bar) # proportion of major allele weigted by read quality
                            Q = Q_bar/(P_bar+Q_bar) # proportion of minor allele weigted by read quality

                            base_number+=1 # count number of SNPs included in data
                            for n in range(NSAMS):

                                # retrieve bases for this particular sample
                                n=n+1
                                myReads = Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                                (bases, indexDelN) = generics.convertSyms(myReads,mySite)
                                myReads = Reads(bases, myReads.base_quality)

                                # filter by quality
                                [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                                myReads=Reads(bases,qualities)
                                # find all indexes of occurances to be filtered out
                                index_of_X=[]
                                index=-1
                                while True:
                                        index=myReads.base.find('X',index+1)
                                        if index == -1:
                                            break  # all occurrences have been found
                                        index_of_X.append(index)
                                myReads.base=myReads.base.replace('X','')
                                # remove all corresponding base qualities
                                count=0
                                for i in index_of_X:
                                    myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                                    count+=1
                                major_count=generics.calcAlleleFreq(major,myReads)
                                minor_count=generics.calcAlleleFreq(minor,myReads)
                                # take a sample of the bases so that the proportion of data used is as required
                                if downsampling<1:
                                    data_prop = math.ceil(len(myReads.base)*downsampling) # calculate how many bases to include for proportion of sample
                                    rand_samp = random.sample(range(0,len(myReads.base)),data_prop)
                                    base = ""
                                    qualities = ""
                                    for r in rand_samp:
                                        base+=myReads.base[r]
                                        qualities+=myReads.base_quality[r]
                                    myReads=Reads(base,qualities)
                                # find sample depth of filtered data
                                sampleDepth = len(myReads.base)
                                NUMSITES[0]+=sampleDepth # count the number of reads for each sample
                                NUMSITES[n]+=sampleDepth
                                sep="\t"
                                content=""
                                content=(mySite.chrom,str(mySite.position),str(n),mySite.reference,str(sampleDepth),alleles[major],alleles[minor],str(major_count),str(minor_count))
                                content=sep.join(content)
                                content+="\t"

                                for ip in ploidy:
                                    Nploid = generics.calcGenoLogLikeN_MajorMinor(ip, myReads, mySite, major, minor)
                                    content2 ="\t".join(map(str,Nploid))
                                    content += content2
                                    content += "\t"
                                content += "\n"

                                # write file of genotype likelihoods
                                with gzip.open(output,'at+') as f:
                                    f.write(content)

                                # end likelihood calc
                            # end for sample
                        # end for if max minor
                    # end if not filtered for global depth
                # end for line
