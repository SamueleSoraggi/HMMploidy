#! usr/bin/python3
import sys
import gzip
import generics
import numpy as np
import math 
import scipy.stats
import random
from statistics import mode

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

ploidy = [1,2,3,4,5,6,7,8]

input = sys.argv[1] #Input file in form of gziped mpileup
output = sys.argv[2] # output file to contain genotype likelihoods
NSAMS = sys.argv[3]
NSAMS = int(NSAMS)
F=list(np.zeros(NSAMS)) # Inbreeding coefficient for each sample currently all taken as 0. Update in later version
phredscale=33
gzip.open(output,'w+') 
no_bases=0
total_bases=0
with gzip.open(input,'rb') as gz:# opens the mpilup. Use mpileup.read() to display content
    for line in gz:
        Data=line.decode().strip('\n')# Convert bytes into string
        l = Data.split('\t')
        #if((NSAMS)!=int((len(l)-3)/3)):
        #    sys.exit("Number of smaples does not match the mpileup (%s)"%((len(l)-3)/3))

        mySite = Site(str(l[0]),int(l[1]),str(l[2]))
        myReads = Reads("","")

        # pooled reads for first level filtering (global depth) and estimation of minor/major alleles and allele frequencies
        for n in range(NSAMS):
            n=n+1
            subReads= Reads(l[(n-1)*3+4],l[(n-1)*3+5])
            myReads.base=str(myReads.base+subReads.base)
            myReads.base_quality=str(myReads.base_quality+subReads.base_quality)

        # convert to bases
        [bases,indexDelN] = generics.convertSyms(myReads,mySite)
        myReads=Reads(bases,myReads.base_quality)
        if len(myReads.base)!=len(myReads.base_quality):
            sys.exit("Conversion not succesful")
        #filter by quality
        #currently just removing X                        
        #find all indexes of occurances to be filtered out
        index_of_X=[]
        index=-1
        while True:
                index=myReads.base.find('X',index+1)
                if index == -1:
                    break  # All occurrences have been found
                index_of_X.append(index)
        myReads.base=myReads.base.replace('X','')
        #Remove all corresponding base qualities
        count=0
        for i in index_of_X:
            myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
            count+=1




        
        globalDepth = len(myReads.base)
        if globalDepth > 0:     
            total_bases+=globalDepth
            no_bases+=1
            #counts of non-major bases
            nonMajorCount = generics.calcNonMajorCounts(myReads)
            nonMajorProp = nonMajorCount/len(myReads.base)
            #filter the site based on global depth
            Set_min_prop = 0.2 # minimum proportion of nonMajorCount
            if nonMajorProp>Set_min_prop: # remove bases where more that 1-Set_min_prop are major allele i.e monomorphic bases
                prob_of_ancestral_allelle_maj=1-nonMajorProp

                haploid = generics.calcGenoLogLike1(myReads,mySite)
                tri_ref = haploid[4] # retrieve reference value for if the base is not triallilic
                haploid=haploid[:4] # remove reference value
                # Keep reference allele as one possible allele so always assume KeepRef=0
                [major,minor,minor2,minor3] = [haploid.index(sorted(haploid,reverse=True)[0]),haploid.index(sorted(haploid,reverse=True)[1]),haploid.index(sorted(haploid,reverse=True)[2]),haploid.index(sorted(haploid,reverse=True)[3])] 
                #remove sites with >0.1 frequency of minor 2 or minor 3 allele to remove non biallilic sites (0.1 error built in for sequencing error)
                minor2_prop=generics.calcAlleleFreq(minor2,myReads)/len(myReads.base) # Calculate allele frequencies of minor2&3 alleles
                minor3_prop=generics.calcAlleleFreq(minor3,myReads)/len(myReads.base)


                if(minor2_prop<0.1 and minor3_prop<0.1):
                    P_bar=0
                    Q_bar=0
                    for read in range(len(myReads.base)):
                        if myReads.base[read]==alleles[major]:
                            P_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))
                        elif myReads.base[read]==alleles[minor]:
                            Q_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))

                    
                    P = P_bar/(P_bar+Q_bar) #proportion of major allele weigted by read quality
                    Q = Q_bar/(P_bar+Q_bar) #proportion of minor allele weigted by read quality

                    
                    for n in range(NSAMS):
                        
                        #retrieve bases for this particular sample
                        n=n+1
                        myReads = Reads(l[(n-1)*3+4],l[(n-1)*3+5])
                        (bases, indexDelN) = generics.convertSyms(myReads,mySite)
                        myReads = Reads(bases, myReads.base_quality)
                        #establish prior probabilities
                        HWE_Prob_hap = [P,Q]
                        HWE_Prob_dip = [(1-F[n-1])*(P**2)+F[n-1]*P,(1-F[n-1])*2*P*Q,(1-F[n-1])*(Q**2)+F[n-1]*Q]
                        HWE_Prob_tri = [(1-F[n-1])*(P**3)+F[n-1]*P,(1-F[n-1])*3*(P**2)*Q,(1-F[n-1])*3*P*(Q**2),(1-F[n-1])*(Q**3)+F[n-1]*Q]
                        HWE_Prob_tetra = [(1-F[n-1])*(P**4)+F[n-1]*P,(1-F[n-1])*4*(P**3)*Q,(1-F[n-1])*6*(P**2)*(Q**2),(1-F[n-1])*4*P*(Q**3),(1-F[n-1])*(Q**4)+F[n-1]*Q]
                        HWE_Prob_pent = [(1-F[n-1])*(P**5)+F[n-1]*P,(1-F[n-1])*5*(P**4)*Q,(1-F[n-1])*10*(P**3)*(Q**2),(1-F[n-1])*10*(P**2)*(Q**3),(1-F[n-1])*5*P*(Q**4),(1-F[n-1])*(Q**5)+F[n-1]*Q]
                        HWE_Prob_hex = [(1-F[n-1])*(P**6)+F[n-1]*P,(1-F[n-1])*6*(P**5)*Q,(1-F[n-1])*15*(P**4)*(Q**2),(1-F[n-1])*20*(P**3)*(Q**3),(1-F[n-1])*15*(P**2)*(Q**4),(1-F[n-1])*6*P*(Q**5),(1-F[n-1])*(Q**6)+F[n-1]*Q]
                        HWE_Prob_hept = [(1-F[n-1])*(P**7)+F[n-1]*P,(1-F[n-1])*7*(P**6)*Q,(1-F[n-1])*21*(P**5)*(Q**2),(1-F[n-1])*35*(P**4)*(Q**3),(1-F[n-1])*35*(P**3)*(Q**4),(1-F[n-1])*21*(P**2)*(Q**5),(1-F[n-1])*7*P*(Q**6),(1-F[n-1])*(Q**7)+F[n-1]*Q]
                        HWE_Prob_oct = [(1-F[n-1])*(P**8)+F[n-1]*P,(1-F[n-1])*8*(P**7)*Q,(1-F[n-1])*28*(P**6)*(Q**2),(1-F[n-1])*56*(P**5)*(Q**3),(1-F[n-1])*70*(P**4)*(Q**4),(1-F[n-1])*56*(P**3)*(Q**5),(1-F[n-1])*28*(P**2)*(Q**6),(1-F[n-1])*8*P*(Q**7),(1-F[n-1])*(Q**8)+F[n-1]*Q]

                        #filter by quality
                        #find all indexes of occurances to be filtered out
                        index_of_X=[]
                        index=-1
                        while True:
                                index=myReads.base.find('X',index+1)
                                if index == -1:
                                    break  # All occurrences have been found
                                index_of_X.append(index)
                        myReads.base=myReads.base.replace('X','')
                        #Remove all corresponding base qualities
                        count=0
                        for i in index_of_X:
                            myReads.base_quality = myReads.base_quality[:i-count] + myReads.base_quality[i+1-count:]
                            count+=1
                        #find sample depth of filtered data    
                        sampleDepth = len(myReads.base)


                        

                        if 1 in ploidy:
                            haploid = generics.calcGenoLogLike1_MajorMinor(myReads,mySite,major,minor)
                            
                        if 2 in ploidy:
                            diploid = generics.calcGenoLogLike2_MajorMinor(myReads,mySite,major,minor)


                        if 3 in ploidy:
                            triploid = generics.calcGenoLogLike3_MajorMinor(myReads,mySite,major,minor)


                        if 4 in ploidy:
                            tetraploid = generics.calcGenoLogLike4_MajorMinor(myReads,mySite,major,minor)


                        if 5 in ploidy:
                            pentaploid = generics.calcGenoLogLike5_MajorMinor(myReads,mySite,major,minor)


                        if 6 in ploidy:
                            hexaploid = generics.calcGenoLogLike6_MajorMinor(myReads,mySite,major,minor)

                        if 7 in ploidy:
                            heptaploid = generics.calcGenoLogLike7_MajorMinor(myReads,mySite,major,minor)


                        if 8 in ploidy:
                            octaploid = generics.calcGenoLogLike8_MajorMinor(myReads,mySite,major,minor)



                        # Write file of genotype likelihoods
                        sep="\t"
                        content=(mySite.chrom,str(mySite.position),str(n),mySite.reference,str(sampleDepth),alleles[major],alleles[minor],"\t".join(map(str,haploid)),"\t".join(map(str,diploid)),"\t".join(map(str,triploid)),"\t".join(map(str,tetraploid)),"\t".join(map(str,pentaploid)),"\t".join(map(str,hexaploid)),"\t".join(map(str,heptaploid)),"\t".join(map(str,octaploid)))
                        content=sep.join(content)
                        content=content+"\n"
                        with gzip.open(output,'at+') as f: 
                            f.write(content)
                        content=""




                        #end likelihood calc
                    # end for sample
            #end if not filtered for global depth
        #end for line   





















