#! usr/bin/python3

import math

phredScale=33

def convertSyms(read,site):
    ''' read must be a Reads type object and site a Sites type object'''
    bases=""
    indexDelN=[]

    i=0
    while i <= len(read.base)-1:
        if read.base[i] in ['.',',']: #reference
            bases=str(bases+site.reference)
        elif read.base[i] in ['A','C','G','T','a','c','g','t']: #alternate
            bases=str(bases+read.base[i].upper())
        elif read.base[i] in ['^']: # start read, index +1 since the following character is mapping quality
            i=i+1
        elif read.base[i] in ['$']: #end read
            i=i #do nothing but keep for clarity
        elif read.base[i] in ['*','N','n','>','<']:  # asterisk is deleted base, N or n is undefined, > and < are reference skips, these will be then filtered out later on
            bases=str(bases+'X')
            indexDelN.append(i)
        elif read.base[i] in ['-','+']: # indel, skip to the the next non-indel base
            lenIndel = int(read.base[i+1])
            try: # Indels longer than 9
                lenIndel=int(str(lenIndel)+str(int(read.base[i+2])))+1
            except:
                pass
            i=i+lenIndel+1
        i=i+1
    return(bases,indexDelN)

def calcNonMajorCounts(read):

    alleles =['A','C','G','T']
    counts = [0,0,0,0]

    if len(read.base)>0:

        for i in range(len(read.base)):
            counts[alleles.index(read.base[i])]+=1
    return sum(counts)-max(counts)

def calcAlleleFreq(Allele,Reads):
    alleles =['A','C','G','T']
    Allele=alleles[Allele]
    tot=0
    if len(Reads.base)>0:
        for i in range(len(Reads.base)):
            if Reads.base[i]==Allele:
                tot+=1
    return tot

def calcGenoLogLike1(reads,site):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0]
    phredScale=33

    #cycle across all possible genotypes
    for j in range(len(alleles)):
        if j == 0:                
            for i in range(len(reads.base)):
            
         
                #get base probability from quality score
                bP = 10**((phredScale-ord(str(reads.base_quality[i])))/10)

                sublike=0.0
            
                if alleles[j]==reads.base[i]:
                    sublike += 1-(bP)
                else:
                    sublike += (bP/3)

                log_likes[j] += math.log(sublike)
                log_likes[4] += math.log(bP/3)
        else:
            for i in range(len(reads.base)):
                #get base probability from quality score
                bP = 10**((phredScale-ord(str(reads.base_quality[i])))/10)
                sublike=0.0
                if alleles[j]==reads.base[i]:
                    sublike += 1-(bP)
                else:
                    sublike += (bP/3)
                log_likes[j] += math.log(sublike)
    return log_likes


# calculate genotype likelihoods (in ln format) in case of haploids for Major and Minor
def calcGenoLogLike1_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0]
    iter=-1
    ploidy=1
    phredScale=33

    # cycle across all possible genotypes
    for j in [major,minor]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            sublike = 0.0
            if alleles[j] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)

# calculate genotype likelihoods (in ln format) in case of diploids for Major and Minor
def calcGenoLogLike2_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0]
    iter=-1
    ploidy=2
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2] in [[major,major],[major,minor],[minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)


# calculate genotype likelihoods (in ln format) in case of triploids for Major and Minor
def calcGenoLogLike3_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=3
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3] in [[major,major,major],[major,major,minor],[major,minor,minor],[minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)


# calculate genotype likelihoods (in ln format) in case of tetraploids for Major and Minor
def calcGenoLogLike4_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=4
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3,j4] in [[major,major,major,major],[major,major,major,minor],[major,major,minor,minor],[major,minor,minor,minor],[minor,minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j4] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)

# calculate genotype likelihoods (in ln format) in case of pentaploids for Major and Minor
def calcGenoLogLike5_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=5
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3,j4,j5] in [[major,major,major,major,major],[major,major,major,major,minor],[major,major,major,minor,minor],[major,major,minor,minor,minor],[major,minor,minor,minor,minor],[minor,minor,minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j4] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
        
            if alleles[j5] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)


# calculate genotype likelihoods (in ln format) in case of sessaploids for Major and Minor
def calcGenoLogLike6_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=6
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3,j4,j5,j6] in [[major,major,major,major,major,major],[major,major,major,major,major,minor],[major,major,major,major,minor,minor],[major,major,major,minor,minor,minor],[major,major,minor,minor,minor,minor],[major,minor,minor,minor,minor,minor],[minor,minor,minor,minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j4] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
        
            if alleles[j5] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j6] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)


# calculate genotype likelihoods (in ln format) in case of settaploids for Major and Minor
def calcGenoLogLike7_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=7
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3,j4,j5,j6,j7] in [[major,major,major,major,major,major,major],[major,major,major,major,major,major,minor],[major,major,major,major,major,minor,minor],[major,major,major,major,minor,minor,minor],[major,major,major,minor,minor,minor,minor],[major,major,minor,minor,minor,minor,minor],[major,minor,minor,minor,minor,minor,minor],[minor,minor,minor,minor,minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j4] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
        
            if alleles[j5] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j6] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j7] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)


# calculate genotype likelihoods (in ln format) in case of octaploids for Major and Minor
def calcGenoLogLike8_MajorMinor(read,site,major,minor):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    iter=-1
    ploidy=8
    phredScale=33

    # cycle across all possible genotypes
    for [j1,j2,j3,j4,j5,j6,j7,j8] in [[major,major,major,major,major,major,major,major],[major,major,major,major,major,major,major,minor],[major,major,major,major,major,major,minor,minor],[major,major,major,major,major,minor,minor,minor],[major,major,major,major,minor,minor,minor,minor],[major,major,major,minor,minor,minor,minor,minor],[major,major,minor,minor,minor,minor,minor,minor],[major,minor,minor,minor,minor,minor,minor,minor],[minor,minor,minor,minor,minor,minor,minor,minor]]:
        iter += 1

        for i in range(len(read.base)):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            
            sublike = 0.0
            if alleles[j1] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j2] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j3] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j4] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
        
            if alleles[j5] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j6] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j7] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy

            if alleles[j8] == read.base[i]:
                sublike += (1-bP)/ploidy
            else:
                sublike += (bP/3)/ploidy
            
            log_likes[iter] += math.log(sublike)
            
    return(log_likes)

def exp_or_zero(x):
    if(x==0):
        x=0.0
    else:
        x=math.exp(x)
    return(x)

def log_or_zero(x):
    if(x==0):
        x=0.0
    else:
        x=math.log(x)
    return(x)

def delta_to_ploidy(delta_prob):
    output=0
    check=0
    for i in range(len(delta_prob)-1): 
        if math.exp(delta_prob[i])-math.exp(delta_prob[i+1])<0: # look for a case of where the jump in delta_prob has increased 
            output=i
            check=1
            break
    if check==0:
        output=1
    else:
        output=output
    return(output)

def dist(ploidies):
    ''' Function to return the distribution of ploidies predicted from the inputted array'''
    number = len(ploidies)
    dist = [ploidies.count(i)/number for i in range(1,9)]
    return dist


def filter(reads,min_quality_score):
    phredScale=33
    bases=""
    qualities=""
    reads_count=len(reads.base)
    for i in range(0,reads_count):
        if ord(reads.base_quality[i])-phredScale>min_quality_score:
            bases+=reads.base[i]
            qualities+=reads.base_quality[i]
    return(bases,qualities)

