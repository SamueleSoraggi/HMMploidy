#! usr/bin/python3
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

ploidy = [1,2,3,4,5,6,7,8]



parser = argparse.ArgumentParser()
parser.add_argument("input",help="file containing the list of basenames for gzipped mpileup files for use in analysis to be used")
parser.add_argument("-i","--Inbreeding",help="Inbreeding coefficients for samples e.g 0.1x3,0.2 = 0.1,0.1,0.1,0.2 ")
parser.add_argument("-d","--downsampling",help="Fraction of data to be used in the calculations",default=1)
parser.add_argument("-m","--min_non_major_freq",type=float,help="Set the minimum frequency of non major alleles for bases to be included in the calculations",default=0.2)
parser.add_argument("-q","--min_quality_score",type=int,help="Set the minimum quality score of a read to be included in the calculation",default=1)
parser.add_argument("-dp","--min_global_depth",type=float,help="Set the minimum global depth of a base to be included in calculations",default=0)
parser.add_argument("-M2","--max_minor2_freq",type=float,help="Set the maximum frequency of third most prolific alleles for bases to be included in the calculations",default=0.1)
parser.add_argument("-M3","--max_minor3_freq",type=float,help="Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations",default=0.1)
args = parser.parse_args()




input = args.input #Input file in form of gziped mpileup
list_of_inputs=[]
with open(input,'rb') as f:# opens the mpilup. Use mpileup.read() to display content
    for line in f:
        line=line.decode().strip('\n')# Convert bytes into string
        #line=line+".mpileup.gz"
        list_of_inputs.append(line)
Nfiles=len(list_of_inputs)
print('%d files found' %Nfiles)
print(list_of_inputs)

for g in list_of_inputs: #output files names
    output = g+".genolikes.gz" # output file to contain genotype likelihoods
    #output_2 = g+".ploids" # output file to contain inferred ploidies and whether aneuploidy or not

for g in list_of_inputs: #formatting input files names
    directory = g.split('/')
    if len(directory) == 1:
        list_of_inputs.append("./"+line+".mpileup.gz")
    else:
        list_of_inputs.append(line+".mpileup.gz")

for g in list_of_inputs:        
    with gzip.open(g) as f:
        first_line = f.readline()
        Data=first_line.decode().strip('\n')# Convert bytes into string
        l = Data.split('\t')
    NSAMS=int((len(l)-3)/3)
    if args.Inbreeding:
        inbreed = args.Inbreeding
    else:
        inbreed = "0x{}".format(str(NSAMS))


    #parse inbreeding coeffients
    F=[]
    temp=inbreed.split(',')
    for t in temp:
        vals=t.split('x')
        F+=list(np.repeat(float(vals[0]),int(vals[1])))
    downsampling=float(args.downsampling) # fraction of data to be used (0-1]. 
    Original_sample_number=NSAMS
    win=50 # window size for calculating ploidy
    phredscale=33
    NUMSITES=np.zeros(NSAMS+1,int)
    NUMSITES_HWE=np.zeros(NSAMS+1,int)
    ExpectedPloidy=[[] for i in range(NSAMS+1)]
    Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float) # (NSAMS+1)xploidies array for probabilities of each ploidy for each sample and overall ploidy probabilities
    Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
    delta_prob=np.zeros((NSAMS+1,len(ploidy)),float)
    counts=np.zeros((NSAMS+1,len(ploidy)),float)+1 # Counts of bases for each ploidy being most likely 
    base_number=0 # count for bases
    list_of_window=[] 
    list_of_window2=[]
    gzip.open(directory+'/'+output,'wt') 
    gzip.open(directory+'/'+output_2,'wt')
    no_bases=0
    total_bases=0
    with gzip.open(directory+'/'+g,'rb') as gz:# opens the mpilup. Use mpileup.read() to display content
        for line in gz:
            Data=line.decode().strip('\n')# Convert bytes into string
            l = Data.split('\t')
            #if((NSAMS)!=int((len(l)-3)/3)):
            #    sys.exit("Number of samples does not match the mpileup (%s)"%((len(l)-3)/3))

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
            [bases,qualities] = generics.filter(myReads,args.min_quality_score)
            myReads=Reads(bases,qualities)

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
            if globalDepth > args.min_global_depth:     
                total_bases+=globalDepth
                no_bases+=1
                #counts of non-major bases
                nonMajorCount = generics.calcNonMajorCounts(myReads)
                nonMajorProp = nonMajorCount/len(myReads.base)
                #filter the site based on global depth
                Set_min_prop = args.min_non_major_freq # minimum proportion of nonMajorCount
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


                    if(minor2_prop<args.max_minor2_freq and minor3_prop<args.max_minor3_freq):
                        P_bar=0
                        Q_bar=0
                        for read in range(len(myReads.base)):
                            if myReads.base[read]==alleles[major]:
                                P_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))
                            elif myReads.base[read]==alleles[minor]:
                                Q_bar+=(1-(10**((phredscale-ord(str(myReads.base_quality[read])))/10)))


                        P = P_bar/(P_bar+Q_bar) #proportion of major allele weigted by read quality
                        Q = Q_bar/(P_bar+Q_bar) #proportion of minor allele weigted by read quality


                        base_number+=1 # Count number of SNPs included in data

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
                            [bases,qualities] = generics.filter(myReads,args.min_quality_score)
                            myReads=Reads(bases,qualities)
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
                            major_count=generics.calcAlleleFreq(major,myReads)
                            minor_count=generics.calcAlleleFreq(minor,myReads)
                            # Take a sample of the bases so that the proportion of data used is as required
                            data_prop = math.ceil(len(myReads.base)*downsampling) # calculate how many bases to include for proportion of sample
                            rand_samp = random.sample(range(0,len(myReads.base)),data_prop)
                            base = ""
                            qualities = ""
                            for r in rand_samp:
                                base+=myReads.base[r]
                                qualities+=myReads.base_quality[r]
                            myReads=Reads(base,qualities)
                            #find sample depth of filtered data    
                            sampleDepth = len(myReads.base)
                            p=[] # list to fill with probabilities for this base



                            if 1 in ploidy:
                                haploid = generics.calcGenoLogLike1_MajorMinor(myReads,mySite,major,minor)
                                haploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,haploid,HWE_Prob_hap)))
                                p.append(haploid_sum)

                            if 2 in ploidy:
                                diploid = generics.calcGenoLogLike2_MajorMinor(myReads,mySite,major,minor)
                                diploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,diploid,HWE_Prob_dip)))
                                p.append(diploid_sum)

                            if 3 in ploidy:
                                triploid = generics.calcGenoLogLike3_MajorMinor(myReads,mySite,major,minor)
                                triploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,triploid,HWE_Prob_tri)))
                                p.append(triploid_sum)

                            if 4 in ploidy:
                                tetraploid = generics.calcGenoLogLike4_MajorMinor(myReads,mySite,major,minor)
                                tetraploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,tetraploid,HWE_Prob_tetra)))
                                p.append(tetraploid_sum)

                            if 5 in ploidy:
                                pentaploid = generics.calcGenoLogLike5_MajorMinor(myReads,mySite,major,minor)
                                pentaploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,pentaploid,HWE_Prob_pent)))
                                p.append(pentaploid_sum)

                            if 6 in ploidy:
                                hexaploid = generics.calcGenoLogLike6_MajorMinor(myReads,mySite,major,minor)
                                hexaploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,hexaploid,HWE_Prob_hex)))
                                p.append(hexaploid_sum)

                            if 7 in ploidy:
                                heptaploid = generics.calcGenoLogLike7_MajorMinor(myReads,mySite,major,minor)
                                heptaploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,heptaploid,HWE_Prob_hept)))
                                p.append(heptaploid_sum)

                            if 8 in ploidy:
                                octaploid = generics.calcGenoLogLike8_MajorMinor(myReads,mySite,major,minor)
                                octaploid_sum = generics.log_or_zero(sum(map(lambda x,y: generics.exp_or_zero(x)*y,octaploid,HWE_Prob_oct)))
                                p.append(octaploid_sum)



                            Overall_Prob_HWE[0]+=p # Add probabilities to overall counter
                            NUMSITES[0]+=sampleDepth # count the number of reads for each sample

                            Overall_Prob_HWE[n]+=p # Add probabilities and depths to sample-wise counter
                            NUMSITES[n]+=sampleDepth 



                            # Write file of genotype likelihoods
                            sep="\t"
                            content=(mySite.chrom,str(mySite.position),str(n),mySite.reference,str(sampleDepth),alleles[major],alleles[minor],str(major_count),str(minor_count),"\t".join(map(str,haploid)),"\t".join(map(str,diploid)),"\t".join(map(str,triploid)),"\t".join(map(str,tetraploid)),"\t".join(map(str,pentaploid)),"\t".join(map(str,hexaploid)),"\t".join(map(str,heptaploid)),"\t".join(map(str,octaploid)))
                            content=sep.join(content)
                            content=content+"\n"
                            with gzip.open(directory+'/'+output,'at+') as f: 
                                f.write(content)
                            content=""


                        # After going through all samples record information for chunks and 
#                        if base_number%win==1:
#                            list_of_window2.append(Overall_Prob_HWE) # make a list of single base probabilities for each sample at every 5th SNP 

#                        if base_number%win==0: #calculate most likely ploidy for each 50 biallelic bases chunks of supercontig
#                            list_of_window.append(Overall_Prob_HWE)

#                            for i in range(NSAMS+1):
#                                max_index=list(Overall_Prob_HWE[i]).index(max(list(Overall_Prob_HWE[i])))
#                                ExpectedPloidy[i].append(max_index+1)      # Make a list of most likely ploidy in each window for each sample                          

#                                Overall_Prob[i]+=Overall_Prob_HWE[i] # Keep a track of the overall probability of each ploidy at each sample


#                            Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float) # reset the probability holder for the next window

                            #end likelihood calc
                        # end for sample
                #end if not filtered for global depth
            #end for line   





    #Ensure enough reads have been passed to for at least one window, if not use all that are available        
#    if base_number<win:
#        list_of_window.append(Overall_Prob_HWE)
#        for i in range(NSAMS+1):
#            Overall_Prob[i]+=Overall_Prob_HWE[i]
#            max_index=list(Overall_Prob_HWE[i]).index(max(list(Overall_Prob_HWE[i])))
#            if max_index==0:
#                ExpectedPloidy[i].append(1)
#            else:
#                ExpectedPloidy[i].append(max_index+1) 


    # Perform bootstap analysis for the overall ploidy level
#    bootstrap_repeats=100 # number of times parameter is estimated
#    bootstrap_sum=list(np.zeros(len(ploidy))) #initiate empty vecotr
#    for i in range(bootstrap_repeats): #cycle through estimations
#        total=list(np.zeros(len(ploidy))) # initiate 
#        for j in range(len(list_of_window)):
#            rand=random.randint(0,len(list_of_window)-1)
#            total=list(map(lambda x,y: x+y ,list(list_of_window[rand][0]),total))
#            est_ploidy=total.index(max(total))
#        bootstrap_sum[est_ploidy]+=1    # bootstrap analysis for ploidy levels
#    print('Bases included: %s' %base_number)
    # Infer most likely ploidy for contig via likelihood
#    Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
#    print("Inferred ploidy is: %s" %Inferred_Ploidy)
#    print("Bootstrap distribtuion of ploidy:")
#    print(bootstrap_sum)




    #Test for aneuploidy     
#    import Delta_anal_poly as DAN
    #import Delta_anal_ANN as DAN # for ANN classifier
#    def aneuploidy_in_sample(x):
#        x=DAN.polynomial_regressor.transform(x)
#        isit=DAN.classifier.predict(x)
#        prob=DAN.classifier.predict_proba(x)
#        return(prob)

#    prob_of_aneu=0


#    H_0=0.0 # place holder for null hypothesis
#    H_1=0.0 # place holder for alternative hypothesis
#    sample_ploidies=list(map(lambda x: list(Overall_Prob[x]).index(max(list(Overall_Prob[x])))+1,[i+1 for i in range(NSAMS) ])) # make a list of sample inferred ploidies
#    Samp_probs=np.zeros((NSAMS,len(ploidy)),float)
#    for i in range(NSAMS*len(list_of_window2)): # cycle through the independently choosen read base
#        rand1=random.randint(1,NSAMS) # choose base from sample at random
#        rand2=random.randint(0,len(list_of_window2)-1) # choose window taken from at random
#        new_part=list_of_window2[rand2][rand1] # record result of this base
#        Samp_probs[rand1-1]+=new_part
#    for i in range(0,NSAMS):
#        H_0+=Samp_probs[i][Inferred_Ploidy-1]
#        H_1+=max(Samp_probs[i])
#    Test=((H_1)-(H_0)) # take the test statistic
#    Ploidy_inferred=Inferred_Ploidy
#    SNPs=NSAMS*len(list_of_window2)
#    average_depth=total_bases/(no_bases*NSAMS*Inferred_Ploidy)
#    Mean_haploid_read_depth=average_depth
#    Normalised_delta=Test/SNPs
#    Samples=NSAMS
#    x=np.array([[Ploidy_inferred,Mean_haploid_read_depth,Normalised_delta,Samples]])
#    prob_of_aneu=aneuploidy_in_sample(x)[0][1]
#    print('Probability of aneuploidy:')
#    print(prob_of_aneu)

    #infer baseline ploidy
#    if base_number>(win*10): 
#        Overall_inferred=[]
#        for i in range(NSAMS):
#            Overall_inferred=Overall_inferred+ExpectedPloidy[i]
#        ploidy_freqs=generics.dist(Overall_inferred)
#        Baseline_Ploidy=ploidy_freqs.index(max(ploidy_freqs))+1
#    else:
#        Baseline_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1

#    samples_included=list(range(1,NSAMS+1))
#    samples_removed=[]
#    while(prob_of_aneu>0.5):
#        if NSAMS>1:
#            # Remove sample least likely to be base ploidy
#            sample_baseline_likelihood=[]
#            for i in samples_included:
#                i=i-1 # sort out indexing
#                sample_baseline_likelihood.append(generics.dist(ExpectedPloidy[i])) # list of likelihood of each sample being baseline sample
#            sample_to_remove=samples_included[sample_baseline_likelihood.index(min(sample_baseline_likelihood))] # Choose sample with lowest prob of baseline sample
#            samples_included.remove(sample_to_remove) # Remove this sample
#            samples_removed.append(sample_to_remove) # Add this sample to list of those removed
#            Overall_Prob[0]=Overall_Prob[0]-Overall_Prob[sample_to_remove] 
            # Revaluate the inferred ploidy without this sample
#            Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
#            NSAMS=NSAMS-1 # Decrease sample size by 1

            #Retest for aneuploidy without removed sample
#            prob_of_aneu=0



#            H_0=0.0 # place holder for null hypothesis
#            H_1=0.0 # place holder for alternative hypothesis
#            sample_ploidies=list(map(lambda x: list(Overall_Prob[x]).index(max(list(Overall_Prob[x])))+1,[i+1 for i in range(NSAMS) ])) # make a list of sample inferred ploidies
#            Samp_probs=np.zeros((NSAMS,len(ploidy)),float)
#            for i in range(NSAMS*len(list_of_window2)): # cycle through the independently choosen read base
#                rand1=random.choice(samples_included)# choose base from sample at random
#                current_index_of_rand1=samples_included.index(rand1)
#                rand2=random.randint(0,len(list_of_window2)-1) # choose window taken from at random
#                new_part=list_of_window2[rand2][rand1] # record result of this base
#                Samp_probs[current_index_of_rand1]+=new_part
#            for i in range(0,NSAMS):
#                H_0+=Samp_probs[i][Inferred_Ploidy-1]
#                H_1+=max(Samp_probs[i])
#            Test=((H_1)-(H_0)) # take the test statistic
#            Ploidy_inferred=Inferred_Ploidy
#            SNPs=NSAMS*len(list_of_window2)
#            average_depth=total_bases/(no_bases*NSAMS*Inferred_Ploidy)
#            Mean_haploid_read_depth=average_depth
#            Normalised_delta=Test/SNPs
#            Samples=NSAMS
#            x=np.array([[Ploidy_inferred,Mean_haploid_read_depth,Normalised_delta,Samples]])
#            prob_of_aneu=aneuploidy_in_sample(x)[0][1]
#            print("Sample being removed {}".format(sample_to_remove))
#            print("New probability of aneuploidy: {}".format(prob_of_aneu))

#        else:
#            samples_removed=list(range(1,Original_sample_number+1))
#            samples_included=[]
#            prob_of_aneu=0

#    print("probability of aneuploidy with samples removed: {}".format(prob_of_aneu))

    # Calculate baseline ploidy by considering windows
    #if base_number>(win*10): 
    #    Overall_inferred=[]
    #    for i in samples_included:
    #        Overall_inferred=Overall_inferred+ExpectedPloidy[i]
    #    ploidy_freqs=generics.dist(Overall_inferred)
    #    Inferred_Ploidy=ploidy_freqs.index(max(ploidy_freqs))+1
    #else:
    #    Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
#    Inferred_Ploidy=Baseline_Ploidy
#    ploidies=""
#    aneuploidy=""
#    content2=""
#    for sample in range(0,Original_sample_number):
#        sample=sample+1
#        if sample in samples_included :
#            sam_ploid=Inferred_Ploidy
#            aneuploid=0
#        else:
#            if base_number>(win*10):
#                sam_ploid = generics.dist(ExpectedPloidy[sample]).index(max(generics.dist(ExpectedPloidy[sample])))+1 # calculate the ploidy of aneuploidy samples
#            else:
#                sam_ploid = list(Overall_Prob[sample]).index(max(list(Overall_Prob[sample])))+1

#            if sam_ploid==Inferred_Ploidy:
#                aneuploid=0
#                samples_included.append(sample)
#                samples_removed.remove(sample)
#            else:
#                aneuploid=1
        #print(list(Overall_Prob[sample]).index(max(list(Overall_Prob[sample])))+1)
        #print(generics.dist(ExpectedPloidy[sample]).index(max(generics.dist(ExpectedPloidy[sample])))+1)
#        ploidies=ploidies+str(sam_ploid)+"\t"
#        aneuploidy=aneuploidy+str(aneuploid)+"\t"
#    ploidies=ploidies.strip("\t")+"\n"
#    aneuploidy=aneuploidy.strip("\t")+"\n"
#    content2=ploidies+aneuploidy

#    print("Samples without aneuploidy and ploidy level")
#    print(samples_included)
#    print(Inferred_Ploidy)

#    print("Samples with aneuploidy")
#    print(samples_removed)
#    with open(directory+'/'+output_2,'wt') as f: 
#        f.write(content2)
#        ploidies=""
#        aneuploidy=""
#        content2=""








