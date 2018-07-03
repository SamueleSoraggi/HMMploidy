import gzip
import generics
import sys
import numpy as np
import argparse
import random
import pickle
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("input",help="file containing the list of basenames for gzipped mpileup files for use in analysis to be used")
args = parser.parse_args()

ploidy=[1,2,3,4,5,6,7,8]
alleles = ['A','C','G','T']


classifier = pickle.load(open('Aneuploidy_Classifier.sav', 'rb'))
polynomial_regressor = pickle.load(open('Polynomial_Regressor.pk1','rb'))



input = args.input #Input file in form of gziped mpileup
list_of_inputs=[]
with open(input,'rb') as f:# opens the mpilup. Use mpileup.read() to display content
    for line in f:
        line=line.decode().strip('\n')# Convert bytes into string
        line=line+".mpileup.gz"
        list_of_inputs.append(line)
Nfiles=len(list_of_inputs)
print('%d files found' %Nfiles)
print(list_of_inputs)



for g in list_of_inputs:
    directory = '/'.join(g.split('/')[:-1])
    if directory =='':
        directory="."
    g = g.split('/')[-1]
    basename='.'.join(g.split('.')[:-2]) # Extract name of file for saving results, -2 to remove .mpileup.gz
    output = basename+".genolikes.gz" # output file to contain genotype likelihoods
    output_2 = basename+".ploids" # output file to contain inferred ploidies and whether aneuploidy or not

    with gzip.open(directory+'/'+g) as f:
        first_line = f.readline()
        Data=first_line.decode().strip('\n')# Convert bytes into string
        l = Data.split('\t')
    NSAMS=int((len(l)-3)/3)
    Original_sample_number=NSAMS

    ExpectedPloidy=[[] for i in range(NSAMS+1)]
    Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float) # (NSAMS+1)xploidies array for probabilities of each ploidy for each sample and overall ploidy probabilities
    Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
    gzip.open(directory+'/'+output_2,'wt')





    Overall_Prob=np.zeros((NSAMS+1,len(ploidy)),float)
    Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float)
    NUMSITES=np.zeros(NSAMS+1,int)
    NUMSITES_HWE=np.zeros(NSAMS+1,int)
    win=50
    list_of_window=[] 
    list_of_window2=[]
    ExpectedPloidy=[[] for i in range(NSAMS+1)]

    base_number=0
    total_bases=0
    no_bases=0
    with gzip.open(directory+'/'+output) as genos:
        for line in genos:
            Data=line.decode().strip('\n')# Convert bytes into string
            l = Data.split('\t')
            sampleDepth=int(l[4])
            n=int(l[2]) # sample number
            total_bases+=int(l[4])
            if n == 1:
                base_number+=1
                no_bases+=1

            p=[]
            if 1 in ploidy:
                haploid = l[9:11]
                haploid_sum = sum(map(lambda x:float(x),haploid))
                p.append(haploid_sum)
            if 2 in ploidy:
                diploid = l[11:14]
                diploid_sum = sum(map(lambda x:float(x),diploid))
                p.append(diploid_sum)
            if 3 in ploidy:
                triploid = l[14:18]
                triploid_sum = sum(map(lambda x:float(x),triploid))
                p.append(triploid_sum)
            if 4 in ploidy:
                tetraploid = l[18:23]
                tetraploid_sum = sum(map(lambda x:float(x),tetraploid))
                p.append(tetraploid_sum)
            if 5 in ploidy:
                pentaploid = l[23:29]
                pentaploid_sum = sum(map(lambda x:float(x),pentaploid))
                p.append(pentaploid_sum)
            if 6 in ploidy:
                hexaploid = l[29:36]
                hexaploid_sum = sum(map(lambda x:float(x),hexaploid))
                p.append(hexaploid_sum)
            if 7 in ploidy:
                heptaploid = l[36:44]
                heptaploid_sum = sum(map(lambda x:float(x),heptaploid))
                p.append(heptaploid_sum)
            if 8 in ploidy:
                octaploid = l[44:53]
                octaploid_sum = sum(map(lambda x:float(x),octaploid))
                p.append(octaploid_sum)
            Overall_Prob_HWE[0]+=p # Add probabilities to overall counter
            NUMSITES[0]+=sampleDepth # count the number of reads for each sa
            Overall_Prob_HWE[n]+=p # Add probabilities and depths to sample-wise counter
            NUMSITES[n]+=sampleDepth 
    # After going through all samples record information for chunks and 
            if base_number%win==1:
                list_of_window2.append(Overall_Prob_HWE) # make a list of single base probabilities for each sample at every 5th SNP 
            if base_number%win==0: #calculate most likely ploidy for each 50 biallelic bases chunks of supercontig
                list_of_window.append(Overall_Prob_HWE)
                for i in range(NSAMS+1):
                    max_index=list(Overall_Prob_HWE[i]).index(max(list(Overall_Prob_HWE[i])))
                    ExpectedPloidy[i].append(max_index+1)      # Make a list of most likely ploidy in each window for each sample                          
                    Overall_Prob[i]+=Overall_Prob_HWE[i] # Keep a track of the overall probability of each ploidy at each sample
                Overall_Prob_HWE=np.zeros((NSAMS+1,len(ploidy)),float) # reset the probability holder for the next window
               #end likelihood calc
            # end for sample
        #end if not filtered for global depth
    #end for line   

    #Ensure enough reads have been passed to for at least one window, if not use all that are available        
    if base_number<win:
        list_of_window.append(Overall_Prob_HWE)
        for i in range(NSAMS+1):
            Overall_Prob[i]+=Overall_Prob_HWE[i]
            max_index=list(Overall_Prob_HWE[i]).index(max(list(Overall_Prob_HWE[i])))
            if max_index==0:
                ExpectedPloidy[i].append(1)
            else:
                ExpectedPloidy[i].append(max_index+1) 


    # Perform bootstap analysis for the overall ploidy level
    bootstrap_repeats=100 # number of times parameter is estimated
    bootstrap_sum=list(np.zeros(len(ploidy))) #initiate empty vecotr
    for i in range(bootstrap_repeats): #cycle through estimations
        total=list(np.zeros(len(ploidy))) # initiate 
        for j in range(len(list_of_window)):
            rand=random.randint(0,len(list_of_window)-1)
            total=list(map(lambda x,y: x+y ,list(list_of_window[rand][0]),total))
            est_ploidy=total.index(max(total))
        bootstrap_sum[est_ploidy]+=1    # bootstrap analysis for ploidy levels
    print('Bases included: %s' %base_number)
    # Infer most likely ploidy for contig via likelihood
    Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
    print("Inferred ploidy is: %s" %Inferred_Ploidy)
    print("Bootstrap distribtuion of ploidy:")
    print(bootstrap_sum)




#TEST FOR ANEUPLOIDY
    #import Delta_anal_ANN as DAN # for ANN classifier
    def aneuploidy_in_sample(x):
        x=polynomial_regressor.transform(x)
        isit=classifier.predict(x)
        prob=classifier.predict_proba(x)
        return(prob)

    prob_of_aneu=0


    H_0=0.0 # place holder for null hypothesis
    H_1=0.0 # place holder for alternative hypothesis
    sample_ploidies=list(map(lambda x: list(Overall_Prob[x]).index(max(list(Overall_Prob[x])))+1,[i+1 for i in range(NSAMS) ])) # make a list of sample inferred ploidies
    Samp_probs=np.zeros((NSAMS,len(ploidy)),float)
    for i in range(NSAMS*len(list_of_window2)): # cycle through the independently choosen read base
        rand1=random.randint(1,NSAMS) # choose base from sample at random
        rand2=random.randint(0,len(list_of_window2)-1) # choose window taken from at random
        new_part=list_of_window2[rand2][rand1] # record result of this base
        Samp_probs[rand1-1]+=new_part
    for i in range(0,NSAMS):
        H_0+=Samp_probs[i][Inferred_Ploidy-1]
        H_1+=max(Samp_probs[i])
    Test=((H_1)-(H_0)) # take the test statistic
    Ploidy_inferred=Inferred_Ploidy
    SNPs=NSAMS*len(list_of_window2)
    average_depth=total_bases/(no_bases*NSAMS*Inferred_Ploidy)
    Mean_haploid_read_depth=average_depth
    Normalised_delta=Test/SNPs
    Samples=NSAMS
    x=np.array([[Ploidy_inferred,Mean_haploid_read_depth,Normalised_delta,Samples]])
    prob_of_aneu=aneuploidy_in_sample(x)[0][1]
    print('Probability of aneuploidy:')
    print(prob_of_aneu)

    #infer baseline ploidy
    if base_number>(win*10): 
        Overall_inferred=[]
        for i in range(NSAMS):
            Overall_inferred=Overall_inferred+ExpectedPloidy[i]
        ploidy_freqs=generics.dist(Overall_inferred)
        Baseline_Ploidy=ploidy_freqs.index(max(ploidy_freqs))+1
    else:
        Baseline_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1

    samples_included=list(range(1,NSAMS+1))
    samples_removed=[]
    while(prob_of_aneu>0.5):
        if NSAMS>1:
            # Remove sample least likely to be base ploidy
            sample_baseline_likelihood=[]
            for i in samples_included:
                i=i-1 # sort out indexing
                sample_baseline_likelihood.append(generics.dist(ExpectedPloidy[i])) # list of likelihood of each sample being baseline sample
            sample_to_remove=samples_included[sample_baseline_likelihood.index(min(sample_baseline_likelihood))] # Choose sample with lowest prob of baseline sample
            samples_included.remove(sample_to_remove) # Remove this sample
            samples_removed.append(sample_to_remove) # Add this sample to list of those removed
            Overall_Prob[0]=Overall_Prob[0]-Overall_Prob[sample_to_remove] 
           # Revaluate the inferred ploidy without this sample
            Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
            NSAMS=NSAMS-1 # Decrease sample size by 1

            #Retest for aneuploidy without removed sample
            prob_of_aneu=0



            H_0=0.0 # place holder for null hypothesis
            H_1=0.0 # place holder for alternative hypothesis
            sample_ploidies=list(map(lambda x: list(Overall_Prob[x]).index(max(list(Overall_Prob[x])))+1,[i+1 for i in range(NSAMS) ])) # make a list of sample inferred ploidies
            Samp_probs=np.zeros((NSAMS,len(ploidy)),float)
            for i in range(NSAMS*len(list_of_window2)): # cycle through the independently choosen read base
                rand1=random.choice(samples_included)# choose base from sample at random
                current_index_of_rand1=samples_included.index(rand1)
                rand2=random.randint(0,len(list_of_window2)-1) # choose window taken from at random
                new_part=list_of_window2[rand2][rand1] # record result of this base
                Samp_probs[current_index_of_rand1]+=new_part
            for i in range(0,NSAMS):
                H_0+=Samp_probs[i][Inferred_Ploidy-1]
                H_1+=max(Samp_probs[i])
            Test=((H_1)-(H_0)) # take the test statistic
            Ploidy_inferred=Inferred_Ploidy
            SNPs=NSAMS*len(list_of_window2)
            average_depth=total_bases/(no_bases*NSAMS*Inferred_Ploidy)
            Mean_haploid_read_depth=average_depth
            Normalised_delta=Test/SNPs
            Samples=NSAMS
            x=np.array([[Ploidy_inferred,Mean_haploid_read_depth,Normalised_delta,Samples]])
            prob_of_aneu=aneuploidy_in_sample(x)[0][1]
            print("Sample being removed {}".format(sample_to_remove))
            print("New probability of aneuploidy: {}".format(prob_of_aneu))
        else:
            samples_removed=list(range(1,Original_sample_number+1))
            samples_included=[]
            prob_of_aneu=0
    print("probability of aneuploidy with samples removed: {}".format(prob_of_aneu))
    #Calculate baseline ploidy by considering windows
    #if base_number>(win*10): 
    #    Overall_inferred=[]
    #    for i in samples_included:
    #        Overall_inferred=Overall_inferred+ExpectedPloidy[i]
    #    ploidy_freqs=generics.dist(Overall_inferred)
    #    Inferred_Ploidy=ploidy_freqs.index(max(ploidy_freqs))+1
    #else:
    #    Inferred_Ploidy=list(Overall_Prob[0]).index(max(list(Overall_Prob[0])))+1
    Inferred_Ploidy=Baseline_Ploidy
    ploidies=""
    aneuploidy=""
    content2=""
    for sample in range(0,Original_sample_number):
        sample=sample+1
        if sample in samples_included :
            sam_ploid=Inferred_Ploidy
            aneuploid=0
        else:
            if base_number>(win*10):
                sam_ploid = generics.dist(ExpectedPloidy[sample]).index(max(generics.dist(ExpectedPloidy[sample])))+1 # calculate the ploidy of aneuploidy samples
            else:
                sam_ploid = list(Overall_Prob[sample]).index(max(list(Overall_Prob[sample])))+1
            if sam_ploid==Inferred_Ploidy:
                aneuploid=0
                samples_included.append(sample)
                samples_removed.remove(sample)
            else:
                aneuploid=1
        print(list(Overall_Prob[sample]).index(max(list(Overall_Prob[sample])))+1)
        print(generics.dist(ExpectedPloidy[sample]).index(max(generics.dist(ExpectedPloidy[sample])))+1)
        ploidies=ploidies+str(sam_ploid)+"\t"
        aneuploidy=aneuploidy+str(aneuploid)+"\t"
    ploidies=ploidies.strip("\t")+"\n"
    aneuploidy=aneuploidy.strip("\t")+"\n"
    content2=ploidies+aneuploidy
    print("Samples without aneuploidy and ploidy level")
    print(samples_included)
    print(Inferred_Ploidy)
    print("Samples with aneuploidy")
    print(samples_removed)
    with open(directory+'/'+output_2,'wt') as f: 
        f.write(content2)
        ploidies=""
        aneuploidy=""
        content2=""
