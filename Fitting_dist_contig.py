#! usr/bin/python3


####Required packages for installation
import numpy as np # organising the depths
import scipy as sc # Opperations such as percentile
import scipy.stats as scs # Statistical functions such as poisson distribution and percentiles
import pandas as pd # Creating and manipulating depthsframes of results
import lmfit # Fitting the models to the depths
import matplotlib.pyplot as plt # plotting the models fitted
import re # Read the results and extract the values
import sys # allow command line arguments of input and output files
import os # To subprocess the extraction of depth files from mpileup files
import argparse


########### Negative Binomial mixture models ############
### see below for parameter definitions

def neg_binomial_1(Neg_Binom_Params):
    u=Neg_Binom_Params.valuesdict()
    model_value=scs.nbinom(u['n1'],u['p1']).pmf(vals)
    residuals=model_value-actual
    return(residuals) 
def neg_binomial_2(Neg_Binom_Params):
    u=Neg_Binom_Params.valuesdict()
    model_value=(u['coef1'])*scs.nbinom(u['n1'],u['p1']).pmf(vals)+(1-u['coef1'])*scs.nbinom(u['n2'],u['p2']).pmf(vals)
    residuals=model_value-actual
    return(residuals) 
def neg_binomial_3(Neg_Binom_Params):
    u=Neg_Binom_Params.valuesdict()
    model_value=(u['coef1'])*scs.nbinom(u['n1'],u['p1']).pmf(vals)+(u['coef2'])*scs.nbinom(u['n2'],u['p2']).pmf(vals)+(1-u['coef1']-u['coef2'])*scs.nbinom(u['n3'],u['p3']).pmf(vals)
    residuals=model_value-actual
    return(residuals) 
def neg_binomial_4(Neg_Binom_Params):
    u=Neg_Binom_Params.valuesdict()
    model_value=(u['coef1'])*scs.nbinom(u['n1'],u['p1']).pmf(vals)+(u['coef2'])*scs.nbinom(u['n2'],u['p2']).pmf(vals)+(u['coef3'])*scs.nbinom(u['n3'],u['p3']).pmf(vals)+(1-u['coef1']-u['coef2']-u['coef3'])*scs.nbinom(u['n4'],u['p4']).pmf(vals)
    residuals=model_value-actual
    return(residuals) 
def neg_binomial_5(Neg_Binom_Params):
    u=Neg_Binom_Params.valuesdict()
    model_value=(u['coef1'])*scs.nbinom(u['n1'],u['p1']).pmf(vals)+(u['coef2'])*scs.nbinom(u['n2'],u['p2']).pmf(vals)+(u['coef3'])*scs.nbinom(u['n3'],u['p3']).pmf(vals)+(u['coef4'])*scs.nbinom(u['n4'],u['p4']).pmf(vals)+(1-u['coef1']-u['coef2']-u['coef3']-u['coef4'])*scs.nbinom(u['n5'],u['p5']).pmf(vals)
    residuals=model_value-actual
    return(residuals) 


###### Parse arguments

parser = argparse.ArgumentParser()
parser.add_argument("input",help="File containing a list of basenames for the .genolike files to be used")
parser.add_argument("-d","--num_dist",help="Fix the number of distributions in the mixture distribution",default=0)
parser.add_argument("-f","--filter",help="comma seperated decimals for upper and lower levels of filtering. Default = 0,1",default='0,1')
args=parser.parse_args()


input=args.input
fix_sample=args.num_dist
l=int(float(args.filter.split(',')[0])*100)
u=int(float(args.filter.split(',')[1])*100)


list_of_inputs=[]
with open(input,'rb') as f:# opens the mpilup. Use mpileup.read() to display content
    for line in f:
        line=line.decode().strip('\n')# Convert bytes into string
        line=line+".genolikes"
        list_of_inputs.append(line)
Nfiles=len(list_of_inputs)
print('%d files found' %Nfiles)
print(list_of_inputs)


###### Creation of depth files
for g in list_of_inputs:
    overall_content=""
    contig='.'.join(g.split('.')[:-1]) # Extract name of file for saving results
    print("Using file: {}".format(g))
    # Extract NSAMS
    with open(g,'rb') as genos:
        old_sam=0
        test=0
        while test==0:
            for line in genos:
                line=line.decode().strip('\n') # Convert to text
                line=line.split('\t') # Sperate into elements
                SAM=int(line[2]) #Extract sample value
                if SAM>old_sam:
                    old_sam=SAM
                else:
                    NSAMS=old_sam
                    test=1
    print("NSAMS = {}".format(NSAMS))
                


    with open(g,'rb') as genos: # For each file
        Depth=[[] for n in range(NSAMS+1)] # Create empty array for depths
        for line in genos: # Cycle through lines
            line=line.decode().strip('\n') # Convert to text
            line=line.split('\t') # Sperate into elements
            SAM=int(line[2]) #Extract sample value
            if SAM>len(Depth):
                print("More samples found than inputted value! Please check NSAMS parameter")
                sys.exit(1)
            depth=int(line[4]) #Extract depth value
            Depth[SAM-1].append(depth) #stores depths in memory, may be better to not
        Depth=[d for d in Depth if d !=[]]
    # Store depths in depth files, seperate for each sample
    genos.close() # close genotype likelihood file to save memory
    Sample=0
    for depths in Depth:
        Sample+=1
        # filter depths
        upper=np.percentile(depths,u) # set upper filter level
        lower=np.percentile(depths,l) # set lower filter level
        depths = list(filter(lambda x : x < upper and x > lower , depths)) # filter the depths
        
        #set up required value holders

        Neg_Binom_Actual=[]
        Neg_Binom_min=[]
        MUS1_NB=[]
        MUS2_NB=[]
        MUS3_NB=[]
        MUS4_NB=[]
        MUS5_NB=[]
        COEF1_NB=[]
        COEF2_NB=[]
        COEF3_NB=[]
        COEF4_NB=[]
        COEF5_NB=[]
        Neg_Binom=[]
        Neg_Binom_Vals=[]
        

        # create parameters for the models
        Neg_Binom_Params_1=lmfit.Parameters()
        Neg_Binom_Params_1.add('n1',value=sum(depths)/(40),min=1)
        Neg_Binom_Params_1.add('p1',value=1-np.percentile(depths,50)/(sum(depths)/(40)),min=0.001,max=.999)
        
        Neg_Binom_Params_2=lmfit.Parameters()
        Neg_Binom_Params_2.add('n1',value=sum(depths)/(2*40),min=1)
        Neg_Binom_Params_2.add('p1',value=1-np.percentile(depths,33)/(sum(depths)/(2*40)),min=0.001,max=.999)
        Neg_Binom_Params_2.add('n2',value=sum(depths)/(2*40),min=1)
        Neg_Binom_Params_2.add('p2',value=1-np.percentile(depths,66)/(sum(depths)/(2*40)),min=0.001,max=.999)        
        Neg_Binom_Params_2.add('coef1',value=0.5,max=1,min=0.001)
         
        Neg_Binom_Params_3=lmfit.Parameters()
        Neg_Binom_Params_3.add('n1',value=sum(depths)/(3*40),min=1)
        Neg_Binom_Params_3.add('p1',value=1-np.percentile(depths,25)/(sum(depths)/(3*40)),min=0.001,max=.999)
        Neg_Binom_Params_3.add('n2',value=sum(depths)/(3*40),min=1)
        Neg_Binom_Params_3.add('p2',value=1-np.percentile(depths,50)/(sum(depths)/(3*40)),min=0.001,max=.999)
        Neg_Binom_Params_3.add('n3',value=sum(depths)/(3*40),min=1)
        Neg_Binom_Params_3.add('p3',value=1-np.percentile(depths,75)/(sum(depths)/(3*40)),min=0.001,max=.999)        
        Neg_Binom_Params_3.add('coef1',value=0.33,max=1,min=0.001)
        Neg_Binom_Params_3.add('coef2',value=0.33,max=1,min=0.001)
         
        Neg_Binom_Params_4=lmfit.Parameters()
        Neg_Binom_Params_4.add('n1',value=sum(depths)/(4*40),min=1)
        Neg_Binom_Params_4.add('p1',value=1-np.percentile(depths,20)/(sum(depths)/(4*40)),min=0.001,max=.999)
        Neg_Binom_Params_4.add('n2',value=sum(depths)/(4*40),min=1)
        Neg_Binom_Params_4.add('p2',value=1-np.percentile(depths,40)/(sum(depths)/(4*40)),min=0.001,max=.999)
        Neg_Binom_Params_4.add('n3',value=sum(depths)/(4*40),min=1)
        Neg_Binom_Params_4.add('p3',value=1-np.percentile(depths,60)/(sum(depths)/(4*40)),min=0.001,max=.999)
        Neg_Binom_Params_4.add('n4',value=sum(depths)/(4*40),min=1)
        Neg_Binom_Params_4.add('p4',value=1-np.percentile(depths,80)/(sum(depths)/(4*40)),min=0.001,max=.999)        
        Neg_Binom_Params_4.add('coef1',value=0.25,max=1,min=0.001)
        Neg_Binom_Params_4.add('coef2',value=0.25,max=1,min=0.001)
        Neg_Binom_Params_4.add('coef3',value=0.25,max=1,min=0.001)
         
        Neg_Binom_Params_5=lmfit.Parameters()
        Neg_Binom_Params_5.add('n1',value=sum(depths)/(5*40),min=1)
        Neg_Binom_Params_5.add('p1',value=1-np.percentile(depths,17)/(sum(depths)/(5*40)),min=0.001,max=.999)
        Neg_Binom_Params_5.add('n2',value=sum(depths)/(5*40),min=1)
        Neg_Binom_Params_5.add('p2',value=1-np.percentile(depths,33)/(sum(depths)/(5*40)),min=0.001,max=.999)
        Neg_Binom_Params_5.add('n3',value=sum(depths)/(5*40),min=1)
        Neg_Binom_Params_5.add('p3',value=1-np.percentile(depths,50)/(sum(depths)/(5*40)),min=0.001,max=.999)
        Neg_Binom_Params_5.add('n4',value=sum(depths)/(5*40),min=1)
        Neg_Binom_Params_5.add('p4',value=1-np.percentile(depths,67)/(sum(depths)/(5*40)),min=0.001,max=.999)
        Neg_Binom_Params_5.add('n5',value=sum(depths)/(5*40),min=1)
        Neg_Binom_Params_5.add('p5',value=1-np.percentile(depths,84)/(sum(depths)/(5*40)),min=0.001,max=.999)          
        Neg_Binom_Params_5.add('coef1',value=0.2,max=1,min=0.001)
        Neg_Binom_Params_5.add('coef2',value=0.2,max=1,min=0.001)
        Neg_Binom_Params_5.add('coef3',value=0.2,max=1,min=0.001)
        Neg_Binom_Params_5.add('coef4',value=0.2,max=1,min=0.001)


        # prepare data
        x = list(np.arange(int(min(depths)),int(max(depths)),1)) # Make list from min to max occurence step of 1
        vals=x
        y=[]    # set an empty list
        for j in x:
            y.append(depths.count(j)/len(depths))  # Add the proportions of the total count that match each value    
        actual=y

        # fit the models

        try:

            Neg_Binom_1_fit =lmfit.minimize(neg_binomial_1,Neg_Binom_Params_1)
            Neg_Binom.append(lmfit.fit_report(Neg_Binom_1_fit))
            Neg_Binom_Vals.append(Neg_Binom_1_fit.residual+actual)
            print("fit NB1")
            counter2=1
            Neg_Binom_2_fit =lmfit.minimize(neg_binomial_2,Neg_Binom_Params_2)
            Neg_Binom.append(lmfit.fit_report(Neg_Binom_2_fit))
            Neg_Binom_Vals.append(Neg_Binom_2_fit.residual+actual)
            print("fit NB2")
            counter2=2
            Neg_Binom_3_fit =lmfit.minimize(neg_binomial_3,Neg_Binom_Params_3)
            Neg_Binom.append(lmfit.fit_report(Neg_Binom_3_fit))
            Neg_Binom_Vals.append(Neg_Binom_3_fit.residual+actual)
            print("fit NB3")
            counter2=3
            Neg_Binom_4_fit =lmfit.minimize(neg_binomial_4,Neg_Binom_Params_4)
            Neg_Binom.append(lmfit.fit_report(Neg_Binom_4_fit))
            Neg_Binom_Vals.append(Neg_Binom_4_fit.residual+actual)
            print("fit NB4")
            counter2=4
            Neg_Binom_5_fit =lmfit.minimize(neg_binomial_5,Neg_Binom_Params_5)
            Neg_Binom.append(lmfit.fit_report(Neg_Binom_5_fit))   
            Neg_Binom_Vals.append(Neg_Binom_5_fit.residual+actual)     
            print("fit NB5")
            counter2=5
        except (OSError,TypeError,ValueError):
            pass

        # Determine best fitting model

        AIC_Neg_Binom=[]

        for k in range(counter2):
            A_NB = re.search('(Akaike\sinfo\scrit.*=\s)(.*)',Neg_Binom[k])
            AIC_Neg_Binom.append(float(A_NB.group(2)))

        ### Retrieve parameters for best model from the fit report.
        ### Ensuring if an error occured when fitting  model that the values for that ad more complex models not fitted are 0.
        if fix_sample==0:
            Best_NB = AIC_Neg_Binom.index(min(AIC_Neg_Binom))
        else:
            Best_NB = min([int(fix_sample),counter2])-1

        #### N_i is the number of failures
        #### P_i is the probability of failure
        
        if Best_NB == 4:
               
            N_1 = re.search('(\n\s*n1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_2 = re.search('(\n\s*n2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_3 = re.search('(\n\s*n3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_4 = re.search('(\n\s*n4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_5 = re.search('(\n\s*n5:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_1 = re.search('(\n\s*p1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_2 = re.search('(\n\s*p2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_3 = re.search('(\n\s*p3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_4 = re.search('(\n\s*p4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_5 = re.search('(\n\s*p5:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])   
            alpha_1=(float(P_1.group(2)))/(1-float(P_1.group(2)))
            alpha_2=(float(P_2.group(2)))/(1-float(P_2.group(2)))
            alpha_3=(float(P_3.group(2)))/(1-float(P_3.group(2)))
            alpha_4=(float(P_4.group(2)))/(1-float(P_4.group(2)))
            alpha_5=(float(P_5.group(2)))/(1-float(P_5.group(2)))
            beta_1=float(N_1.group(2))
            beta_2=float(N_2.group(2))
            beta_3=float(N_3.group(2))
            beta_4=float(N_4.group(2))
            beta_5=float(N_5.group(2))
            MUS1_NB.append(scs.nbinom.mean(float(N_1.group(2)),float(P_1.group(2))))
            MUS2_NB.append(scs.nbinom.mean(float(N_2.group(2)),float(P_2.group(2))))
            MUS3_NB.append(scs.nbinom.mean(float(N_3.group(2)),float(P_3.group(2))))
            MUS4_NB.append(scs.nbinom.mean(float(N_4.group(2)),float(P_4.group(2))))
            MUS5_NB.append(scs.nbinom.mean(float(N_5.group(2)),float(P_5.group(2))))
            Co_1 = re.search('(\n\s*coef1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_2 = re.search('(\n\s*coef2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_3 = re.search('(\n\s*coef3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_4 = re.search('(\n\s*coef4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_5 = re.search('(\n\s*coef5:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])            
            COEF1_NB.append(float(Co_1.group(2)))
            COEF2_NB.append(float(Co_2.group(2)))
            COEF3_NB.append(float(Co_3.group(2)))
            COEF4_NB.append(float(Co_4.group(2)))
            COEF5_NB.append(1-float(Co_1.group(2))-float(Co_2.group(2))-float(Co_3.group(2))-float(Co_4.group(2)))
        elif Best_NB == 3:
            
            N_1 = re.search('(\n\s*n1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_2 = re.search('(\n\s*n2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_3 = re.search('(\n\s*n3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_4 = re.search('(\n\s*n4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_1 = re.search('(\n\s*p1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_2 = re.search('(\n\s*p2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_3 = re.search('(\n\s*p3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_4 = re.search('(\n\s*p4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            alpha_1=(float(P_1.group(2)))/(1-float(P_1.group(2)))
            alpha_2=(float(P_2.group(2)))/(1-float(P_2.group(2)))
            alpha_3=(float(P_3.group(2)))/(1-float(P_3.group(2)))
            alpha_4=(float(P_4.group(2)))/(1-float(P_4.group(2)))
            alpha_5=0
            beta_1=float(N_1.group(2))
            beta_2=float(N_2.group(2))
            beta_3=float(N_3.group(2))
            beta_4=float(N_4.group(2))
            beta_5=0
            MUS1_NB.append(scs.nbinom.mean(float(N_1.group(2)),float(P_1.group(2))))
            MUS2_NB.append(scs.nbinom.mean(float(N_2.group(2)),float(P_2.group(2))))
            MUS3_NB.append(scs.nbinom.mean(float(N_3.group(2)),float(P_3.group(2))))
            MUS4_NB.append(scs.nbinom.mean(float(N_4.group(2)),float(P_4.group(2))))
            MUS5_NB.append(0)
            Co_1 = re.search('(\n\s*coef1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_2 = re.search('(\n\s*coef2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_3 = re.search('(\n\s*coef3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_4 = re.search('(\n\s*coef4:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])            
            COEF1_NB.append(float(Co_1.group(2)))
            COEF2_NB.append(float(Co_2.group(2)))
            COEF3_NB.append(float(Co_3.group(2)))
            COEF4_NB.append(1-float(Co_1.group(2))-float(Co_2.group(2))-float(Co_3.group(2)))
            COEF5_NB.append(0)
        elif Best_NB == 2:
            
            N_1 = re.search('(\n\s*n1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_2 = re.search('(\n\s*n2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_3 = re.search('(\n\s*n3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_1 = re.search('(\n\s*p1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_2 = re.search('(\n\s*p2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_3 = re.search('(\n\s*p3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            alpha_1=(float(P_1.group(2)))/(1-float(P_1.group(2)))
            alpha_2=(float(P_2.group(2)))/(1-float(P_2.group(2)))
            alpha_3=(float(P_3.group(2)))/(1-float(P_3.group(2)))
            alpha_4=0
            alpha_5=0
            beta_1=float(N_1.group(2))
            beta_2=float(N_2.group(2))
            beta_3=float(N_3.group(2))
            beta_4=0
            beta_5=0
            MUS1_NB.append(scs.nbinom.mean(float(N_1.group(2)),float(P_1.group(2))))
            MUS2_NB.append(scs.nbinom.mean(float(N_2.group(2)),float(P_2.group(2))))
            MUS3_NB.append(scs.nbinom.mean(float(N_3.group(2)),float(P_3.group(2))))
            MUS4_NB.append(0)
            MUS5_NB.append(0)
            Co_1 = re.search('(\n\s*coef1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_2 = re.search('(\n\s*coef2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_3 = re.search('(\n\s*coef3:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])           
            COEF1_NB.append(float(Co_1.group(2)))
            COEF2_NB.append(float(Co_2.group(2)))
            COEF3_NB.append(1-float(Co_1.group(2))-float(Co_2.group(2)))
            COEF4_NB.append(0)
            COEF5_NB.append(0)   
        elif Best_NB == 1:
            
            N_1 = re.search('(\n\s*n1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            N_2 = re.search('(\n\s*n2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])  
            P_1 = re.search('(\n\s*p1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            P_2 = re.search('(\n\s*p2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            alpha_1=(float(P_1.group(2)))/(1-float(P_1.group(2)))
            alpha_2=(float(P_2.group(2)))/(1-float(P_2.group(2)))
            alpha_3=0
            alpha_4=0
            alpha_5=0
            beta_1=float(N_1.group(2))
            beta_2=float(N_2.group(2))
            beta_3=0
            beta_4=0
            beta_5=0
            MUS1_NB.append(scs.nbinom.mean(float(N_1.group(2)),float(P_1.group(2))))
            MUS2_NB.append(scs.nbinom.mean(float(N_2.group(2)),float(P_2.group(2))))
            MUS3_NB.append(0)
            MUS4_NB.append(0)
            MUS5_NB.append(0)
            Co_1 = re.search('(\n\s*coef1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            Co_2 = re.search('(\n\s*coef2:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])            
            COEF1_NB.append(float(Co_1.group(2)))
            COEF2_NB.append(1-float(Co_1.group(2)))
            COEF3_NB.append(0)
            COEF4_NB.append(0)
            COEF5_NB.append(0)
        else:
            
            N_1 = re.search('(\n\s*n1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])  
            P_1 = re.search('(\n\s*p1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB]) 
            alpha_1=(float(P_1.group(2)))/(1-float(P_1.group(2)))
            alpha_2=0
            alpha_3=0
            alpha_4=0
            alpha_5=0
            beta_1=float(N_1.group(2))
            beta_2=0
            beta_3=0
            beta_4=0
            beta_5=0
            MUS1_NB.append(scs.nbinom.mean(float(N_1.group(2)),float(P_1.group(2))))
            MUS2_NB.append(0)
            MUS3_NB.append(0)
            MUS4_NB.append(0)
            MUS5_NB.append(0)
            Co_1 = re.search('(\n\s*coef1:\s*)(\d*\.?\d*e?\+?\d*)(\s*\+/-\s*)(\d*.\d*e?\+?\d*)',Neg_Binom[Best_NB])            
            COEF1_NB.append(1)
            COEF2_NB.append(0)
            COEF3_NB.append(0)
            COEF4_NB.append(0)
            COEF5_NB.append(0)          



        #### take basic stats on which model is best and compare
        #Neg_Binom_min.append(min(AIC_Neg_Binom)) # minimum AIC value for the Negative Binomial models
        #Best_Fit_Neg_Binom=(AIC_Neg_Binom.index(min(AIC_Neg_Binom))+1) # Index of best fit model + 1 so that have the number of curves in best fitted model
        #Neg_Binom_Actual.append(Best_Fit_Neg_Binom)
        #### extract lists of the fitted values of best distributions for plotting 
        #try:
        #    Neg_Binom_1_Vals = list(Neg_Binom_Vals[0])
        #    Neg_Binom_2_Vals = list(Neg_Binom_Vals[1])
        #    Neg_Binom_3_Vals = list(Neg_Binom_Vals[2])
        #    Neg_Binom_4_Vals = list(Neg_Binom_Vals[3])
        #    Neg_Binom_5_Vals = list(Neg_Binom_Vals[4])
        #    Neg_Binom_Best_Model_Vals = list(Neg_Binom_Vals[AIC_Neg_Binom.index(min(AIC_Neg_Binom))])   
        #except (OSError,TypeError,ValueError,IndexError):
        #    pass


        #### refresh content to be written
        alphas=""
        betas=""
        content=""


        #### alphas and betas

        alphas="\t".join((str(alpha_1),str(alpha_2),str(alpha_3),str(alpha_4),str(alpha_5)))
        betas="\t".join((str(beta_1),str(beta_2),str(beta_3),str(beta_4),str(beta_5)))    

        content="\n".join((alphas,betas))
        content=content+"\n"



        

        overall_content+=content




    #choose correct file to print results into 
    filename="%s.par" %contig

    with open(filename,'wt') as h:
        h.write(overall_content)
        h.close()
