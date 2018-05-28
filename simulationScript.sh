#!/bin/bash

#script to generate simulated ploidy levels

#folders with julia and ngsPoly
JULIA=~/julia/bin/julia 
NGSPOLY=../../BdData/ngsJulia/ngsPoly/ngsPolyLite.jl

ploidy=(2 5 4 2) #sequence of ploidy levels
depth=(3 8) #haploid depth (one simulation for each depth) 
sites=1000 #number of loci for a ploidy level
samples=(5 10) #individuals to simulate (one simulation for each value)
FOLDER=provaScript #folder where to save simulations (must exist)
BASENAME=poliploidyGenome #name of the output file

#generate data for all combinations of parameters (sample and depth)
for SAM in ${samples[@]}
do
    for DP in ${depth[@]}
    do
	NAME=$FOLDER/$BASENAME.DP${DP}.NIND${SAM} #basename for the file
	echo "GENERATING FILE: " $NAME
	offset=1 #loci counter
	rm -f $NAME.mpileup.gz $NAME.txt $NAME.fai
	printf 'Chrom\tStart\tend\n'
	
	for PL in ${ploidy[@]} #sequentially generate ploidy levels
	do

	    A=`Rscript -e "cat($DP*$PL)"` #ploidy level depth
	    #echo $A
	    Rscript ngsPoly/simulMpileup.R --out test.DP${DP}.NIND${SAM}.txt --copy ${PL}x${SAM} --sites $sites --depth $A --qual 20 --ksfs 1 --ne 10000 --offset $offset | gzip -c -f > $NAME.BUFFER.mpileup.gz

	    printf 'copy_%dx%d\t%d\t%d\n' "$PL" "$SAM" "$offset" "$(($offset + $sites - 1))"
	    printf 'copy_%dx%d\t%d\t%d\n' "$PL" "$SAM" "$offset" "$(($offset + $sites - 1))" >> $NAME.fai

	    #keep track of number of loci
	    offset=$(($offset + $sites))

	    #concatenate ploidy levels in a buffer file
	    zcat $NAME.BUFFER.mpileup.gz >> $NAME.BUFFER.txt
	    
	done

	#handle some stuff
	rm -f test.DP${DP}.NIND${SAM}.txt $NAME.BUFFER.mpileup.gz
	cat $NAME.BUFFER.txt | gzip -c > $NAME.mpileup.gz
	rm -f $NAME.BUFFER.txt

	#genotype likelihoods for multiple ploidies
	$JULIA $NGSPOLY --fin $NAME.mpileup.gz --fglikes $NAME.genolikes.gz --nSamples $SAM --minNonMajorCount 2 --minQ 10 --verbose 1
    
	gunzip -f $NAME.mpileup.gz
	gunzip -f $NAME.genolikes.gz

	#calculate allele frequencies
	#------some filtering options:
	#minInd: min number of individuals with data for the locus to be used
	#        the same option is present in the HMM script I made
	#minIndDepth: min depth for each individual to be used
	./angsd/angsd -pileup $NAME.mpileup -GL 1 -out $NAME -doMaf 8 -fai $NAME.fai -nind $SAM -doCounts 1 -p 1 -doMajorMinor 1 -minInd 1 -minIndDepth 1

	gunzip -f $NAME.mafs.gz
	rm -f $NAME.fai
	
	echo "done :)"

    done
done

