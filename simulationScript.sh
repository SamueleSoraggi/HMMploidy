#!/bin/bash

POSITIONAL=()

while [[ $# -gt 0 ]]
do
    key="$1"
    minDp=0
    freqs=(0 0 0)
    minInd=1
    FOLDER=.
    SCRIPTFOLDER=.
    BASENAME=out
    
    case $key in
	-s|--simulatorFolder)
	    SCRIPTFOLDER="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-f|--folder)
	    FOLDER="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-o|--out)
	    BASENAME="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-p|--ploidy)
	    ploidy1="$2"
	    IFS=',' read -ra ploidy <<< "$ploidy1"	
	    shift # past argument
	    shift # past value
	    ;;
	-d|--depth)
	    depth1="$2"
	    IFS=',' read -ra depth <<< "$depth1"	
	    shift # past argument
	    shift # past value
	    ;;
	-n|--nSamples)
	    samples1="$2"
	    IFS=',' read -ra samples <<< "$samples1"
	    shift # past argument
	    shift # past value
	    ;;
	-l|--loci)
	    sites="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-g|--minGlobalDepth)
	    minDp="$2"
	    #IFS=',' read -ra minDp <<< "$minDp1"
	    shift # past argument
	    shift # past value
	    ;;
	-m|--minInd)
	    minInd="$2"
	    #IFS=',' read -ra minInd <<< "$minInd1"
	    shift # past argument
	    shift # past value
	    ;;
	-q|--minorFreq)
	    freqs1="$2"
	    IFS=',' read -ra freqs <<< "$freqs1"
	    shift # past argument
	    shift # past value
	    ;;
	--default)
	    DEFAULT=YES
	    shift # past argument
	    ;;
	*)    # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

FILELIST=names.${BASENAME}.filelist
rm -f $FILELIST

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
	    Rscript ${SCRIPTFOLDER}/simulMpileup.R --out test.DP${DP}.NIND${SAM}.txt --copy ${PL}x${SAM} --sites $sites --depth $A --qual 20 --ksfs 1 --ne 10000 --offset $offset | gzip -c -f > $NAME.BUFFER.mpileup.gz

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
	echo $NAME >> $FILELIST

    done
done

echo "End of simulations :)"

#python3 ${SCRIPTFOLDER}/Genotype_Likelihoods.py ${FILEFOLDER} -dp $minDp -m ${freqs[0]} -M2 ${freqs[1]} -M3 ${freqs[2]} ${FILEFOLDER}
