#!/bin/bash

POSITIONAL=()

while [[ $# -gt 0 ]]
do
    key="$1"
    minDp=0
    freqs=(0 0 0)
    BASENAME=./out
    pvar=1
    
    case $key in
	-o|--out)
	    BASENAME="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-p|--ploidyFile)
	    ploidy="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-d|--depth)
	    depth1="$2"
	    IFS=',' read -ra depth <<< "$depth1"	
	    shift # past argument
	    shift # past value
	    ;;
	#-n|--nSamples)
	#    samples1="$2"
	#    IFS=',' read -ra samples <<< "$samples1"
	#    shift # past argument
	#    shift # past value
	#    ;;
	-l|--loci)
	    sites1="$2"
	    IFS=',' read -ra sites <<< "$sites1"
	    shift # past argument
	    shift # past value
	    ;;
	-v|--verbose)
	    VERBOSE=TRUE
	    shift # past argument
	    shift # past value
	    ;;
	-r|--pvar)
	    pvarInput="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-g|--rndGenerator)
	    seedInput="$2"
	    shift # past argument
	    shift # past value
	    ;;
	#-m|--minInd)
	#    minInd="$2"
	#    #IFS=',' read -ra minInd <<< "$minInd1"
	#    shift # past argument
	#    shift # past value
	#    ;;
	#-q|--minorFreq)
	#    freqs1="$2"
	#    IFS=',' read -ra freqs <<< "$freqs1"
	#    shift # past argument
	#    shift # past value
	#    ;;
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

#reformat to assign absolute path, so that the list of files
#can be used from anywhere 
BBB=`realpath $BASENAME`
BASENAME=`dirname $BBB`/`basename $BASENAME`

FILELIST=${BASENAME}.filelist
rm -f $FILELIST

#FIND SOURCE FOLDER OF THE SCRIPTS
SCRIPTFOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#NLINES=`cat $ploidy | wc -l` does not count last line if there is not a \n
NLINES=`grep -c "" $ploidy` 
echo $NLINES

#generate data for all combinations of parameters (sample and depth)
for DP in ${depth[@]}
do

    SITESCOUNTER=0
    NAME=$BASENAME.DP${DP} #basename for the file
    echo "GENERATING FILE: " "\"" $NAME "\""
    offset=1 #loci counter
    rm -f $NAME.mpileup.gz $NAME.txt
    printf 'Chrom\tStart\tend\n'
    
    for CFGLINE in `seq 1 $NLINES`
    do
	#A=`sed -n -e ${CFGLINE}p $ploidy`
	#echo $A
	#A=`Rscript -e "cat($DP*$PL)"` #ploidy level depth
	#echo $A
	Rscript ${SCRIPTFOLDER}/simulMpileup.R --out test.DP${DP}.txt --copy `sed -n -e ${CFGLINE}p $ploidy` --sites ${sites[$SITESCOUNTER]} --depth $DP --qual 20 --ksfs 1 --ne 10000 --offset $offset --seed $seedInput --pvar $pvarInput | gzip -c -f > $NAME.BUFFER.mpileup.gz
	
	printf '%s\t%d\t%d\n' `sed -n -e ${CFGLINE}p $ploidy` "$offset" "$(($offset + ${sites[$SITESCOUNTER]} - 1))"
	    
	#printf 'copy_%dx%d\t%d\t%d\n' "$PL" "$SAM" "$offset" "$(($offset + ${sites[$SITESCOUNTER]} - 1))" >> $NAME.fai

	#keep track of number of loci
	offset=$(($offset+${sites[$SITESCOUNTER]}))
	
	#concatenate ploidy levels in a buffer file
	zcat $NAME.BUFFER.mpileup.gz >> $NAME.BUFFER.txt

	SITESCOUNTER=$(($SITESCOUNTER + 1))
    done

    #handle some stuff
    rm -f test.DP${DP}.txt $NAME.BUFFER.mpileup.gz
    cat $NAME.BUFFER.txt | gzip -c > $NAME.mpileup.gz
    rm -f $NAME.BUFFER.txt
    echo $NAME >> $FILELIST
    
done

echo "output file list " "\"" $FILELIST "\""
