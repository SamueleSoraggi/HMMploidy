# Hidden Markov Ploidy
Tools for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the followings:

* python 3, with packages `gzip, numpy, scipy`
* R, with packages `pracma, data.table, Rcpp`

## Calculate genotype likelihoods from mpileup files

Overview: calculate genotype likelihoods

`Genotype_Likelihoods.py names.filelist`

### Input

* `Input`: name of a text file containing the prefix of each `.mpileup.gz` file (that is, the name of each file excluding the `.mpileup.gz` extension)

### Options

* `-o` or `--outFolder`: Output folder. Default: the folder of each input files
* `-i` or `--Inbreeding`: Inbreeding coefficients for each sample accepted as a comma seperated list e.g `0.3,0.2,0.1` alternatively can take in the format `0.2x3,0.4` which is equivilent to `0.2,0.2,0.2,0.4`. All values must be between 0 and 1. Default value is `0xNSAMS`
* `-d` or `--downsampling`: Fraction of the data to be included included in the calculation of genotype likelihoods and aneuploidy inference. That is for a value `v` in (0,1] for each read there is a `vx100%` chance the base is included in the calculations. this can be used to speed up calculations for high coverage samples. Be careful using this argument for low coverage data. Default: `1`
* `-m` or `--min_non_major_freq`: Set the minimum frequency of non major alleles for bases to be included in the calculations. Default: `0.2`
* `-M2` or `--max_minor2_freq`: Set the maximum frequency of third most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `-M3` or `--max_minor3_freq`: Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `-dp` or `--min_global_depth`: Set the minimum global depth of a base to be included in calculations. All bases with more than this number of reads, after filtering for other conditions mentioned above, across all bases will be included.
* `-dpInd` or `--min_ind_depth`: Set the minimum depth of a base for each sample to included those in the calculations. A locus is not considered if one or more samples have depth lower than the minimum. Default: `0`.

### Output

A `.genolikes.gz` file for each prefix in the input file. The columns of the file represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, major allele counts, minor allele counts, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 8 (9 columns)

### Syntax Example

```Shell
python3 Genotype_Likelihoods.py names.filelist -i 0.1x7,0.15,0.1x2 -d 0.9 -m 0.2 -M2 0.15 -M3 -0.1 -dp 5
```

## Simulation of polyploid data

Overview: simulate poliploidy organisms and output them in an `mpileup.gz` format

`simulationScript.sh -p $FILE -d $DEPTH -l $LOCI -o $OUTNAME`

### Options 

* `-p` or `--ploidyFile`: file containing the desired simulated data (see below for the sintax). 
Each line of the `ploidyFile` contains the ploidy number and the number of individuals for $K$ adjacent segments of a number $J$ of genomes. Example: $K=3$ adjacent segments and $J=20$ genomes. The first segment is diploid for all $20$ genomes, the second with is diploid for 10 genomes and tetraploid for 10 genomes, and the third is triploid for all genomes.
```Shell
2x20
2x10,4x10
3x20
```
The length of each segment is expressed by the `--loci` parameter below.
* `-l` or `--loci`: sequence of $K$ number of loci $l_1,l_2,...,l_K$ used in each simulated segment of the genomes.
* `-d` or `--depth`: vector of haploid depths to simulate. For each depth, a compressed `mpileup` file is created.
* `-o` or `--out`: output file name (without extension). Default: out 

### Output

One `.mpileup.gz` file for each haploid depth, number of individuals $N$, with the name `$OUTNAME.DP$depth.mpileup.gz`.

### Syntax Example

```Shell
$SOFTWAREPATH/simulationScript.sh -p ploidy.file -d 10,20 -l 1000,1000,1000 -o outFile
```

with the file `ploidy.file` containing
```Shell
2x20
2x10,4x10
3x20
```

The output will consists of two files: `outFile.DP10.mpileup.gz` and `outFile.DP20.mpileup.gz`

## Inference of ploidy levels

Overview: infer ploidy levels of each individual in a `.genolikes` file (uncompressed file from the output of the `Genotype_Likelihoods.py` script) using sequencing coverage and genotype likelihoods.

````Shell
Rscript hiddenMarkovPloidyShell.R fileList=$FILELIST wind=100
````
### Input

Text file containing the names of `.genolikes` files to analyze

### Options

* `fileList`: list of base names of files with formats `.genolikes`, uncompressed from the output of `Genotype_Likelihoods.py`. Alternatively, use `file` and write directly the prefix of the desired file (that is, without extension) that has format `.genolikes`
* `wind`: windows size, i.e. nr of loci whose means and genotype likelihoods are summed together
* `maxPloidy`: max number of ploidy levels (default `6`) 
* `nameFile`: file with names of all the individuals in the genotype likelihoods file (default: names are ind1, ind2, ...)
* `chosenInd`: comma separated indices of individuals to analyze (default: analyzes all individuals)
* `quantileTrim`: comma separated values of 2 quantiles to trim depth values. (default `0,1` = keep all data)
* `minInd`: min number of individuals with data for which a locus is usedconsider loci (default 1)
* `eps`: sequencing/mapping error rate (default 0.0005)

Note: 

* if working on simulations, the list of basenames is already given in output by the simulation script.
* the `.genolikes.gz` files given in output by the `Geotype_Likelihoods.py` script must be gunzipped.

### Output

For each base name, there are two outputs:
* a `.pdf` file with inferred ploidy numbers for each individual
* a `.hiddenMarkovPloidy` file where, for each individual, results are arranged on lines as it follows:
   * File name and sample name
   * starting probability vector $\delta$ of inferred ploidy numbers
   * transition matrix $A$ printed on one line (column by column)
   * $\alpha$ vector of parameters for the depth distribution
   * $\beta$ vector of parameters for the depth distributions
   * final loglikelihood of the model
   * inferred ploidy numbers of the model
   * posterior probabilities for the inferred states (column by column)
   * proportion of posterior probability for each inferred ploidy number
   * starting locus of each window of SNPs
   * ending locus of each window of SNPs
   * inferred ploidy level in each of the windows defined in the two lines above

## Application Example: Analyze ploidy numbers from simulations

Simulate 5 genomes ina  file called `poliploidyGenome`, with 4 segments having ploidy numbers 2-5-4-2 for all genomes, and simulating two different haploid depths, 3X and 8X. Consider `1000` simulated loci for each segment. Let `SCRIPTPATH` be the folder containing the scripts of this repository.

Simulate the dataset
```Shell
$SCRIPTPATH/simulationScript.sh -p ploidy.file -d 3,8 -l 1000,1000,1000,1000 -o poliploidyGenome
```

with the file `ploidy.file` containing:
```Shell
2x5
5x5
4x5
2x5
```

There are two simulated `.mpileup.gz` files in output: `poliploidyGenome.DP3.mpileup.gz, poliploidyGenome.DP8.mpileup.gz, and a file containing the list of prefixes for calculating the genotype likelihoods: `poliploidyGenome.filelist`

Calculate genotype likelihoods without applying any filter (see the script options for filtering details):
```Shell
python3 $SCRIPTPATH/Genotype_Likelihoods.py poliploidyGenome.filelist
```

Run the analysis of ploidy numbers for the two simulated datasets. Use window size 100, analyze all the individuals, consider the max ploidy number being 5, do not trim the data, and use loci where there is data for >=2 individuals
```Shell
gunzip *.genolikes.gz #the R script needs gunzipped genolikes files
Rscript $SCRIPTPATH/hiddenMarkovPloidyShell.R  fileList=poliploidyGenome.filelist  maxPloidy=5  wind=100  minInd=2
```

In output you will find the pdf files `poliploidyGenome.DP4.pdf`,`poliploidyGenome.DP8.pdf` with a fancy plot of the inferred ploidy for each individual and the mean depth of each window. The text output is contained in the files `poliploidyGenome.DP4.HMMploidy`,`poliploidyGenome.DP8.HMMploidy`

# PLANNING/EXPECTATION FOR THE NEXT UPDATE:

* SCRIPT FOR COMPARING LOGLIKELIHOOD DIFFERENCES BASED ON A "REFERENCE" HIDDEN MARKOV MODEL
* FILTERING BASED ON THE PROPORTION OF POSTERIOR PROBABILITY FOR EACH PLOIDY
* FILTERING BASED ON OVERLAPPING PARAMETERS OF THE 
