# Hidden Markov Ploidy
Tools for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the following tool
* [ANGSD](https://github.com/ANGSD/angsd). Once ANGSD is downloaded, substitute the file `abcFreq.cpp` with the one in this github repository, and thereafter compile ANGSD. This allows to model allele frequencies for poliploidy genomes and not only for the diploid case.
* [NGSPOLY](https://github.com/ImperialCollegeLondon/ngsJulia/tree/master/ngsPoly)
* [JULIA vers >= 0.4.7](https://julialang.org/downloads/), with `GZip` and `ArgParse` packages

## Generating genotype likelihoods

Overview: calculate genotype likelihoods

`Genotype_Likelihoods.py names.filelist`

### Input

* `Input`: name of a text file containing the suffix of each `.mpileup.gz` file

Additionally the following options are available:

* `--Inbreeding [-i]`: Inbreeding coefficients for each sample accepted as a comma seperated list e.g `0.3,0.2,0.1` alternatively can take in the format `0.2x3,0.4` which is equivilent to `0.2,0.2,0.2,0.4`. All values must be between 0 and 1. Default value is `0xNSAMS`
* `--downsampling [-d]`: Fraction of the data to be included included in the calculation of genotype likelihoods and aneuploidy inference. That is for a value `v` in [0,1] for each read there is a `vx100%` chance the base is included in the calculations. this can be used to speed up calculations for high coverage samples. Be careful using this argument for low coverage data. Default: `1`
* `--min_non_major_freq [-m]`: Set the minimum frequency of non major alleles for bases to be included in the calculations. Default: `0.2`
* `--max_minor2_freq [-M2]`: Set the maximum frequency of third most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `--max_minor3_freq [-M3]`: Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `--min_global_depth [-dp]`: Set the minimum global depth of a base to be included in calculations. All bases with more than this number of reads, after filtering for other conditions mentioned above, across all bases will be included.

### Output

A `.genolikes` file for each prefix in the input file. The columns of the file represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, major allele counts, minor allele counts, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 8 (9 columns)

### Syntax Example

```Shell
python Genotype_Likelihoods.py names.filelist -i 0.1x7,0.15,0.1x2 -d 0.9 -m 0.2 -M2 0.15 -M3 -0.1 -dp 5
```

## Simulation of polyploid data mpileup files

Overview: simulate poliploidy `mpileup.gz` files

`simulationScript.sh -s $SCRIPTFOLDER -f $FOLDER -o $OUTNAME -p $PLOIDY -d $DEPTH -n $INDIVIDUALS -l $LOCI`

### Input options 

* `-s` or `--simulatorFolder`: folder containing the script 
* `-f` or `--folder`: folder of the output
* `-o` or `--out`: prefix of the output name to be saved in the folder defined by `-f`
* `-p` or `--ploidy`: vector of $K$ ploidy numbers in the output in the format $p_1,p_2,p_3,...,p_K$
* `-d` or `--depth`: vector of $J$ haploid depths, each used for a simulation $d_1,d_2,d_3,...,d_J$
* `-n` or `--nSamples`: vector of $M$ haploid depths, each used for a simulation $n_1,n_2,n_3,...,n_M$
* `-l` or `--loci`: number of loci used in each chromosome

### Output: 

In the output folder (option `-f`), one `.mpileup.gz` file for each depth $D$, number of individuals $N$, with the name
`$OUTNAME.DP$D.NIND$N.mpileup.gz`.

### Syntax Example

```Shell
$PATH/simulationScript.sh -s $PATH -f folder -o prefix -p 2,4,5,2 -d 10,20 -n 5,10 -l 1000
```

## Application Example: Analyze ploidy numbers from simulations

Simulate a genome called `poliploidyGenome` with sequence of ploidy numbers 2-5-4-2 for two different haploid depths, 3X and 8X, and two different amount of individuals, 5 and 10. Consider `1000` simulated loci for each ploidy level. Let `PATH` be the folder containing the scripts of this repository, and `FOLDER` the folder containing the output.

Simulate the dataset
```Shell
$PATH/simulationScript.sh -s $PATH -f $FOLDER -o poliploidyGenome -p 2,4,5,2 -d 3,8 -n 5,10 -l 1000
```

In `$FOLDER`, for each base name of a simulated dataset, there are four simulated `.mpileup.gz` files: `poliploidyGenome.DP3.NIND10.genolikes, poliploidyGenome.DP8.NIND10.genolikes, poliploidyGenome.DP3.NIND5.genolikes, poliploidyGenome.DP8.NIND5.genolikes`, and a file containing the list of prefixes for calculating the genotype likelihoods: `names.poliploidyGenome.filelist`

Calculate genotype likelihoods without applying any filter (see the script options for filtering details):
```Shell
$PATH/python3 Genotype_Likelihoods.py test.mpileup.gz
```

Run the analysis of ploidy numbers for the four simulated datasets. Use window size 100, analyze all the individuals, consider the max amount of ploidys being 5, do not trim the data, and use loci where there is data for >=3 individuals

```Shell
Rscript hiddenMarkovPloidyShell.R  fileList=names.poliploidyGenome.filelist  maxPloidy=5  wind=100  minInd=3
```

For each basename, there are two outputs:
* a `.pdf` file with inferred ploidy numbers for each individual
* a `.hiddenMarkovPloidy` file where, for each individual, results are arranged on lines as it follows:
   * File name and individual index
   * starting probabilities of inferred ploidies
   * transition matrix printed on one line
   * alpha parameters depth distributions
   * beta parameters depth distributions
   * final loglikelihood of the model
   * inferred ploidy numbers
   * posterior probabilities for the inferred states printed on one line
   * empty line









--------------------------------------

CHANGED BY SAMUELE UNTIL HERE - NEW SCRIPT

--------------------------------------

#### Options

* `fileList`: list of base names for group of files with formats `.genolikes`,`.par`. Alternatively, use `file` and write directly the basename of the desired files in format `.genolikes`,`.par`
* `maxPloidy`: max number of ploidy levels (default 6) 
* `chosenInd`: comma separated indices of individuals to analyze (default NA = analyzes all individuals)
* `quantileTrim`: comma separated values of 2 quantiles to trim depth values. (default 0,1 = keep all data)
* `minInd`: min number of individuals with data for which a locus is usedconsider loci (default 1)
* `eps`: sequencing/mapping error rate (default 0.005)


   
## Inferring alpha and beta values from .genolikes files 
 
Initial estimates for alpha and beta to be fed into `hiddenMarkovPloidyShell.R` can be calculated by using `Fitting_dist_contig.py`.

Currently this is performed by running the script with a list of the input files

### Options

Along with the input file the following options are available:
* `--num_dist [-d]`: Set the number of distributions in the mixture distribtion fitted to the samples. Note not all samples may be able to support higher numbers of distributions and when they can't be fitted the highest supported will be taken instead. Takes integer values 0-5, default = 0 where 0 means number is not fixed and most likely number is choosen for each sample
* `--filter [-f]`: Comma seperated values for the upper lower and upper bound for filtering the depths of data. Default = `0,1` where each value is between 0 and 1 and lower bound must be less than upper bound.

#### Output

Output is in the format of a `filename.par` file for each input file in the file list. Each `filename.par` file will have two tab seperated lines for each individual. The first line will contain the alpha values for individual i and the second line will contain the beta values for individual i.

E.g. an example file of M individuals would look as follows:
```
individual1.alpha1   individual1.alpha2   individual1.alpha3 ... individual1.alphaN
individual2.beta1   individual2.beta2   individual2.beta3 ... individual2.betaN
.
.
.
individualM.alpha1   individualM.alpha2  individualM.alpha3 ... individualM.alphaN
individualM.beta1   individualM.beta2  individualM.beta3 ... individualM.betaN
```
### Example

Run the script as follows:

```Shell
python Fitting_dist_contig.py basenames.filelist --num_dist 4 --filter 0.1,0.9
```
Where basenames.filelist is the file containg the list of input files.


## Generating genolikelihoods and inferring probability of aneuploidy

Genolikelihoods can be calculated and recorded as a `.genolikes` file by using `Genotype_Likelihoods.py`. This will also return useful information about the probability of aneuploidy with a set of samples to the screen and record an inferred ploidy for each individual in the form of a `.ploids` file

### Inputs

`Genotype_Likelihoods.py` is run with the required arguments of:
* `Input`: A gzipped mpileup file containing the data for the samples to be analysed

Additionally the following options are available:
* `--Inbreeding [-i]`: Inbreeding coefficients for each sample accepted as a comma seperated list e.g `0.3,0.2,0.1` alternatively can take in the format `0.2x3,0.4` which is equivilent to `0.2,0.2,0.2,0.4`. All values must be between 0 and 1. Default value is `0xNSAMS`
* `--downsampling [-d]`: Fraction of the data to be included included in the calculation of genotype likelihoods and aneuploidy inference. That is for a value `v` in [0,1] for each read there is a `vx100%` chance the base is included in the calculations. this can be used to speed up calculations for high coverage samples. Be careful using this argument for low coverage data. Default: `1`
* `--min_non_major_freq [-m]`: Set the minimum frequency of non major alleles for bases to be included in the calculations. Default: `0.2`
* `--max_minor2_freq [-M2]`: Set the maximum frequency of third most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `--max_minor3_freq [-M3]`: Set the maximum frequency of fourth most prolific alleles for bases to be included in the calculations. Used to determine strengh of confidence on bases being biallelic. Default: `0.1`
* `--min_global_depth [-dp]`: Set the minimum global depth of a base to be included in calculations. All bases with more than this number of reads, after filtering for other conditions mentioned above, across all bases will be included.

### Outputs

On screen outputs will give you the number of bases included in the calculations, the most likely ploidy based on genoptype likelihoods alone (Inferred ploidy), The bootstrap distribution of this inferred ploidy by taking random bases, The probability there ius aneuploidy in the entire dataset, the list of samples likely to not have aneuploidy and their inferred baseline ploidy, the list of samples likely to have aneuploidy and the probability of aneuploidy within the dataset after the samples inferred tro have aneuploidy are removed.


In addition to the on-screen outputs the two output files are as follows:

#### `.genolikes`
where the columns represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, major allele counts, minor allele counts, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 8 (9 columns)

#### `.ploids`
Where the first row represents the inferred ploidy for each individual in the data set and the second row is a list of 1 and 0 where 1 represents the sample was inferred to have aneuploidy and 0 that it was not inferred to have aneuploidy.


### Example

With an mpileup file `test.mpileup.gz`, the program can be run as follows:


```Shell
python Genotype_Likelihoods.py test.mpileup.gz -i 0.1x7,0.15,0.1x2 -d 0.9 -m 0.2 -M2 0.15 -M3 -0.1 -dp 5
```






