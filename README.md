# Hidden Markov Ploidy
Tool for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the following tool
* [ANGSD](https://github.com/ANGSD/angsd). Once ANGSD is downloaded, substitute the file `abcFreq.cpp` with the one in this github repository, and thereafter compile ANGSD. This allows to model allele frequencies for poliploidy genomes and not only for the diploid case.
* [NGSPOLY](https://github.com/ImperialCollegeLondon/ngsJulia/tree/master/ngsPoly)
* [JULIA vers >= 0.4.7](https://julialang.org/downloads/), with `GZip` and `ArgParse` packages

## Poliploid data simulations
`simulationScript.sh`: simulate poliploidy data: set the ploidy numbers, haploid depths and number of individuals editing the options in the file (ps: I will make this as a proper script with inputs from command line).

#### Input options 
* `ploidy`: vector with the desired sequence of ploidy numbers  
* `depth`: vector with the desired haploid depths
* `sites`: number of sites for each ploidy number
* `samples`: number of individuals to simulate
* `BASENAME`: base name for the output files 
* `FOLDER`: folder where to save the outputs
* In the ANGSD command called in the script, see the explanation for the `minInd` and `minIndDepth` filtering options.

#### Output: 

In the folder `$FOLDER`, there are four groups of data used in ploidy inference, one for each combination of depth and number of individuals (4 combinations in this example) in the following formats and with same basename:
* `.genolikes` where the columns represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 6 (7 columns)
* `.mafs` output from ANGSD calculating allele frequencies at each locus. Each column represents: chromosome, site number, major allele, minor allele, reference allele, estimated frequency, number of individuals with data.
* `.mpileup` mpileup file of the simulated genome
* `.par` initial estimate for ploidy levels parameters. Each two lines represent the alpha and beta parameters of negative binomial distributions modelled on observed depths of one individual.



## Inference of ploidy numbers from data

Inference of ploidy numbers is performed through the R script `hiddenMarkovPloidyShell.R`. This can be used in two different ways:

1) giving in input the base name of a single file and its initial parameters alpha and beta
```Shell
   Rscript hiddenMarkovPloidyShell.R file=#FILE  maxPloidy=$MAXP alpha=$ALPHAS  beta=$BETAS  wind=$WS   chosenInd=$IND  quantileTrim=$Q  minInd=$M
```
2) giving in input a list of base names

```Shell
   Rscript hiddenMarkovPloidyShell.R  fileList=$LIST  maxPloidy=$MAXP  wind=$WS  chosenInd=$IND  quantileTrim=$Q   minInd=$M
```

#### Options

`file`: base name of a group of files with formats `.genolikes`,`.mafs`,`.mpileup`,`.par`.
`fileList`: list of base names for group of files with formats `.genolikes`,`.mafs`,`.mpileup`,`.par`.
`maxPloidy`: max number of ploidy levels (default 6) 
`alpha`: comma separated values of alphas for the depth distributions (default NA)
`beta`: comma separated values of alphas for the depth distributions (default NA)
`chosenInd`: comma separated indices of individuals to analyze (default NA = analyzes all individuals)
`quantileTrim`: comma separated values of 2 quantiles to trim depth values. (default 0,1 = keep all data)
`minInd`: min number of individuals with data for which a locus is usedconsider loci 

### Example: Analyze ploidy numbers from simulations

Simulate a genome called `poliploidyGenome` with sequence of ploidy numbers 2-5-4-2 for two different haploid depths, 3X and 8X, and two different amount of individuals, 5 and 10. Consider `1000` simulated loci for each ploidy level. 
Inside the script, edit the options as it follows
* `ploidy=(2 5 4 2)` 
* `depth=(3 8)`
* `sites=1000`
* `samples=(5 10)`
* `BASENAME=poliploidyGenome`
* `FOLDER=ploidySim`

Run the script:
```Shell
bash simulationScript.sh
```

In `$FOLDER`, for each base name of a simulated dataset, there are four files used for our analysis. For example, for `poliploidyGenome.DP3.NIND10`, one can find
* `poliploidyGenome.DP3.NIND10.genolikes`
* `poliploidyGenome.DP3.NIND10.mafs`
* `poliploidyGenome.DP3.NIND10.mpileup`
* `poliploidyGenome.DP3.NIND10.par`

Create a file `basenames.filelist` with the list of basenames for the simulated genomes. 
Continuing the example from simulated data, the file would for example contain the following basenames:
```
ploidySim/poliploidyGenome.DP3.NIND10
ploidySim/poliploidyGenome.DP3.NIND5
ploidySim/poliploidyGenome.DP8.NIND10
ploidySim/poliploidyGenome.DP8.NIND5
```

Run the analysis of ploidy numbers for the four simulated datasets. Use window size 100, analyze all the individuals, consider the max amount of ploidys being 5, do not trim the data, and use loci where there is data for >=3 individuals

```Shell
Rscript hiddenMarkovPloidyShell.R  fileList=basenames.filelist  maxPloidy=5  wind=100  minInd=3
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
   
