# Hidden Markov Ploidy
Tool for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the following tool
* [ANGSD](https://github.com/ANGSD/angsd)
* [NGSPOLY](https://github.com/ImperialCollegeLondon/ngsJulia/tree/master/ngsPoly)
* [JULIA vers >= 0.4.7](https://julialang.org/downloads/)

Once ANGSD is downloaded, substitute the file abcFreq.cpp with the one in this github repository, and thereafter compile ANGSD. This allows to model allele frequencies for poliploidy genomes and not only for the diploid case.

## Poliploid data simulations
`simulationScript.sh`: simulate poliploidy data: set the ploidy numbers, haploid depths and number of individuals editing the options in the file `simulationScript.sh`.

Options in the script are the following as it follows
* `ploidy`: vector with the desired sequence of ploidy numbers  
* `depth`: vector with the desired haploid depths
* `sites`: number of sites for each ploidy number
* `samples`: number of individuals to simulate
* `BASENAME`: base name for the output files 
* `FOLDER`: folder where to save the outputs
* In the ANGSD command called in the script, see the explanation for the `minInd` and `minIndDepth` filtering options.

Output: 

In you folder `$FOLDER`, there are four groups of data used in ploidy inference, one for each combination of depth and number of individuals (4 combinations in this example) in the following formats and with same basename:
* `.genolikes` where the columns represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 6 (7 columns)
* `.mafs` output from ANGSD calculating allele frequencies at each locus. Each column represents: chromosome, site number, major allele, minor allele, reference allele, estimated frequency, number of individuals with data.
* `.mpileup` mpileup file of the simulated genome
* `.initials` initial estimate for ploidy levels parameters. Each two lines represent the alpha and beta parameters of negative binomial distributions modelled on observed depths of one individual.

#### Example:

simulate a genome called `poliploidyGenome` with sequence of ploidy numbers 2-5-4-2 for two different haploid depths, 3X and 8X, and two different amount of individuals, 5 and 10. Consider `1000` simulated loci for each ploidy level. 
Inside the script, edit the options as it follows
* `ploidy=(2 5 4 2)` 
* `depth=(3 8)`
* `sites=1000`
* `samples=(5 10)`
* `BASENAME=poliploidyGenome`

Run the script:
```Shell
bash simulationScript.sh
```



## Infer ploidy levels from simulated data

Create a file `basenames.filelist` with the list of basenames for the simulated genomes. 
Continuing the example from simulated data, the file would for example contain the following basenames:
```
poliploidyGenome.DP3.NIND10
poliploidyGenome.DP3.NIND5
poliploidyGenome.DP8.NIND10
poliploidyGenome.DP8.NIND5
```

Run the R script from the shell



## Infer ploidy levels from `mpileup` file







