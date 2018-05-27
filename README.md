# Hidden Markov Ploidy
Tool for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the following tool
* [ANGSD](https://github.com/ANGSD/angsd)
* [NGSPOLY](https://github.com/ImperialCollegeLondon/ngsJulia/tree/master/ngsPoly)
* [JULIA vers >= 0.4.7](https://julialang.org/downloads/)

Once ANGSD is downloaded, substitute the file abcFreq.cpp with the one in this github repository, and thereafter compile ANGSD. This allows to model allele frequencies for poliploidy genomes and not only for the diploid case.

## Simulate data

Simulate poliploidy data: set the ploidy numbers, haploid depths and number of individuals editing the options in the file `simulationScript.sh`.

### Example:

simulate a genome called `poliploidyGenome` with sequence of ploidy numbers 2-5-4-2 for two different haploid depths, 3X and 8X, and two different amount of individuals, 5 and 10. Consider `1000` simulated loci for each ploidy level. 
Inside the script, edit the options as it follows
* `ploidy=(2 5 4 2)` 
* `depth=(3 8)`
* `sites=1000`
* `samples=(5 10)`
* `BASENAME=poliploidyGenome`

Note: Eventually. change the variable `FOLDER` to save the output in a specific folder that must exist when the script is executed.

Run the script:
```shell
bash simulationScript.sh
```

In you folder `$FOLDER`, there are four groups of data used in ploidy inference, one for each combination of depth and number of individuals (4 combinations in this example) in the following formats and with same basename:
* `.genolikes` where the columns represent: chromosome name, site number, individual number, ref.allele, site coverage, major allele, minor allele, genotype likelihoods at ploidy 1 (2 columns), genotype likelihoods at ploidy 2 (3 columns), ..., genotype likelihoods at ploidy 6 (7 columns)
* `.mafs` output from ANGSD calculating allele frequencies at each locus. Each column represents: chromosome, site number, major allele, minor allele, reference allele, estimated frequency, number of individuals with data.
* `.mpileup` mpileup file of the simulated genome



## Infer ploidy levels from simulated data



## Infer ploidy levels from `mpileup` file







