# Hidden Markov Ploidy
Tool for inferring ploidy levels, testing for aneuploidy and other stuff.
Calculating its allele frequencies and genotype likelihoods requires to download the following tool
* [ANGSD](https://github.com/ANGSD/angsd)
* [NGSPOLY](https://github.com/ImperialCollegeLondon/ngsJulia/tree/master/ngsPoly)
* [JULIA vers >= 0.4.7](https://julialang.org/downloads/)

## Simulate data

Simulate poliploidy data: set the ploidy numbers, haploid depths and number of individuals editing the options in the file `simulationScript.sh`.

### Example:
simulate a genome called `poliploidyGenome` with sequence of ploidy numbers 2-5-4-2 for two different haploid depths, 3X and 8X, and two different amount of individuals, 5 and 10. Consider `1000` simulated loci for each ploidy level. 
Inside the script, edit the options as it follows
* `ploidy=(2 5 4 2)` 
* `depth=(3 8)`
* `sites=1000`
* `samples=(5 10)`
* BASENAME=`poliploidyGenome`

Note: Eventually. change the variable `FOLDER` to save the output in a specific folder that must exist when the script is executed.

Run the script:
 bash simulationScript.sh



## Script examples
