#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
                                        
##### INPUTS #####
baseName=args[1]
#lth=as.numeric(args[2])
NIND=as.numeric(args[2])
truePl=as.numeric(args[3])
#totLth=as.numeric(args[4])

simFolder = './ploidy_ref/'

print("Reading reference simulations from the simulated data at depth 50X")




minorFFiltered = list()
for(i in 1:5){
    fileName = paste(simFolder, '/ploidy_', i, '/sim.DP50.genolikes', sep='')
    cat("\t",fileName,"\n")
    GL <- fread(input=fileName,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE,select=c(2,5,8,9),col.names=c("sites","tot","major","minor"))
    
    cat(fileName,"\n")
    sites <- unlist(GL[,'sites'])
    ##for(NINT in 0:(ceiling(totLth/lth) - 1)){

    sites <- unlist(GL[,'sites'])
    minor <- unlist(GL[,'minor'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
    major <- unlist(GL[,'major'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
    tot <- major+minor
    freqs <- minor/tot
        
    sitesSNP <- (freqs>.1 & freqs<.9)

    if(i==1)
	minorFFiltered[[i]] <- freqs
    else
	minorFFiltered[[i]] <- freqs[sitesSNP & !is.na(sitesSNP)]
    print("simulation")
    print(i)
}





fileName = paste(baseName,".genolikes",sep="")
GL <- fread(input=fileName,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE,select=c(2,5,8,9),col.names=c("sites","tot","major","minor"))
    
cat(fileName,"\n",sep="")
sites <- unlist(GL[,'sites'])
chosen=c()
##for(NINT in 0:(ceiling(totLth/lth) - 1)){

sites <- unlist(GL[,'sites'])
minor <- unlist(GL[,'minor'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
major <- unlist(GL[,'major'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
tot <- major+minor
tot <- tot#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
minorF <- minor/tot
minorFFiltered_data = list()

freqs <- minorF
sitesSNP <- (freqs>.1 & freqs<.9)
freqsSNP <- freqs[sitesSNP]

for(i in 1:NIND){
    cat("Individual ",i," denoised. Ploidy estimated by ploidyNGS-like KStest: ")
    minorInd <- minorF[seq(i,length(minor),NIND)]
    minorIndSNP <- (minorInd>.1 & minorInd<.9)
    minorFFiltered_data[[i]] <- minorInd[minorIndSNP & !is.na(minorIndSNP)]
    if(is.null(length(minorFFiltered_data[[i]])) | length(minorFFiltered_data[[i]])<100)
	minorFFiltered_data[[i]] <- minorInd[!is.na(minorIndSNP)]
    
    chosen[i]=0
    dist=c(1,1,1,1,1)
    for(pp in 1:5){
    	   x = minorFFiltered_data[[i]] #data
	   y = minorFFiltered[[pp]] #simulated data
	   #print(x)
	   #print(y)
    	   dist[pp] = try({ ks.test(x, y)$statistic })
    }

    
    chosen[i] = which(dist==min(dist))[1]
    print(length(which(dist==min(dist))))
    
    if(length(which(dist==min(dist)))>1)
	chosen[i]=0
    cat(chosen[i],"\n")
    print(dist)
}

print(chosen)
rate = (sum(chosen==truePl)/length(chosen))
cat(rate,file=paste(baseName,".ploidyNGS",sep=""))
cat(baseName," done\n",sep="")
##}

