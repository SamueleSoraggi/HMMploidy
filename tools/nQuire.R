#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

                                        
library(data.table)
library(Rcpp)

 
cppFunction('NumericVector alleleFrequencies(NumericVector major,NumericVector minor, int nInd, int minInd, double eps){

int totCountsNorm = 0;
int sites = major.size()/nInd;
NumericVector out( sites );
int indWithData;
NumericVector pis( nInd );
NumericVector wis( nInd );
int ni = 0;
int nt = 0;
int normC;

for(int s=0;s<sites;s++){
    totCountsNorm = 0;
    indWithData = nInd;

    normC = 0;
    for(int i=0;i<nInd;i++){
        totCountsNorm += major[s*nInd+i] + minor[s*nInd+i];
    }
    //if the site is variable
    for(int i=0;i<nInd;i++){
       ni = minor[s*nInd+i];
       nt = major[s*nInd+i] + minor[s*nInd+i]; 
       if(nt==0){//if we dont have any reads for individual i
         indWithData--;
         pis[i] = 0;
         wis[i] = 0;
         continue;
       }
       pis[i] = (ni-eps*nt)/(nt*(1-2*eps));
       wis[i] = (double)nt/totCountsNorm; //weights infinite ploidy
     }

  if(indWithData < minInd){
     out[s] = -1;
  }
  else{
     out[s] = 0;
     for(int i=0;i<nInd;i++)
        out[s] += wis[i]*pis[i];     
  }

}
return(out);

}'
)


#p=runif(10000)
#g0=runif(10000)
#g1=rnorm(10000,0.25,.05)
#g2=rnorm(10000,0.5,.05)
#g3=rnorm(10000,0.75,.05)
#prop=c(.1,.3,.3,.3)
#inter=cumsum(prop)
#res=(p<inter[1])*g0 + (p>inter[1] & p<inter[2])*g1 + (p>inter[2] & p<inter[3])*g2 + (p>inter[3] & p<inter[4])*g3


##function for gmm. if unifNoise is TRUE, then the first component is Uniform(0,1) noise, and the other gaussian components. k=number of gaussians.
##alpha0 must includes a proportion  
gmmuEM <- function( x, k=NULL, alpha0=NULL, mu0=NULL, sigma0=NULL, unifNoise=FALSE, fixMu=FALSE, fixSigma=FALSE, fixAlpha=FALSE, fixNoise=FALSE, printout=FALSE, maxIter=1000, tries=100, tol=1e-3, seed=NULL ){
    
    
    L <- length(x)
    if(!is.null(seed))
        set.seed(seed)

    if(is.null(mu0)) #mu0 <- sample(x, k, replace=FALSE, prob=rep(1/L,L))
        mu0 <- quantile(x,seq(0,1,length.out=k+2))[-c(1,k+2)]
    if(is.null(sigma0)){
        #mixtools-like method for binning. Does not help a lot, but pretty stable.
        x2 = sort(x)
        x2.bin = list()
        for (j in 1:k) 
            x2.bin[[j]] <- x2[max(1, floor((j - 1) * L/k)):ceiling(j * L/k)]
        s.hyp = as.vector(sapply(x2.bin, sd))
        if (any(s.hyp == 0)) 
            s.hyp[which(s.hyp == 0)] = runif(sum(s.hyp == 0), 0, sd(x2))
        s = 1/rexp(k, rate = s.hyp)
        
        sigma0 = s^2
    }
    if(is.null(alpha0))
        alpha0 <- rep(1/(k+unifNoise),k+unifNoise)

    sigma <- sigma0
    mu <- mu0
    alpha <- alpha0
    
    if(unifNoise){
        k=k+1 
        sigma <- c(0,sigma0)
        mu <- c(0,mu0)
    }
    
    
    cycleFlag <- TRUE

    nu <- matrix(0, nrow=L, ncol=k)
    if(unifNoise)
        nu[,1] <- alpha[1]*dunif(x) 
    else
        nu[,1] <- alpha[1]*dnorm(x,mu[1],sqrt(sigma[1]))
    
    for(j in 2:k)
        nu[,j] <- alpha[j]*dnorm(x,mu[j],sqrt(sigma[j]))
    
    llkNew <- sum(log(rowSums( nu )))
    
    #normalize
    normConst <- rowSums(nu)
    for(j in 1:k)
        nu[,j] <- nu[,j]/normConst
    
    mixSize <- colSums(nu) #mixture size
    stepCounter <- 0
  
    while(cycleFlag){

        llkOld <- llkNew
        stepCounter <- stepCounter+1

        if(any(sigma<1e-5)){#check for crazy sigmas
            sigmaBuffer <- sigma
            if(unifNoise)
                sigmaBuffer <- sigmaBuffer[-1]
            sigmaBuffer[sigmaBuffer < 1e-5] <- 1e-5
            if(unifNoise)
                sigma <- c(0,sigmaBuffer)
            else
                sigma <- sigmaBuffer
        }
    
        #calc responsibilities
        if(unifNoise)
            nu[,1] <- alpha[1]*dunif(x) 
        else
            nu[,1] <- alpha[1]*dnorm(x,mu[1],sqrt(sigma[1]))
            
        for(j in 2:k)
            nu[,j] <- alpha[j]*dnorm(x,mu[j],sqrt(sigma[j]))
        
        llkOld <- sum(log(rowSums( nu )))
        
        #normalize
        normConst <- rowSums(nu)
        for(j in 1:k)
            nu[,j] <- nu[,j]/normConst

        mixSize <- colSums(nu) #mixture size

        #estimate new parameters
        sigmaOld <- sigma
        muOld <- mu
        alphaOld <- alpha

        if(unifNoise){
            if(!fixNoise)
                alpha[1] <- mixSize[1]/L
            if(!fixNoise & fixAlpha)
                alpha <- alpha/sum(alpha)
        }
        else{
            if(!fixAlpha)
                alpha[1] <- mixSize[1]/L
            if(!fixMu)
                mu[1] <-  sum(nu[, 1]*x)/mixSize[1]
            if(!fixSigma)
                sigma[1] <- sum(nu[, 1]*(x-mu[1])^2)/mixSize[1]
        }

        for(j in 2:k){
            if(!fixAlpha)
                alpha[j] <- mixSize[j]/L
            if(!fixMu)
                mu[j] <-  sum(nu[, j]*x)/mixSize[j]
            if(!fixSigma)
                sigma[j] <- sum(nu[, j]*(x-mu[j])^2)/mixSize[j]
        }

        #calc responsibilities
        if(unifNoise)
            nu[,1] <- alpha[1]*dunif(x)
        else
            nu[,1] <- alpha[1]*dnorm(x,mu[1],sqrt(sigma[1]))
            
        for(j in 2:k)
            nu[,j] <- alpha[j]*dnorm(x,mu[j],sqrt(sigma[j]))
        
        llkNew <- sum(log(rowSums( nu )))
        if(unifNoise)
            llkNoNoise <- sum(log(rowSums( as.matrix(nu[,-1]) )))
        
        #normalize
        normConst <- rowSums(nu)
        for(j in 1:k)
            nu[,j] <- nu[,j]/normConst

        mixSize <- colSums(nu) #mixture size

        epsilon <- abs(llkNew-llkOld)/abs(llkOld)
        if(epsilon < tol){
            cycleFlag = FALSE
            if(printout)
                cat("Convergence in ",stepCounter," iterations\n",s)
        }
        if(stepCounter==maxIter){
            cycleFlag = FALSE
            if(printout)
                cat("No convergence in ",stepCounter," iterations\n",s)
        }

        if(printout)
            cat("---step ",stepCounter,"\nepsilon: ",epsilon,"\nalpha: ",alpha,"\nmu: ",mu,"\nsigma",sigma,"\nllkOld",llkOld,"\nllkNew",llkNew,"\n")
    }
    if(!unifNoise)
        llkNoNoise <- NULL
    
    return(list(alpha=alpha,llk=llkNew,mu=mu,variance=sigma,stdev=sqrt(sigma),alpha0=alpha0,mu0=mu0,sigma0=sigma0,nu=nu,llkNoNoise=llkNoNoise))
}




baseName=args[1]
#lth=as.numeric(args[2])
NIND=as.numeric(args[2])
truePl=as.numeric(args[3])
#totLth=as.numeric(args[4])

fileName = paste(baseName,".genolikes",sep="")
GL <- fread(input=fileName,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE,select=c(2,5,8,9),col.names=c("sites","tot","major","minor"))
    
cat(fileName,"\n",sep="")
sites <- unlist(GL[,'sites'])
chosen=c()
    #for(NINT in 0:(ceiling(totLth/lth) - 1)){

sites <- unlist(GL[,'sites'])
minor <- unlist(GL[,'minor'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
major <- unlist(GL[,'major'])#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
tot <- major+minor
tot <- tot#[ sites>NINT*lth & sites<=(NINT+1)*lth ]
minorF <- minor/tot
        
freqs <- alleleFrequencies(major,minor,nInd=as.numeric(NIND),minInd=1,eps=0)
sitesSNP <- (freqs>.1 & freqs<.9)
freqsSNP <- freqs[sitesSNP]

for(i in 1:NIND){
    cat("Individual ",i," denoised. Ploidy estimated by nQuire: ")
    minorInd <- minorF[seq(i,length(minor),NIND)]
    minorIndSNP <- minorInd[sitesSNP]
    minorIndFiltered <- minorIndSNP[minorIndSNP>.1 & minorIndSNP<.9 & !is.na(minorIndSNP)]

    ###
    free2=try( gmmuEM(x=minorIndFiltered, k=1, unifNoise=TRUE, mu0=c(.5), alpha0=c(.1,.9), fixMu=TRUE, fixAlpha=TRUE), TRUE )
    if(class(free2)=="try-error")
        free2$'llkNoNoise'= +Inf
    free3=try( gmmuEM(x=minorIndFiltered, k=2, unifNoise=TRUE, mu0=c(.33,.66), alpha0=c(.1,.45,.45), fixMu=TRUE, fixAlpha=TRUE), TRUE )
    if(class(free3)=="try-error")
        free3$'llkNoNoise'= +Inf
    free4=try( gmmuEM(x=minorIndFiltered, k=3, unifNoise=TRUE, mu0=c(.25,.5,.75), alpha0=c(.1,.3,.3,.3), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free4)=="try-error")
        free4$'llkNoNoise'= +Inf
    free5=try( gmmuEM(x=minorIndFiltered, k=4, unifNoise=TRUE, mu0=c(.2,.4,.6,.8), alpha0=c(.1,.9/4,.9/4,.9/4,.9/4), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free5)=="try-error")
        free5$'llkNoNoise'= +Inf
    free6=try( gmmuEM(x=minorIndFiltered, k=5, unifNoise=TRUE, mu0=c(.16,.33,.6,.66,.83), alpha0=c(.1,.9/5,.9/5,.9/5,.9/5), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free6)=="try-error")
        free6$'llkNoNoise'= +Inf
            
    freeALL=try( gmmuEM(x=minorIndFiltered, k=4, unifNoise=TRUE, mu0=c(.2,.4,.6,.8) ), TRUE )
    if(class(freeALL)=="try-error")
        freeALL$'llkNoNoise'= +Inf
            
    llkFixed <- c(free2$llkNoNoise,free3$llkNoNoise,free4$llkNoNoise,free5$llkNoNoise)
        
    chosen[i]=0
    if( !is.nan(freeALL$llkNoNoise) & !is.na(freeALL$llkNoNoise) & !is.infinite(freeALL$llkNoNoise) & (class(freeALL)!="try-error") & sum(is.infinite(llkFixed) | is.nan(llkFixed) | is.na(llkFixed)) != 4 ){
        diffLlk <- abs(llkFixed - freeALL$llkNoNoise)
        #print(diffLlk)
        chosen[i] <- which.min(diffLlk)+1}
    cat(chosen[i],"\n")
}

rate = (sum(chosen==truePl)/length(chosen))
cat(rate,file=paste(baseName,".denoised.gmmu",sep=""))
cat(baseName," done\n",sep="")
#}



#################################
#################################
#################################


for(i in 1:NIND){
    cat("Individual ",i," denoised. Ploidy estimated by nQuire: ")
    minorInd <- minorF[seq(i,length(minor),NIND)]
    minorIndSNP <- minorInd[sitesSNP]
    minorIndFiltered <- minorIndSNP[minorIndSNP>.1 & minorIndSNP<.9 & !is.na(minorIndSNP)]

    ###
    free2=try( gmmuEM(x=minorIndFiltered, k=1, unifNoise=FALSE, mu0=c(.5), alpha0=c(.1,.9), fixMu=TRUE, fixAlpha=TRUE), TRUE )
    if(class(free2)=="try-error")
        free2$'llk'= +Inf
    free3=try( gmmuEM(x=minorIndFiltered, k=2, unifNoise=FALSE, mu0=c(.33,.66), alpha0=c(.1,.45,.45), fixMu=TRUE, fixAlpha=TRUE), TRUE )
    if(class(free3)=="try-error")
        free3$'llk'= +Inf
    free4=try( gmmuEM(x=minorIndFiltered, k=3, unifNoise=FALSE, mu0=c(.25,.5,.75), alpha0=c(.1,.3,.3,.3), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free4)=="try-error")
        free4$'llk'= +Inf
    free5=try( gmmuEM(x=minorIndFiltered, k=4, unifNoise=FALSE, mu0=c(.2,.4,.6,.8), alpha0=c(.1,.9/4,.9/4,.9/4,.9/4), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free5)=="try-error")
        free5$'llk'= +Inf
    free6=try( gmmuEM(x=minorIndFiltered, k=5, unifNoise=FALSE, mu0=c(.16,.33,.6,.66,.83), alpha0=c(.1,.9/5,.9/5,.9/5,.9/5), fixMu=TRUE,fixAlpha=TRUE), TRUE )
    if(class(free6)=="try-error")
        free6$'llk'= +Inf
            
    freeALL=try( gmmuEM(x=minorIndFiltered, k=4, unifNoise=FALSE, mu0=c(.2,.4,.6,.8) ), TRUE )
    if(class(freeALL)=="try-error")
        freeALL$'llk'= +Inf
            
    llkFixed <- c(free2$llk,free3$llk,free4$llk,free5$llk)
        
    chosen[i]=0
    if( !is.nan(freeALL$llk) & !is.na(freeALL$llk) & !is.infinite(freeALL$llk) & (class(freeALL)!="try-error") & sum(is.infinite(llkFixed) | is.nan(llkFixed) | is.na(llkFixed)) != 4 ){
        diffLlk <- abs(llkFixed - freeALL$llk)
        #print(diffLlk)
        chosen[i] <- which.min(diffLlk)+1}
    cat(chosen[i],"\n")
}

rate = (sum(chosen==truePl)/length(chosen))
cat(rate,file=paste(baseName,".noisy.gmmu",sep=""))
cat(baseName," done\n",sep="")
