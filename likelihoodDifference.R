library(pracma) #strcmp function library
library(data.table) #fread function library
library(Rcpp)
set.seed(041205)

##############################
### script input interface ###
##############################

l<-commandArgs(TRUE)
getArgs<-function(x,l)
    unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
    if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
        cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
        q("no")
    }
    arguments<-list()
    for(a in names(args))
        arguments[[a]]<-getArgs(a,l)
    
    if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
        cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
        q("no")
    }
    for(a in names(args))
        if(is.null(arguments[[a]]))
            arguments[[a]]<-args[[match(a,names(args))]]
    
    arguments
}

print.args<-function(args,des){
    if(missing(des)){
        des<-as.list(rep("",length(args)))
        names(des)<-names(args)
    }
    cat("->  Needed arguments:\n")
    mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
    cat("->  Optional arguments (defaults):\n")
    mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
    q("no")
}


## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(file = NA, #single basename of file to analize (does not need the list 'filelist' for multiple files)
           fileList = NA, #list of basenames for GUNZIPPED .genolike, .mafs and .par files
           wind = NA, #size of window for depth and genotype likelihoods. we work on a chromosome-basis.
           minInd = 1, #min ind having reads 
           CNVInd = NA, #which Individual to consider (one at the time for now)
           quantileTrim ="0,1", #quantiles for trimming
           maxPloidy=6,
           eps = .0005 #effect of sequencing and mapping error
           )


des<-list(fileList="[string] list of .genolike and phat.mafs and eventual .windows files",
          wind="[integer] Size of window for depth and genotype likelihoods (NA)",
          minInd="[integer] min Nr of individuals per locus having data (1)",
          CNVInd ="[integers] which Individuals might have CNV",
          beta="[numerics] beta parameters comma separated (NA, read from .par file)",
          quantileTrim="[integers] comma-separated quantiles for trimming (0,1)",
          maxploidy="max ploidy",
          eps="sequencing/mapping error rate "
          )



##get arguments and add to workspace
##do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
    cat(" Arguments: output prefix\n")
    q("no")
}


##read file names from the prefix list and create inputs/outputs
if(is.na(file) & is.na(fileList)){
    cat("You MUST input either file= or fileList= when using the R script\n")
    q("no")
}
if(!is.na(file))
    filez <- file
if(is.na(file))
    filez <- unlist( read.table(fileList, header=FALSE, as.is=T)  )

fileVector <- c(); hmmVector <- c(); outPdf <- c(); outTxt <- c();
BASENAMEFILE <- c();
for(i in 1:length(filez)){
    fileVector[i] <- paste(filez[i],".genolikes",sep="")
    hmmVector[i] <- paste(filez[i],".hiddenMarkovPloidy",sep="")
    splittedName <- unlist(strsplit(filez[i],split="/"))
    BASENAMEFILE[i] <- splittedName[length(splittedName)]
    outPdf[i] <- paste(filez[i],".CNV.pdf",sep="")
    outTxt[i] <- paste(filez[i],".hiddenMarkovCNV",sep="")
}


isNumericCNVInd <- all(!is.na(CNVInd)) #check for choice of individuals

##individuals chosen for analysis (one by one at the moment)
##if chosenInd id NA it will be assigned as all individuals later
if(isNumericCNVInd)
    CNVInd <- eval( parse( text=paste("c(",CNVInd,")",sep="") ) )




wind <- as.numeric(wind)
minInd <- as.numeric(minInd)

##print Rscript input
cat("----------\nfileList: ", fileList, " wind: ", wind," minInd: ", minInd, " CNVInd: ", CNVInd ," quantileTrim: ", quantileTrim, " eps: ", eps,  "\n-----------\n" )


cppFunction('NumericVector alleleFrequencies(NumericVector major, NumericVector minor, int nInd, int minInd, double eps){

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


llkCalc <- function(count,delta,TRANS,alpha,beta,genolike){

    N <- nrow(TRANS)
    Total <- dim(count)[1]
    forwrd <- matrix(0, Total, N)
    forwrd2 <- matrix(0, Total, N)
    bckwrd <- matrix(0, Total, N)
    bckwrd2 <- matrix(0, Total, N)
    dens <- matrix(0, Total, N)
    scale <- rbind(1, matrix(0, nrow=Total-1))
    scale2 <- rbind(1, matrix(0, nrow=Total-1))
    
    cm <- max(count)
    if(cm > 50000){
        dnorm <- as.matrix(lgamma(count + 1))
    } else {
        tmp <- cumsum(rbind(0, log(as.matrix(1:max(count)))))
        dnorm <- as.matrix(tmp[count+1])
    }
    
    
    densLog <- matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha)) - count %*% log(1+beta) + lgamma(count %*% matrix(1, ncol=N) + matrix(1, nrow=Total) %*% alpha) - dnorm %*% matrix(1, ncol=N) 
    
    dens2 <- matrix(0,nrow=nrow(densLog),ncol=dim(genolike)[2])
    for(ii in 1:dim(genolike)[2]) dens2[,ii] <- densLog[,ii] + genolike[,ii]
    
    #dens <- exp( densLog )
    dens2 <- exp( dens2 )

    #dens <- apply(dens, 2, function(x) {x[x==0 | is.na(x) | is.nan(x)] <- .Machine$double.xmin; x})        
    dens2 <- apply(dens2, 2, function(x) {x[x==0 | is.na(x) | is.nan(x)] <- .Machine$double.xmin; x})
    
    
    #forwrd[1,] <- delta*dens[1,];
    forwrd2[1,] <- delta*dens2[1,]
		
    for(t in 2:Total){
                                        #forwrd[t,] <- (forwrd[t-1,] %*% TRANS) * dens[t,]
        forwrd2[t,] <- (forwrd2[t-1,] %*% TRANS) * dens2[t,]
        #scale[t] <- sum(forwrd[t,])
        scale2[t] <- sum(forwrd2[t,])
        #forwrd[t,] <- forwrd[t,] / scale[t]
        forwrd2[t,] <- forwrd2[t,] / scale2[t]
        if(is.nan(forwrd2[t,]))
            forwrd2[t,] <- .Machine$double.xmin
    }
				
    llk <- log(sum(forwrd2[Total,])) + sum(log(scale2))
    
    
    return(llk)
}

##logarithmic normalization of a vector
logRescale <- function(v){
    L <- length(v)
    m <- max(v)
    w <- which.max(v)
    diffVec <- v[-w] - m
    if(any(diffVec < -700)){
        idx <- which(diffVec < -700)
        tooSmall <- diffVec[idx]
        sortVec <- sort(tooSmall, index.return=TRUE)
        rescaled <- seq(-700, -800, length.out=length(tooSmall))
        tooSmall[sortVec$ix] <- rescaled
        diffVec[idx] <- tooSmall
    }
    res <- m + log( 1 + sum( exp( diffVec ) ) )
    return( v - res )
}



##sum of values in a windows. When lociSNP=loci all values in the windows are used.
##avg=TRUE performs average instead of sum. ws=window size and dp=vector of data.
sumGeno <- function(dp,ws,loci,lociSNP=loci,findSNP=1:length(loci),avg=FALSE){   
    L <- length(dp)    
    S <- seq(loci[1],loci[length(loci)],ws)
    S <- c(S, loci[length(loci)] )
    
    res <- sapply(1:(length(S)-1), function(ll){
        if(ll<(length(S)-1))
            idx <- which(lociSNP>=S[ll] & lociSNP<S[ll+1])
        if(ll==(length(S)-1))
            idx <- which(lociSNP>=S[ll] & lociSNP<=S[ll+1])
        if(length(idx)==0)
            return(c())
        v <- dp[findSNP[idx]]
        if(!avg)
            return( sum( v ) )
        if(avg)
            return( mean(v) )
    })
    return(unlist(res))
}

sumGenoAll <- function(dp,ws=1,loci,lociSNP=loci,avg=FALSE){   
    L <- length(dp)    
    S <- c( seq(loci[1],loci[length(loci)-1],ws), loci[length(loci)] )
    res <- sapply(1:(length(S)-1), function(ll){
        idx <- which(lociSNP>S[ll] & lociSNP<S[ll+1])
        v <- dp[idx]
        if(!avg)
            return( sum( v ) )
        if(avg)
            return( mean(v) )
    })
    return(res)
}

##sum of logarithm on rows and cols of a matrix
rowSumsLog <- function(X){
    res <- apply(X, 1, function(t){
        m <- max(t)
        w <- which.max(t)
        diffVec <- t[-w] - m
        diffVec[diffVec < -700] <- -700
        diffVec[diffVec == 0] <- -0.0001
        return( m + log( 1 + sum( exp( diffVec ) ) ) )
    })
    return( res )
}

colSumsLog <- function(X){
    res <- apply(X, 2, function(t){
        m <- max(t)
        w <- which.max(t)
        diffVec <- t[-w] - m
        diffVec[diffVec < -700] <- -700
        diffVec[diffVec == 0] <- -0.0001
        return( m + log( 1 + sum( exp( diffVec ) ) ) )
    })
    return( res )
}
   

##Likelihood of f=data vector given genotype. gl=genotype likelihoods vector. h=inbreeding coefficient.
pGenoData <- function(f,winL,gl,nInd=1,findSNP=1:sum(winL),h=0){   
    y = ncol(gl)-1
    Lf = sum(winL)
    winIdx=cumsum(c(0,winL))

    X <- rep(0,Lf)
    for(l in 1:length(f)){

        freq <- f[l]
        idx <- (winIdx[l]+1):(winIdx[l+1])
        p <- dbinom(0:y,y,freq,log=TRUE)
        p[is.infinite(p)]=-1000
        #glSum <- gl[idx,] + matrix( rep( p, nInd*length(idx) ), nrow=nInd*length(idx), byrow=T )
        glSum <- apply( gl[idx,], 1, function(r) r + p )
        X[(winIdx[l]+1):(winIdx[l+1])] <- colSumsLog( glSum )
        #print(X[(winIdx[l]+1):(winIdx[l]+2)])
    }
    return( X )    
}

pGenoDataAll <- function(f,gl,nInd=1,h=0){   
    y = ncol(gl)-1
    Lf = length(f)
    nInd = dim(gl)[1] / Lf
    
    matrix( rep(  dbinom(0:y,y,f,log=TRUE), nInd ), nrow=nInd, byrow=T )

    X <- c()
    for(l in 1:length(f)){

        freq <- f[l]
        idx <- ((l-1)*nInd+1):(l*nInd)
        p <- dbinom(0:y,y,freq,log=TRUE)
        p[is.infinite(p)]=-1000
        glSum <- gl[idx,] + matrix( rep(  dbinom(0:y,y,freq,log=TRUE), nInd ), nrow=nInd, byrow=T )
        X[l] <- sum( rowSumsLog( glSum ) ) 
        
    }
    return( X )    
}

pGenoDataSingle <- function(f,gl,h=0){  #use one individual at a time
    y = ncol(gl)-1
    Lf = dim(gl)[1]
    fVector=rep(f,Lf)
    nInd = 1
    
    matrix( rep(  dbinom(0:y,y,fVector,log=TRUE), nInd ), nrow=nInd, byrow=T )
    
    X <- c()
    for(l in 1:Lf){

        freq <- f[l]
        idx <- ((l-1)*nInd+1):(l*nInd)
        p <- dbinom(0:y,y,freq,log=TRUE)
        p[is.infinite(p)]=-1000
        glSum <- gl[idx,] + matrix( rep(  dbinom(0:y,y,freq,log=TRUE), nInd ), nrow=nInd, byrow=T )
        X[l] <- sum( rowSumsLog( glSum ) ) 
        
    }
    #find sites hvor f passer bedst - ordering?
    return( X )    
}


##read genotype likelihood at a certain site, for a certain dataset GL,
##given ploidy and number of individuals
readGL <- function(site,ploidy,nInd=1,GL){
    col = cumsum( c(1:(ploidy+1)) )
    if(nInd > 1){
        res <- as.vector( sapply(site, function(x) return( c( ((x-1)*nInd+1):(x*nInd) ) ) ) )
        site <- res
    }   
    if(length(site)==1){
        return(   GL[site, c(col[ploidy]:(col[ploidy+1]-1))]   )}
    else{
        return(  GL[site, c(col[ploidy]:(col[ploidy+1]-1))]   )}
}      


#need genolikes >>> need frequencies. need to import counts as well.


for(i in 1:length(fileVector)){ #loop over input files
    
    cat("==> Analyze ", filez[i], "\n",sep="")
    ##read in the data from .mafs and .genolikes files
    GL <- fread(input=fileVector[i],sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)
    com <- paste("sed -n '6~8p' ",hmmVector[i], sep="")
    LLKclean <- unlist(fread(com,sep="\n",showProgress=TRUE,header=FALSE,data.table=FALSE))

    sortLLK <- sort(LLKclean,index.return=TRUE,decreasing=TRUE)
    for(j in 1:length(LLKclean)){
        if( !any(sortLLK$ix[j] == CNVInd) ){
            bestIdx <- j
            bestLLK <- sortLLK$x[j]
            break
        }
    }
    cat("\tBest LlK ",sortLLK$x[j]," in individual ",sortLLK$ix[j],"\n",sep="")
    
    com <- paste("sed -n '", (bestIdx-1)*8+4  , "p' ",hmmVector[i], sep="")
    alphaBest <- na.omit(unlist(fread(com,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)))
    com <- paste("sed -n '", (bestIdx-1)*8+5  , "p' ",hmmVector[i], sep="")
    betaBest <- na.omit(unlist(fread(com,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)))
    
    com <- paste("sed -n '", (bestIdx-1)*8+2  , "p' ",hmmVector[i], sep="")
    deltaBest <- na.omit(unlist(fread(com,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)))
    com <- paste("sed -n '", (bestIdx-1)*8+3  , "p' ",hmmVector[i], sep="")
    TRANSbest <- na.omit(unlist(fread(com,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)))
    TRANSbest <- matrix(TRANSbest, ncol=length(TRANSbest)/2, nrow=length(TRANSbest)/2)

    com <- paste("sed -n '", (bestIdx-1)*8+7  , "p' ",hmmVector[i], sep="")
    statesBest <- na.omit(unlist(fread(com,sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)))

    if(length(statesBest)<length(deltaBest)){
        L <- length(deltaBest)
        dL <- length(deltaBest)-length(statesBest)
        sortDelta <- sort(deltaBest,index.return=TRUE,decreasing=TRUE)
        rmIdx <- sortDelta$ix[ seq(L-dL+1, L, 1) ]
        deltaBest <- deltaBest[-c(rmIdx)]
        TRANSbest <- as.matrix( TRANSbest[-c(rmIdx),-c(rmIdx)] )
    }

    
    rowsGL <- dim(GL)[1]
    nInd <- length( unique( GL[,3] ) )
    sites <- unique( GL[ ,2] )    
    DP <- GL[ ,5]

    #calculate allele frequencies
    majorReads <- GL[,8]
    minorReads <- GL[,9]
    freqs <- alleleFrequencies(majorReads,minorReads,nInd,minInd,eps) #frequencies (used for SNPs on a single individual)
    
    GL <- GL[ ,-c(1:9)]
    

    ###############################
    ## begin of FOR loop to      ##
    ## read one genome at a time ##
    ###############################

    diffLLK <- rep(NA,length(CNVInd))
    
    fileCounter = 1
    ##open pdf plot connection
    pdf(outPdf[i])
    title <- sprintf("Difference of LLK calculated on the best scoring model\n%s",BASENAMEFILE[i])
    
    for(whichInd in 1:nInd){ #loop over individuals
        
        idxSingle <- seq(whichInd,rowsGL,nInd)
        majorSingle <- majorReads[idxSingle]
        minorSingle <- minorReads[idxSingle]
        DPsingle <- DP[idxSingle]; GLsingle <- GL[idxSingle, ] #individual depth/genolikes
    ##trim the depth at the chosen quantile
        quantiles <- eval( parse( text=paste("c(",quantileTrim,")",sep="") ) )
        q <- quantile( DPsingle, quantiles )  
        idx <- which( DPsingle<=as.numeric(q[2]) & DPsingle>=as.numeric(q[1]) )
        DPsingle <- DPsingle[idx] #individual filtered data
        GLsingle <- GLsingle[idx, ] #...""
        sitesIndiv <- sites[idx] #......""
        freqsIndiv <- freqs[idx] #......""
        majorSingle <- majorReads[idx] #......""
        minorSingle <- minorReads[idx] #......""
        idxTot = as.vector( sapply(idx, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        GLfiltered <- GL[idxTot, ] #all data filtered
        DPfiltered <- DP[idxTot] #......""

        
        ##find SNPs with thresholds .1<f<.9 and data in the individual
        findSNP <- which( freqsIndiv>.1 & freqsIndiv<.9 )
        freqsSNP <- freqsIndiv[findSNP]
        sitesSNP <- sitesIndiv[findSNP]
        totSNP <- as.vector( sapply(findSNP, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        ##DPSNP = DPsingle[totSNP] I think it is not needed

        
                                        #frequencies over windows
        #winAnalysis <- freqsSingle( majorSingle, minorSingle, ws=wind, sitesIndiv, sitesSNP, findSNP)
        #winFreq <- winAnalysis$winF
        #winLth <- winAnalysis$winL
        
        #geno2 <- matrix(0, nrow=maxPloidy, ncol=sum(winLth))
        geno2 <- matrix(0, nrow=maxPloidy, ncol=length(freqsSNP))
        for(pp in 1:maxPloidy) #change ploidy
            geno2[pp,] <- pGenoDataAll( f=freqsSNP, gl=readGL( findSNP, pp, nInd=1, GLsingle ), nInd=1 )
        #print(geno2)
        #for(pp in 1:maxPloidy) #change ploidy
            #geno2[pp,] <- pGenoData( f=winFreq, winL=winLth, gl=readGL( 1:sum(winLth), pp, nInd=1, GLsingle ), findSNP=findSNP, nInd=1 )

        ##...and per window
        #print(dim(geno2))
        #geno <- apply( geno2, 1, function(x) sumGeno(x,wind,sitesIndiv,sitesSNP,findSNP) )
        geno <- apply( geno2, 1, function(x) sumGenoAll(x,wind,sitesIndiv,sitesSNP) )
        
        DPmean <- sumGeno( DPsingle, wind, sitesIndiv, sitesIndiv, 1:length(sitesIndiv), avg=TRUE )
        ##clean from NA, NaN or infinite values
        keepSites <- apply( geno, 1, function(x) sum(is.na(x) | is.nan(x) | is.infinite(x))==0 )
        #print(which(keepSites))
        DPmean <- DPmean[which(keepSites)]
        geno <- geno[which(keepSites), ]
        #keepSites <- which( !is.na(DPmean) & !is.nan(DPmean) & !is.infinite(DPmean) )
        #print(keepSites)
        #DPmean = DPmean[keepSites]
        #geno = geno[keepSites,]
        ##rescale likelihood of the data (avoids underflow)
        genoResc <- t( apply( geno , 1, logRescale ) )
        genoResc[genoResc>-.00001]=-.00001
        ##some initial parameters
        delta=rep(1/maxPloidy,maxPloidy) #i think it is ok without prior info
        Pi0=matrix(1/maxPloidy,nrow=maxPloidy,ncol=maxPloidy) #tridiagonal makes more sense?
        count <- matrix(DPmean,ncol=1)

        #choose columns
        llkData <- llkCalc(count,deltaBest,TRANSbest,alphaBest,betaBest,as.matrix(genoResc[,statesBest]))
        diffLLK[fileCounter] <- abs(bestLLK - llkData)

        cat("\tFile: ",BASENAMEFILE[i]," individual ",whichInd," ",sep="")
        cat( "\tLLK abs diff ", diffLLK[fileCounter],sep="","\n")

        cat("File: ",BASENAMEFILE[i]," individual  ",whichInd,"\n",sep="",file=outTxt[i],append=!(fileCounter==1))
        cat( diffLLK[fileCounter],"\n",file=outTxt[i],append=TRUE)

        fileCounter=fileCounter+1        

    }

    colour <- rep("green",nInd); colour[CNVInd] <- "orange"
    
    barplot( diffLLK, col=colour, xlab="Individual", ylab="bestLLK - IndividualLLK", main=title )
    dev.off()
}
