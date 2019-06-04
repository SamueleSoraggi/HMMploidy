library(pracma) #strcmp function library
library(data.table) #fread function library
library(Rcpp)

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
           fileList = NA, #list of basenames for GUNZIPPED .genolike files
           nameList = NA, #bed file to choose chromosomes and sites
           useGeno = "yes", #use genotype likelihoods
           SNPtrim = "0.1,0.9",
           wind = NA, #size of window for depth and genotype likelihoods. we work on a chromosome-basis.
           maxPloidy = 4, #maximum ploidy. Must change with choice of potential ploidies (e.g. haploid might be excluded a priori by users)
           minInd = 1, #min ind having reads
           minDepth=1,
           maxDepth=Inf,
           chosenInd = NA, #which Individual to consider (one at the time for now)
           #Must implement option for: all ind together or one at the time ?
                                        #isSims = FALSE, #data is simulated. if TRUE fileList will also refer to file(s) with true ploidies
           truePl = NA, #true ploidies (as argument c(truePl) of an R vector)
           alpha = NA, #alpha parameters comma separated
           beta = NA, #beta parameters comma separated
           quantileTrim ="0.02,0.98", #quantiles for trimming
           outSuffix = NA,#path + prefix for output filename
           eps = .0005 #effect of sequencing and mapping error
           )

#if no argument aree given prints the need arguments and the optional ones with default
des<-list(fileList="[string] list of .genolike files",
          wind="[integer] Size of window for depth and genotype likelihoods (NA)",
          nameList = "[NA] List of names for the samples",
          SNPtrim = "[integers] comma-separated freq for SNP trimming (0.1,0.9)",
          useGeno = "yes=use genotype likelihoods, otherwise not",
          minInd="[integer] min Nr of individuals per locus having data (1)",
          minDepth="min average depth per window (1)",
          maxDepth="max average depth per window (Inf)",
          maxPloidy="[integer] Maximum ploidy allowed (6)", #have to implement case where ploidies are chosen
          chosenInd ="[integers] which Individual to consider (NA=all)",
          isSims="[bool] data is simulated using the simulation script (FALSE)",
          alpha ="[numerics] alpha parameters comma separated (NA, read from .par file)",
          beta="[numerics] beta parameters comma separated (NA, read from .par file)",
          quantileTrim="[integers] comma-separated quantiles for trimming (0,1)",
          eps="sequencing/mapping error rate ",
          outSuffix="suffix for the output file (Default: empty name)"
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

angsdVector <- c(); fileVector <- c(); outPdf <- c(); outTxt <- c();
BASENAMEFILE <- c();
for(i in 1:length(filez)){
    fileVector[i] <- paste(filez[i],".genolikes",sep="")
    angsdVector[i] <- paste(filez[i],".mafs",sep="")
    splittedName <- unlist(strsplit(filez[i],split="/"))
    BASENAMEFILE[i] <- splittedName[length(splittedName)]
    outPdf[i] <- paste(filez[i],".pdf",sep="")
    if(is.na(outSuffix))
        outTxt[i] <- paste(filez[i],"HMMploidy",sep=".")
    if(!is.na(outSuffix))
        outTxt[i] <- paste(filez[i],outSuffix,"HMMploidy",sep=".")
    
}
##numeric conversion of inputs
wind <- as.numeric(wind)
minInd <- as.numeric(minInd)
maxPloidy <- as.numeric(maxPloidy)
eps <- as.numeric(eps)
#flags
isNumericChosenInd <- all(!is.na(chosenInd)) #check for choice of individuals

##individuals chosen for analysis (one by one at the moment)
##if chosenInd id NA it will be assigned as all individuals later
if(isNumericChosenInd)
    chosenInd <- eval( parse( text=paste("c(",chosenInd,")",sep="") ) )


##print Rscript input
cat("----------\nfileList: ", fileList, " wind: ", wind," minInd: ", minInd, " chosenInd: ", chosenInd ," maxPloidy: ", maxPloidy, " alpha: ", alpha, " beta: ", beta, " quantileTrim: ", quantileTrim, " eps: ", eps,  "\n-----------\n" )

############################
### Supporting functions ###
############################


hmmPlotting <- function(hmm, V, axLabel=TRUE, truePl=NA, main="Inferred ploidies", propStates, CNV){
    options(warn=-1)
    loci <- hmm$lociSNP
    chrName <- unique(hmm$chrSNP)
    chr <- hmm$chrSNP
    
    if(length(chrName)==1)
        borderVal <- round( seq(min(loci), max(loci), length.out=min(30,length(V$y)) ) )
    if(length(chrName)>1){
        lth <- max(2,round(30/length(chrName)))
        borderVal <- round( seq(min(loci[chr==chrName[1]]), max(loci[chr==chrName[1]]), length.out=lth ) )
        for(nn in chrName[-1])
            borderVal <- c(borderVal, round( seq(min(loci[chr==nn]), max(loci[chr==nn]), length.out=lth ) ))
    }
    
    xlabels=c()
    if(max(loci)>=1e+6){
        XLAB="Position (Mb)"
        for(i in 1:(length(borderVal)))
            xlabels[i] <- sprintf("%.1f", borderVal[i]/(1e+6))
    }
    if(max(loci)<1e+6){
        XLAB="Position (Kb)"
        for(i in 1:(length(borderVal)))
            xlabels[i] <- sprintf("%.1f", borderVal[i]/(1e+3))
    }

    
    layout(matrix(c(1,1,1,1,2,2), nrow = 3, ncol = 2, byrow = TRUE))

    oldTruePl <- truePl
    if(any(is.na(V$y)))
        V$y[is.na(V$y)] = 0
    
    if(is.na(truePl))
        truePl <- V$y

    
    plot( x=1:length(c(V$y)) ,c(V$y), pch=15, lwd=.75, col="navyblue", main=main, xaxt="n", yaxt="n", bty="n", ylab="ploidy", xlab=XLAB, ylim=c(min(V$y,truePl)-.5, max(V$y,truePl)+1 ), cex=.75, cex.main=1.4, cex.lab=1.2)

    #print(CNV)

    #print(V$y[CNV])
    #points( x=which(CNV), y = V$y[CNV] -.075 , pch=24, lwd=.75, col='coral', cex=.75)
    
    if(!is.na(oldTruePl))
        points( c(truePl)-.075 , pch=15, lwd=.75, col="coral", cex=.75)
   
    abline( h=seq( min(V$y,truePl), max(V$y,truePl) ), col="gray" )
    if( length(borderVal)<50 )
        axis( side=1, at=seq(1,length(V$y),length.out=length(borderVal)), labels=xlabels, las=2, cex=.8 )
    axis( side=2, at=seq(min(V$y), max(V$y)) )

    postProb <- hmm$postprob
    counter=1
    for(yValue in intersect(hmm$states,seq(min(V$y,truePl), max(V$y,truePl)))){
        polygon( x=c( length(V$y),1, seq(1,length(V$y)), length(V$y) ), y= yValue + 0.025 + c( 0, 0, postProb[,counter], 0 )/2, col="deepskyblue1", border=NA)
        lines( x=c(length(V$y)+.4,length(V$y)+.4), y=c(yValue+.025,yValue+0.525), col="deepskyblue1" )
        text(labels="0", x=length(V$y)+.4, y=yValue-.1)
        text(labels="1", x=length(V$y)+.4, y=yValue+.6)
        mtext(text=rep("Posterior",length(yValue)), side=4, at=yValue+.25, cex=.6)
        lines( x=rep(length(V$y)+1,2), y=c(yValue+.025,yValue+propStates[counter]/2), col="deepskyblue1", lwd=6 )
        counter=counter+1
    }
    
    legendCol = c("navyblue","deepskyblue1")
    legendPch = c(15,NA)
    legendTxt = c( "Inferred Ploidy", "Posterior Prob.")
    if(is.null(truePl)){
        legendCol = c("navyblue","deepskyblue1")
        legendPch = c(15,NA)
        legendTxt = c( "Inferred Ploidy", "Posterior Prob.")
    }
    
    legend(x=1, y = max(V$y,truePl)+1, legend=legendTxt, col = legendCol, lwd = c(NA,10,NA), pch = legendPch, bty = "n", ncol = length(legendPch), cex=1.4)

    ##PLOT DEPTH
    plot(hmm$count[,1], type="p", lwd=2, col="deepskyblue1", xlab=XLAB, main="Window-Mean Depth", xaxt="n", ylab="depth", bg=3, ylim=c(min(hmm$count[,1])-.1*min(hmm$count[,1]),max(hmm$count[,1])+.2*max(hmm$count[,1])), cex.lab=1.2, cex.main=1.3)
    if( length(borderVal)<50 )
        axis( side=1, at=seq(1,length(V$y),length.out=length(borderVal)), labels=xlabels, las=2, cex=1 )
    abline(h=hmm$mu, col="coral")
    legend(x=1, y = max(hmm$count[,1])+.25*max(hmm$count[,1]), legend=c("Mean Depth", "Neg.Bin. Mean"), col = c("deepskyblue1","coral"), lwd = rep(3,3), lty=c(NA,1), pch = c(20,NA), bty = "n", ncol = 2, cex=1.4)
options(warn=0)
}


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

# nbh_init  Initalize parameters for nbh_em
# 			Function nbm_em (NB mixture model) is used to find alpha,
# 			beta, and wght (mixprop); wght (1xN) is repeated N times row-wise
# 			to represent the initial TRANS for the subsequent nbh_em training
# 			Use: nbh0 <- nbh_init(count, K)
#			nbh0: list(TRANS, alpha, beta)

nbm_em <- function(count, alpha, beta, wght, NBM_NIT_MAX=250, NBM_TOL=1e-2){
	
	# Data length
	Total <- length(count)	
	
	count <- matrix(count, nrow=Total, 1)
	# Number of mixture components
	N <- length(alpha)
	wght <- matrix(wght, 1, N)
	alpha <- matrix(alpha, 1, N)
	beta <- matrix(beta, 1, N)
	
	
	# Save initial alpha and beta in case error occurs in the first EM
	wght0 <- wght
	alpha0 <- alpha
	beta0 <- beta


	# Compute log(count!), the second solution is usually much faster
	# except if max(count) is very large
	cm <- max(count)
	if(cm > 50000){
		dnorm <- as.matrix(lgamma(count + 1))
	} else {
		tmp <- cumsum(rbind(0, log(as.matrix(1:max(count)))))
		dnorm <- as.matrix(tmp[count+1])
	}
	
	# Variables
	logl <- matrix(0, ncol=NBM_NIT_MAX)
	postprob <- matrix(0, Total, N)
	
		
	# Main loop of the EM algorithm
	for(nit in 1:NBM_NIT_MAX){
		
		# 1: E-Step, compute density values
		postprob <- exp( matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha))
				- count %*% log(1+beta) + lgamma(count %*% matrix(1, ncol=N) + matrix(1, nrow=Total) %*% alpha)
				- dnorm %*% matrix(1, ncol=N) )
    
                # set zero value to the minimum double to avoid -inf when applying log
                # due to large dnorm (or essential large count)
		postprob <- apply(postprob, 2, function(x) {x[x==0] <- .Machine$double.xmin; x})
        
		postprob <- postprob * (matrix(1, Total, 1) %*% wght)
    																
		# Compute log-likelihood
		logl[nit] <- sum(log(apply(postprob, 1, sum)))				
		postprob <- postprob / (apply(postprob, 1, sum) %*% matrix(1, ncol=N))

		# 4: M-Step, reestimation of the mixture weights
		wght <- apply(postprob, 2, sum)		
		wght <- wght / sum(wght)
												
		# 5: CM-Step 1, reestimation of the inverse scales beta with alpha fixed
		eq_count <- apply(postprob, 2, sum)		
		mu <- (t(count) %*% postprob) / eq_count
		beta <- alpha / mu
		
		
		# 5: CM-Step 2, reestimation of the shape parameters with beta fixed
		# Use digamma and trigamma function to perfom a Newton step on
		# the part of the intermediate quantity of EM that depends on alpha
		# Compute first derivative for all components
		# Use a Newton step for updating alpha
		# Compute first derivative		
		grad <- eq_count * (log(beta / (1+beta)) - digamma(alpha)) +
					apply(postprob * digamma(count %*% matrix(1,ncol=N) +
					matrix(1,nrow=Total) %*% alpha), 2, sum)
					
		# and second derivative
		hess <- -eq_count * trigamma(alpha) +
    			apply(postprob * trigamma(count %*% matrix(1,ncol=N) +
    			matrix(1, nrow=Total) %*% alpha), 2, sum)
    
    
		# Newton step
		tmp_step <- - grad / hess
		tmp <- alpha + tmp_step
		
		# erroneous update occurs, give up and return the previous trained parameters
		if(any(is.na(tmp))) {	
			warning("Updated alpha becomes NA probably due to bad initial alpha or insuff. data")
			return(list(wght=wght0, alpha=alpha0, beta=beta0, 
							logl=logl, postprob=postprob))
			
		}
		
		# When performing the Newton step, one should check that the intermediate
		# quantity of EM indeed increases and that alpha does not become negative. In
		# practise this is almost never needed but the code below may help in some
		# cases (when using real bad initialization values for the parameters for
		# instance)
		while (any(tmp <= 0)){
			warning(sprintf("Alpha (%.4f) became negative! Try smaller (10%s) Newton step ...\n", tmp_step,"%"))
			tmp_step <- tmp_step/10
			tmp <- alpha + tmp_step
		}
		alpha <- tmp
		
		# stop iteration if improvement in logl is less than TOL (default 10^-5)
		if(nit > 1 && abs((logl[nit] - logl[nit-1])/logl[nit-1]) < NBM_TOL){
			logl <- logl[1:nit]			
			break
		}
	}

	list(wght=wght, alpha=alpha, beta=beta, 
		logl=logl, postprob=postprob)
				
}


##Estep for the conditional EM optimization
EStep <- function(count,delta,TRANS,alpha,beta,genolike){

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
    
    dens2 <- exp( dens2 )
    dens2 <- apply(dens2, 2, function(x) {x[x==0 | is.na(x) | is.nan(x)] <- .Machine$double.xmin; x})
    
    forwrd2[1,] <- delta*dens2[1,]
		
    for(t in 2:Total){
        forwrd2[t,] <- (forwrd2[t-1,] %*% TRANS) * dens2[t,]
        scale2[t] <- sum(forwrd2[t,])
        forwrd2[t,] <- forwrd2[t,] / scale2[t]
    }
				
    llk <- log(sum(forwrd2[Total,])) + sum(log(scale2))

    bckwrd2[Total,] <- matrix(1, ncol=N)
    
    for(t in (Total-1):1) {			
        bckwrd2[t,] <- (bckwrd2[t+1,] * dens2[t+1,]) %*% t(TRANS)
        bckwrd2[t,] <- bckwrd2[t,] / scale2[t]
    }

    bckwrd2[ is.nan(bckwrd2) | is.infinite(bckwrd2) | is.na(bckwrd2)] = .0001
    forwrd2[ is.nan(forwrd2) | is.infinite(forwrd2) | is.na(forwrd2)] = .0001
    
    ni <- forwrd2 * bckwrd2
    ni <- ni / ( apply(ni, 1, sum) %*% matrix(1,ncol=N) )

    return(list(forwrd=forwrd,bckwrd=bckwrd,ni=ni,llk=llk,dens=dens,forwrd2=forwrd2,bckwrd2=bckwrd2,dens2=dens2))
}

##Mstep for the conditional EM optimization
MStep <- function(E,count,TRANS,alpha,beta){
    
    remStates = FALSE
    N <- nrow(TRANS)
    Total <- dim(count)[1]
    bckwrd <- E$bckwrd; bckwrd2 <- E$bckwrd2
    forwrd <- E$forwrd; forwrd2 <- E$forwrd2
    ni <- E$ni
    dens <- E$dens; dens2 <- E$dens2
    TRANS0=TRANS
    
    delta <- ni[1,]
    TRANS <- TRANS * (t(forwrd2[1:(Total-1),]) %*% (dens2[2:Total,] * bckwrd2[2:Total,]))
    TRANS <- TRANS / (apply(TRANS, 1, sum) %*% matrix(1,ncol=N))

    for(j in 1:nrow(TRANS)){
        if(sum(is.nan(TRANS[j,]))>1){
            remStates = TRUE
            break;
        }
    }
    
    if(remStates)
        return(list(delta=delta,TRANS=TRANS0,alpha=alpha,beta=beta,remStates=remStates))
    
    eq_count <- apply(ni, 2, sum)
    beta <- alpha / ( (t(count) %*% ni) / eq_count )
    
    grad <- eq_count * (log(beta / (1+beta)) - digamma(alpha)) + apply(ni * digamma(count %*% matrix(1,ncol=N) + matrix(1,nrow=Total) %*% alpha), 2, sum)
    
    hess <- -eq_count * trigamma(alpha) + apply(ni * trigamma(count %*% matrix(1,ncol=N) + matrix(1, nrow=Total) %*% alpha), 2, sum)
    
    tmp_step <- - grad / hess
    tmp <- alpha + tmp_step
    
    if(any(is.na(tmp)))
        return(list(delta=delta,TRANS=TRANS,alpha=alpha,beta=beta,remStates=remStates))
    
    countTooMany <- 1
    while (any(tmp <= 0) & countTooMany<50){
        warning(sprintf("Alpha (%.4f)<0 ! Try smaller (10%s) Newton step ...\n", tmp_step,"%"))
        tmp_step <- tmp_step/10
        tmp <- alpha + tmp_step
        countTooMany <- countTooMany + 1
    }
    
    alpha <- tmp
    
    return(list(delta=delta,TRANS=TRANS,alpha=alpha,beta=beta,remStates=remStates))
}


##EM algorithm for HMM optimization
nbHMM <- function(count, delta, TRANS, alpha, beta, genolike=0, ws=1, PLOIDYMAX=maxPloidy, NBH_NIT_MAX=10000, NBH_TOL=1e-5, MAXALPHA=1e7, MAXBETA=1e7, keepStates=NULL){

    L <- length(delta)
    stateVec <- 1:PLOIDYMAX
    if(!is.null(keepStates))
        stateVec <- keepStates
 
    geno <- matrix(genolike[,stateVec], ncol=L)

    bicIter <- 0
    Total <- dim(count)[1]	
    if(any(count < 0)) stop("Data contains values <0")
    
    N <- nrow(TRANS)
    alpha <- matrix(alpha, 1, N)
    beta <- matrix(beta, 1, N)	

    logl <- matrix(0, ncol=NBH_NIT_MAX)
    
    alpha0 <- alpha
    beta0 <- beta
    TRANS0 <- TRANS
    postprob0 <- matrix(1, Total, N)/N
    logl0 <- logl
    dens0 <- matrix(0, Total, N)
    rate <- +Inf
    checkBIC <- TRUE
    
    ##Main loop of the EM algorithm
    for(nit in 1:NBH_NIT_MAX) {

        bicIter <- bicIter + 1

        ##check after some steps for some convergence and try to remove 1 state
        if(rate<=.0001 && bicIter>=30 && checkBIC && length(alpha)>1){

            if(length(alpha)==2){ #if removing one state leaves only one ploidy
                alpha=matrix(alpha,ncol=length(alpha))
                beta=matrix(beta,ncol=length(beta))
                cat("\t==>reduction to 1 state start\n")
                compLL <- c()
                LLtest <- c()
                for(kk in 1:length(stateVec)){#try out one state at a time
                    diff <- +Inf
                    contDiff <- 0
                    gg <- geno[,kk]
                    aa <- alpha[,kk]; aa <- matrix(aa,ncol=length(aa))
                    bb <- beta[,kk]; bb <- matrix(bb,ncol=length(bb))
                    nanFlag <- FALSE
                    
                    while(diff >= 0.0001 & contDiff < 150){
                        contDiff <- contDiff + 1
                        resSingle1 <- MStepSingle(count,aa,bb,gg)
                        llk1 <- resSingle1$llk
                        if(is.nan(llk1)){
                            nanFlag=TRUE
                            break
                            }
                        aa <- resSingle1$alpha
                        bb <- resSingle1$beta
                        contDiff <- contDiff + 1
                        resSingle2 <- MStepSingle(count,aa,bb,gg)
                        llk2 <- resSingle2$llk
                        if(is.nan(llk2)){
                            nanFlag=TRUE
                            break
                            }
                        aa <- resSingle2$alpha
                        bb <- resSingle2$beta
                        diff <- abs( (llk1 - llk2)/llk1 )
                    }

                    if(nanFlag==TRUE)
                        compLL <- c( compLL, -Inf )#models' AIC
                    if(nanFlag==FALSE)
                        compLL <- c( compLL, 2*llk2 )#models' AIC

                    
                }
                compLL[is.nan(compLL)] = min(compLL[!is.nan(compLL)]) - 1

                diff=+Inf
                LLtestOld=c(+Inf,+Inf)
                contDiff <- 0
                nanFlag=FALSE
                ##get the AIC using 2 states
                while(diff >= 0.0001 & contDiff < 150){
                    contDiff = contDiff + 1
                    E1 <- EStep(count,delta,TRANS,alpha,beta,geno)
                    if(is.nan(E1$llk)){
                        nanFlag=TRUE
                        break
                    }
                    LLtestOld[1] <- E$llk
                    M1 <- MStep(E1,count,TRANS,alpha,beta)
                    
                    E2 <- EStep(count,M1$delta,M1$TRANS,M1$alpha,M1$beta,geno)
                    if(is.nan(E2$llk)){
                        nanFlag=TRUE
                        break
                        }
                    LLtestOld[2] <- E2$llk
                    M2 <- MStep(E2,count,M1$TRANS,M1$alpha,M1$beta)
                    TRANS <- M2$TRANS; alpha <- M2$alpha; beta <- M2$beta
                    diff <- abs( (LLtestOld[1] - LLtestOld[2])/LLtestOld[1] )
                    if(is.nan(diff))
                        diff=0
                }

                K <- length(alpha)
                if(nanFlag==TRUE)
                    oldBIC <- -Inf
                if(nanFlag==FALSE)
                    oldBIC <- 2*E2$llk - 2*((K+1)*K - 1)
                newBIC <- max(compLL)
                whichBIC <- which.max( compLL )
                if(is.nan(newBIC))
                    newBIC <- -Inf

                #cat("\tBest new BIC ", newBIC," old BIC", oldBIC, "reduce states ", (newBIC>oldBIC), "\n") 
                cat("\t   reduce states:", (newBIC>oldBIC), "\t",sep="")
                if(oldBIC > newBIC){
                    checkBIC <- FALSE
                    cat("States: ", stateVec,"\n")
                    break
            }
            else{
                geno <- geno[,whichBIC ]
                alpha <- alpha[,whichBIC]
                beta <-  beta[,whichBIC]
                stateVec <- stateVec[ whichBIC ]
                bckwrd <- matrix(1,nrow(count),1)
                TRANS <- matrix(TRANS[whichBIC, whichBIC], 1, 1)            
            }
                cat("States: ", stateVec,"\n")
                #cat("\talpha: ", alpha, " beta: ", beta," mu: ", alpha/beta, "\n")
                break
        }

            alpha <- matrix(alpha,ncol=length(alpha))
            beta <- matrix(beta,ncol=length(beta))
            K <- dim(alpha)[2]
            ##if removing one state leaves at least other two, then what follows will happen
            #cat("\t==>reduction to ",K-1,"states start\n\t\talpha: ", as.vector(alpha),"\n\t\tbeta: ",as.vector(beta),"\n",sep=" ")
            cat("\t==>reduction to ",K-1,"states start\n")
            bicIter <- 0
            viterbiIter <- 0
            combMat <- combs( 1:K, K-1 )
            #sortIdx <- apply(combMat,1,is.sorted)
            #combIdx <- which(sortIdx)#
            combIdx <- 1:(dim(combMat)[1])
            compLL <- c()
            Lcomb = length(combIdx)
            TRANS2 <- TRANS; delta2 <- delta; alpha2 <- alpha; beta2 <- beta; geno2 <- geno;
            ##try each combination of states and genotype likelihoods on ploidies
            for( kk1 in combIdx ){
                for( kk2 in combIdx ){
                    geno <- geno2[ ,combMat[kk1,] ]
                    delta <- delta2[ combMat[kk2,] ]
                    alpha <- matrix(alpha2[ ,combMat[kk2,] ], nrow=1)
                    beta <- matrix(beta2[ ,combMat[kk2,] ], nrow=1)
                    diff=+Inf
                    LLtest=c(+Inf,+Inf)
                    #cat("states: ",stateVec[ combMat[kk1,] ]," -- ","alpha: ",alpha," beta: ",beta,"\n",sep="")
                    TRANS <- TRANS2[ combMat[kk2,], ]
                    TRANS <- TRANS[ ,combMat[kk2,] ]
                    contDiff = 0
                    nanFlag=FALSE
                    #optimize when one state is removed and get the AIC.
                    #nan check useful for messy real data.
                    while(diff >= 0.001 & contDiff < 250){
                        contDiff = contDiff + 1

                        E1 <- EStep(count,delta,TRANS,alpha,beta,geno)
                        if(is.nan(E1$llk)){
                            nanFlag=TRUE
                            break
                            }
                        LLtest[1] <- E$llk
                        M1 <- MStep(E1,count,TRANS,alpha,beta)
                        
                        E2 <- EStep(count,M1$delta,M1$TRANS,M1$alpha,M1$beta,geno)
                        if(is.nan(E2$llk)){
                            nanFlag=TRUE
                            break
                            }
                        LLtest[2] <- E2$llk
                        M2 <- MStep(E2,count,M1$TRANS,M1$alpha,M1$beta)
                        TRANS <- M2$TRANS; alpha <- M2$alpha; beta <- M2$beta
                        
                        diff <- abs( (LLtest[1] - LLtest[2])/LLtest[1] )
                        if(is.nan(diff))
                            diff=0
                        contDiff = contDiff+1
                    }
                    if(nanFlag==TRUE)
                        compLL <- c( compLL, -Inf )
                    if(nanFlag==FALSE)
                        compLL <- c( compLL, 2*E2$llk - ((K-1)*K - 1)*2 )
                    
                    geno <- geno2
                    delta <- delta2
                    alpha <- alpha2
                    TRANS <- TRANS2

                }
            }
            compLL[is.nan(compLL)] = min(compLL[!is.nan(compLL)]) - 1
            geno <- geno2
            delta <- delta2
            alpha <- alpha2
            beta <- beta2
            diff=+Inf
            LLtestOld=c(+Inf,+Inf)
            contDiff <- 0
            
            while(diff >= 0.0001 & contDiff < 250){
                
                contDiff = contDiff + 1
                E1 <- EStep(count,delta,TRANS,alpha,beta,geno)
                if(is.nan(E1$llk))
                    break
                LLtestOld[1] <- E$llk
                M1 <- MStep(E1,count,TRANS,alpha,beta)
                        
                E2 <- EStep(count,M1$delta,M1$TRANS,M1$alpha,M1$beta,geno)
                if(is.nan(E2$llk))
                    break
                LLtestOld[2] <- E2$llk
                M2 <- MStep(E2,count,M1$TRANS,M1$alpha,M1$beta)
                TRANS <- M2$TRANS; alpha <- M2$alpha; beta <- M2$beta
                
                diff <- abs( (LLtest[1] - LLtest[2])/LLtest[1] )
                if(is.nan(diff))
                    diff=0
            }
            
            oldBIC <- 2*E2$llk - 2*((K+1)*K - 1)
            newBIC <- max(compLL)
            whichBIC <- which.max( compLL )
            if(is.nan(newBIC))
                newBIC <- -Inf
            if(is.nan(oldBIC))
                oldBIC <- -Inf
            #cat("\t\tBest new BIC ", newBIC," old BIC", oldBIC, "reduce states:", (newBIC>oldBIC), "\n") 
            cat("\t   reduce states:", (newBIC>oldBIC), "\t",sep="") 
            if(oldBIC > newBIC){
                checkBIC <- FALSE
                delta <- delta2
                TRANS <- TRANS2
                alpha <- alpha2
                beta <- beta2
                geno <- geno2

            }
            else{
                comb1 <- whichBIC %/% Lcomb + 1
                comb2 <- whichBIC %% Lcomb
                if( (whichBIC %% Lcomb) == 0 )
                    comb1 = (comb1-1)
                if( comb2==0 )
                    comb2 = Lcomb
                geno <- geno2[ ,combMat[combIdx[comb1],] ]
                if( (K-1) > 1){
                    delta <- delta2[ combMat[ combIdx[comb2] ,] ]
                    TRANS <- TRANS2[ combMat[combIdx[comb2],], ]
                    TRANS <- TRANS[ ,combMat[combIdx[comb2],] ]
                    }
                alpha <- alpha2[ ,combMat[combIdx[comb2],] ] 
                beta <-   beta2[ ,combMat[combIdx[comb2],] ] 
                stateVec <- stateVec[ combMat[combIdx[comb1],] ]
                sortIdx <- sort(alpha/beta, index.return=TRUE)$ix
                alpha=alpha[sortIdx]
                beta=beta[sortIdx]
            }
            cat("States Relation : ", stateVec,"\n")
            #cat("\t==>alpha: ", alpha, "\n\t\tbeta: ", beta, "\n\t\tmu: ", alpha/beta, "\n")
            #alpha <- matrix(alpha,ncol=length(alpha))
           
        }

        N <- nrow(TRANS)        
        E <- EStep(count,delta,TRANS,alpha,beta,geno)

        logl[nit] <- E$llk
        
        #message( sprintf( 'Iter %d:\tLLK=%.3f\tploidies=%d', (nit-1), logl[nit], dim(alpha)[2] ) )	
        if(is.nan(logl[nit])) {			
            warning("NaN logl data detected. Returning the previous training results")
            logl[nit] <- 0; TRANS=TRANS0; alpha=alpha0; beta=beta0; logl=logl0; bckwrd=postprob0; dens=dens0
            #break
        }
        
        M <- MStep(E,count,TRANS,alpha,beta)
        alphaPrev <- alpha #backup
        betaPrev <- beta
        TRANSprev <- TRANS
        deltaPrev <- delta
        TRANS <- M$TRANS #new values
        alpha <- M$alpha
        beta <- M$beta
        bckwrd <- E$ni
        dens <- E$dens2
        delta <- M$delta

        if(any(is.na(M$alpha))) {
            warning(sprintf("Updated alpha becomes NA probably %s","due to bad initial alpha or insuff. data"))
            return(list(count=count, delta=deltaPrev, TRANS=TRANSprev, alpha=alphaPrev, beta=betaPrev, logl=logl, postprob=bckwrd, dens=dens, mu = alphaPrev/betaPrev, sigma=alphaPrev/betaPrev+alphaPrev/betaPrev^2, states=stateVec)	)
        }
        if(any(is.na(M$beta))) {
            warning(sprintf("Updated beta becomes NA probably %s","due to bad initial alpha orE$dens2 insuff. data"))
            return(list(count=count, delta=deltaPrev, TRANS=TRANSprev, alpha=alphaPrev, beta=betaPrev, logl=logl, postprob=bckwrd, dens=dens, mu = alphaPrev/betaPrev, sigma=alphaPrev/betaPrev+alphaPrev/betaPrev^2, states=stateVec)	)
        }
        
        if(any(M$alpha > MAXALPHA) || any(M$beta > MAXBETA)) {
            warning(sprintf("Updated alpha (%f) or beta (%f) becomes too large probably %s",alpha, beta, "due to bad initial alpha or large count"))
            TRANS=TRANS0; alpha=alpha0; beta=beta0; logl=logl0; bckwrd=postprob0; dens=dens0
            break
        }
	############ Check END ############
        if(nit > 1)
            if(logl[nit] < logl[nit-1])
                bicIter <- 30
      
	rate <- abs((logl[nit] - logl[nit-1])/logl[nit-1]) 
        if(nit > 1 && rate < NBH_TOL && !checkBIC){
            logl <- logl[1:nit]
            break
        }		

        TRANS0 <- TRANS
        alpha0 <- alpha
        beta0 <- beta
        postprob0 <- bckwrd
        dens0 <- dens
        logl0 <- logl
    }

                                   
    list(count=count, delta=delta, TRANS=TRANS, alpha=alpha, beta=beta, logl=logl[logl<0], postprob=bckwrd, dens=dens, mu = alpha/beta, sigma=alpha/beta+alpha/beta^2, states=stateVec, geno=geno)				
}


MStepSingle <- function(count,alpha,beta,geno){

    N <- 1
    Total <- dim(count)[1]
    bckwrd <- rep(1,Total)
    forwrd <- rep(1,Total)
    ni <- rep(1,Total)

    dens <- matrix(0, Total, N)
    scale <- rbind(1, matrix(0, nrow=Total-1))
    
    cm <- max(count)
    if(cm > 50000){
        dnorm <- as.matrix(lgamma(count + 1))
    } else {
        tmp <- cumsum(rbind(0, log(as.matrix(1:max(count)))))
        dnorm <- as.matrix(tmp[count+1])
    }

    densLog <-  matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha)) - count %*% log(1+beta) + lgamma(count + matrix(1,nrow=Total) %*% alpha) - dnorm  + geno

    dens <- exp( densLog )

    dens <- apply(dens, 2, function(x) {x[x==0 | is.na(x) | is.nan(x)] <- .Machine$double.xmin; x})        	

    delta <- 1
    TRANS <- 1
    								
    eq_count <- sum(ni)
    beta <- alpha / ( sum(count) / eq_count)
		
    grad <- eq_count * (log(beta / (1+beta)) - digamma(alpha)) + apply(ni * digamma(count %*% matrix(1,ncol=N) + matrix(1,nrow=Total) %*% alpha), 2, sum)
				
    hess <- -eq_count * trigamma(alpha) + apply(ni * trigamma(count %*% matrix(1,ncol=N) + matrix(1, nrow=Total) %*% alpha), 2, sum)

    tmp_step <- - grad / hess
    tmp <- alpha + tmp_step

    if(any(is.na(tmp)))
        return(list(delta=delta,TRANS=TRANS,alpha=alpha,beta=beta,llk=NaN))
    
    while (any(tmp <= 0)){
        warning(sprintf("Alpha (%.4f)<0 ! Try smaller (10%s) Newton step ...\n", tmp_step,"%"))
        tmp_step <- tmp_step/10
        tmp <- alpha + tmp_step
    }
		
    alpha <- tmp
    return(list(delta=1,TRANS=1,alpha=alpha,beta=beta,llk=sum(densLog)))
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
sumGeno <- function(dp,ws,loci,lociSNP=loci,findSNP=1:length(loci),avg=FALSE){  #useless?
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
        v <- dp[idx]
        if(!avg)
            return( median( v ) )
        if(avg)
            return( mean(v) )
    })
    return(unlist(res))
}

sumGenoAll <- function(dp,ws=1,loci,lociSNP=loci,avg=FALSE){  #useless?
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

#########WINDOWIZED FUNCTIONS
sumGeno <- function(g,wind,lociSNP){ #at least a SNP is present in each window with our windowBuilder function
    L <- length(g)    
    res <- sapply(1:dim(wind)[1], function(ll){
        idx <- which(lociSNP>=wind[ll,1] & lociSNP<=wind[ll,2])
        v <- g[idx]
        return( sum(v, na.rm=TRUE) )
    })
    return(res)
}

meanGeno <- function(dp,wind,loci){ #NO SNPS USED ON DEPTH AVERAGE
    L <- length(dp)    
    res <- sapply(1:dim(wind)[1], function(ll){
        idx <- which(loci>=wind[ll,1] & loci<=wind[ll,2])
        v <- dp[idx]
        return( mean(v, na.rm=TRUE) )
    })
    return(res)
}

windowsBuilder <- function(minSize, loci, lociSNP){ #do windows and merge them when one is empty
    if(length(loci)>2)
        S <- c( seq(loci[1], loci[length(loci)-1], minSize), loci[length(loci)] )
    if(length(loci)==2)
        S <- c( loci[1], loci[2] )
    if( max(loci) - min(loci) < minSize)
        S <- c( min(loci), max(loci) )
    res <- sapply(1:(length(S)-1), function(ll){
        idx <- which(lociSNP>S[ll] & lociSNP<=S[ll+1])
        if(length(idx)==0 && ll<(length(S)-1))
            return(ll+1)
        if(length(idx)==0 && ll==(length(S)-1))
            return(ll)

    })
    
    if(!is.null(unlist(res)))
        S <- S[-c(unique(unlist(res)))]

    windTable <- matrix(nrow=length(S)-1,ncol=2)
    windTable[,1] <- S[-length(S)]
    windTable[,2] <- S[-1]
    return(windTable)
}
#########



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
pGenoData <- function(f,winL,gl,nInd=1,findSNP=1:sum(winL),h=0){   #useless?
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
    }
    return( X )    
}

cppFunction('NumericVector pGenoDataAll(NumericVector f, NumericMatrix gl){
    int y = gl.ncol()-1;
    int Lf = f.size();
    
    NumericVector res(Lf);

    for(int l=0; l<Lf; l++){

        NumericVector p = dbinom( seq(0,y), y, f[l], true );
        NumericVector gg = gl(l,_);
        LogicalVector infIdx = is_infinite(p);
        
        for(int k=0; k<=y; k++)
            if(infIdx[k])
               p[k] = -1000;

        NumericVector X = gg + p;

        //logarithmic sum
        double m = max(X);
        X.erase(which_max(X));

        NumericVector diffVec = X - m;     
        diffVec = clamp(-700, diffVec, -0.0001);     
        res(l) = m + log( 1 + sum( exp( diffVec ) ) );      
    }

    return(res);
}'
)


freqsSingle <- function(major,minor,ws,loci,lociSNP=loci,findSNP=1:length(loci)){#useless?   
    L <- length(major)    
    S <- seq(loci[1],loci[length(loci)],ws)
    S <- c(S, loci[length(loci)] )
    winLength <- rep(0, length(S)-1)
    
    res <- sapply(1:(length(S)-1), function(ll){
        if(ll<(length(S)-1)){
            idx <- which(lociSNP>=S[ll] & lociSNP<S[ll+1])
            outL <- sum(loci>=S[ll] & loci<S[ll+1])
            }
        if(ll==(length(S)-1)){
            idx <- which(lociSNP>=S[ll] & lociSNP<=S[ll+1])
            outL <- sum(loci>=S[ll] & loci<=S[ll+1])
            }
        if(length(idx)==0)
            return(c())
        num <- minor[findSNP[idx]]
        den <- num + major[findSNP[idx]]
        return( c(mean(num/den, na.rm=TRUE), outL) )
    })
    return( list(winF=res[1,], winL=res[2,]) )
}


pGenoDataSingle <- function(f,gl,h=0){  #use one individual at a time #useless?
    y = ncol(gl)-
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


##Viterbi algorithm
Viterbi <- function(hmm){
    n <- dim(hmm$count)[1]
    m <- nrow(hmm$TRANS)
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    nu[1, ] <- log(hmm$delta) + log(hmm$dens[1,])
    logPi <- log(hmm$TRANS)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi, 2, max) + log(hmm$dens[i,])
    }
    if (any(nu[n, ] == -Inf)) 
        nu[is.infinite(nu)] = -700
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])
    return(list(y=y,nu=nu))
}

wind <- as.numeric(wind) #window size
minInd <- as.numeric(minInd)
maxPloidy <- as.numeric(maxPloidy)
#truePl <- eval(parse(text=truePl))

#################################################
### Begin file-by-file analysis #################
#################################################

for(i in 1:length(fileVector)){ #loop over input files
    cat("==> Analyze ", filez[i], "\n",sep="")
    ##read in the data from .mafs and .genolikes files
    GL <- fread(input=fileVector[i],sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)
    rowsGL <- dim(GL)[1]
    chr <- GL[,1]
    chrName = unique(chr) #keep unchanged
    nInd <- length( unique( GL[,3] ) )
    sites <- c() #same loci in different chromosome would get cancelled. do one chromosome at a time
    for(chrN in chrName)
        sites <- c(sites, unique( GL[chr==chrN,2] ) )

    #for(chrN in chrName)
    #    cat("chr ", chrN, "length ",sum(chr==chrN),"\n")
    
    DP <- GL[ ,5]

    #calculate allele frequencies
    majorReads <- GL[,8]
    minorReads <- GL[,9]
    freqs <- alleleFrequencies(majorReads,minorReads,nInd,minInd,eps) #frequencies (used for SNPs on a single individual)
    
    GL <- GL[ ,-c(1:9)]
    
    if(!isNumericChosenInd)
        chosenInd <- 1:nInd

    inputNames <- c()
    if(!is.na(nameList))
        inputNames <- unlist( read.table(nameList, header=FALSE, as.is=T)  )
    if(is.na(nameList))
        inputNames <- 1:nInd

    
    ###############################
    ## begin of FOR loop to      ##
    ## read one genome at a time ##
    ###############################

    fileCounter = 1
    ##open pdf plot connection
    pdf(outPdf[i])

    for(whichInd in chosenInd){ #loop over individuals
        
        chrNameVar = chrName #vector names of chromosomes that can change in the for loop
        chrVar = chr
        
        ##select single individual depth and genolikes
        idxSingle <- seq(whichInd,rowsGL,nInd)
        chrSingle <- chr[idxSingle]
        majorSingle <- majorReads[idxSingle]
        minorSingle <- minorReads[idxSingle]
        DPsingle <- DP[idxSingle]; GLsingle <- GL[idxSingle, ] #individual depth/genolikes
        ##trim the depth at the chosen quantile
        quantiles <- eval( parse( text=paste("c(",quantileTrim,")",sep="") ) )
        q <- quantile( DPsingle, quantiles )  
        idx <- which( DPsingle<=as.numeric(q[2]) & DPsingle>=as.numeric(q[1]) )
        DPsingle <- DPsingle[idx] #individual filtered data
        GLsingle <- GLsingle[idx, ] #...""
        chrSingle <- chrSingle[idx]
        sitesIndiv <- sites[idx] #......""
        freqsIndiv <- freqs[idx] #......""
        majorSingle <- majorReads[idx] #......""
        minorSingle <- minorReads[idx] #......""
        idxTot = as.vector( sapply(idx, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        GLfiltered <- GL[idxTot, ] #all data filtered
        DPfiltered <- DP[idxTot] #......""
        chrVar <- chrVar[idxTot]
        
        ##find SNPs with thresholds .1<f<.9 (or different from input) in the individual
        SNPfilter <- eval( parse( text=paste("c(",SNPtrim,")",sep="") ) )
        findSNP <- which( freqsIndiv>SNPfilter[1] & freqsIndiv<SNPfilter[2] )
        freqsSNP <- freqsIndiv[findSNP]
        sitesSNP <- sitesIndiv[findSNP]
        chrSNP <- chrSingle[findSNP]
        totSNP <- as.vector( sapply(findSNP, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        chrVar <- chrVar[totSNP]
        chrNameVar = unique(chrVar) #remember to update chromosome names because you might remove a chromosome
        
        ##remove contigs without loci or make smaller window if not enough loci
        rmChrSNP <- c()
        rmChrSingle <- c()
        rmChrTot <- c()
        rmChrName <- c()
        if(length(chrNameVar)>1){
            for(nn in chrNameVar){
                ctgSites <- sitesSNP[chrSNP==nn]
                #print(nn)
                if(length(ctgSites)==0)
                    ctgSites=1 #to avoid warnings           
                if( max(ctgSites) - min(ctgSites) == 0 ){
                    #cat(max(ctgSites) - min(ctgSites), " in ", nn,"\n")
                    rmChrSNP <- c(rmChrSNP, which(chrSNP==nn))
                    rmChrSingle <- c(rmChrSingle, which(chrSingle==nn))
                    rmChrTot <- c(rmChrTot, which(chrVar==nn))
                    rmChrName <- c(rmChrName, nn)
                }
            }
            cat("Removed ", length(rmChrName), " contigs/chromosomes that had no SNPs\n", sep="")
            if(!is.null(rmChrSNP)){
                findSNP <- findSNP[ -c(rmChrSNP)  ]
                freqsSNP <- freqsSNP[ -c(rmChrSNP) ]
                sitesSNP <- sitesSNP[ -c(rmChrSNP) ]
            }
            if(!is.null(rmChrSingle)){
                sitesIndiv <- sitesIndiv[ -c( rmChrSingle ) ]
                DPsingle <- DPsingle[-c( rmChrSingle )] #individual filtered data
                ##GLsingle <- GLsingle[-c( rmChrSingle ), ] #...""
                freqsIndiv <- freqs[-c( rmChrSingle )] #......""
                majorSingle <- majorReads[-c( rmChrSingle )] #......""
                minorSingle <- minorReads[-c( rmChrSingle )] #......""      
                chrSingle <- chrSingle[ -c( rmChrSingle ) ]
            }
            if(!is.null(rmChrSNP))
                chrSNP <- chrSNP[ -c( rmChrSNP ) ]
            if(!is.null(rmChrName))
                chrNameVar <- chrNameVar[ !(chrNameVar %in% rmChrName) ] 
        }

        
        if(length(chrNameVar)==1)
            windTable <- windowsBuilder(wind, sitesIndiv, sitesSNP)
        if(length(chrNameVar)>1){
            chrTable <- c()
            rowIndex <- 1
            ##windTable <- matrix(0,nrow=length(sitesSNP),ncol=2)
            subSitesIndiv <- sitesIndiv[chrSingle==chrNameVar[1]]
            subSitesSNP <- sitesSNP[chrSNP==chrNameVar[1]]
            windTable <- windowsBuilder(wind, subSitesIndiv, subSitesSNP)
            chrTable <- c(chrTable, rep(chrNameVar[1], dim(windTable)[1]))
            for(chrIdx in 2:length(chrNameVar)){
                subSitesIndiv <- sitesIndiv[chrSingle==chrNameVar[chrIdx]]
                subSitesSNP <- sitesSNP[chrSNP==chrNameVar[chrIdx]]
                windTable <- rbind(windTable, windowsBuilder(wind, subSitesIndiv, subSitesSNP))
                chrTable <- c(chrTable, rep(chrNameVar[chrIdx], dim(windTable)[1]-length(chrTable)))
            }
        }
        
        geno2 <- matrix(0, nrow=maxPloidy, ncol=length(freqsSNP))
        ##if(strcmp(useGeno,"yes"))
        for(pp in 1:maxPloidy)
            for(nn in chrNameVar)
                geno2[pp,chrSNP==nn] = pGenoDataAll( f=freqsSNP[chrSNP==nn], gl=as.matrix(readGL( findSNP[chrSNP==nn], pp, nInd=1, GLsingle )) )
            
        ##...and per window
        if(length(chrNameVar)==1)
            geno <- apply( geno2, 1, function(x) sumGeno(x,windTable,sitesSNP) )
        if(length(chrNameVar)>1){
            geno <- apply( geno2, 1, function(x) sumGeno(x,matrix(windTable[chrTable==chrNameVar[1],],ncol=2),sitesSNP[chrSNP==chrNameVar[1]]) )
            for(chrIdx in 2:length(chrNameVar))
                geno <- rbind(geno, apply( geno2, 1, function(x) sumGeno(x,matrix(windTable[chrTable==chrNameVar[chrIdx],], ncol=2),sitesSNP[chrSNP==chrNameVar[chrIdx]]) ))
        }
        
        if(length(chrNameVar)==1)
            DPmean <- meanGeno( DPsingle, windTable, sitesIndiv)
        if(length(chrNameVar)>1){
            DPmean <- meanGeno( DPsingle[chrSingle==chrNameVar[1]], matrix(windTable[chrTable==chrNameVar[1],],ncol=2), sitesIndiv[chrSingle==chrNameVar[1]])
            for(chrIdx in 2:length(chrNameVar))
                DPmean <- c( DPmean, meanGeno( dp=DPsingle[chrSingle==chrNameVar[chrIdx]], wind=matrix(windTable[chrTable==chrNameVar[chrIdx],],ncol=2), sitesIndiv[chrSingle==chrNameVar[chrIdx]]) )
        }
        
        ##clean from NA, NaN or infinite values
        keepSites <- apply( geno, 1, function(x) sum(is.na(x) | is.nan(x) | is.infinite(x))==0 )
        DPmean <- DPmean[which(keepSites)]
        geno <- geno[which(keepSites), ]

        ##rescale likelihood of the data (avoids underflow)
        genoResc = geno
        genoResc <- t( apply( geno , 1, logRescale ) )
        genoResc[genoResc>-.00001]=-.00001
        
        ##some initial parameters
        delta=rep(1/maxPloidy,maxPloidy) #i think it is ok without prior info
        Pi0=matrix(1/maxPloidy,nrow=maxPloidy,ncol=maxPloidy)
        count <- matrix(DPmean,ncol=1)

        cat(sprintf("    Initialize parameters for %d ploidy numbers\n", maxPloidy))
        alpha <- tail(quantile(count[count[,1]>0,1], probs=seq(0, 1, 1/maxPloidy), na.rm=TRUE), maxPloidy)
        beta <- rep(1, maxPloidy)	
        delta <- rep(0.5, maxPloidy)
        nbm <- nbm_em(count, alpha, beta, delta)	
        alpha <- nbm$alpha	
        beta <- nbm$beta	
        Pi0 <- matrix(rep(nbm$wght, maxPloidy), ncol=maxPloidy, byrow=TRUE)
        ## Order the parameters s.t. they are in increasing order of mean values
        myorder <- order(alpha/beta)	
        alpha <- matrix(alpha[myorder],ncol=length(myorder))
        beta <- matrix(beta[myorder],ncol=length(myorder))
        Pi0 <- Pi0[myorder, myorder]		
        
        cat("    N.samples ",nInd,"\n    alpha0: ",alpha,"\n    beta0: ",beta,"\n",sep=" ")        
        ##start the HMM
        genoMatrix = matrix(0, nrow(genoResc), ncol(genoResc))
        if(strcmp(useGeno,"yes"))
            genoMatrix = genoResc

        hmmRes <- nbHMM(count, alpha=alpha, beta=beta, TRANS=Pi0, delta=delta, genolike=genoMatrix, PLOIDYMAX=maxPloidy)
        hmmRes$'lociSNP' = sitesSNP;
        ##Clean the matrix
        cat("    Cleaning the transition matrix\n")
        mat = hmmRes$TRANS
        keepIdx = c()
        rmv = 1:dim(mat)[1]
        p = .05
        mat2 <- mat
        print(mat)
        if(length(mat)[1] > 1 & any(diag(mat)<p)){
            for(l in 1:dim(mat)[1])
                if(mat[l,l]>p)
                    keepIdx <- c(keepIdx,l)
            mat2 <- matrix(mat[keepIdx,keepIdx], nrow=length(keepIdx))
            for(l in 1:length(keepIdx))
                mat2[l,] = mat2[l,]/sum(mat2[l,])
            ##print(mat2)
            ## Cleaned data restarted estimation

            L <- nrow(mat2)
            delta <- rep(1/L, L) #i think it is ok without prior info
            Pi0 <- mat2
            keptStates <- hmmRes$states[keepIdx]
            
            cat(sprintf("    Reinitialize parameters for %d ploidy numbers\n", length(keptStates) ))
            alpha <- tail(quantile(count[count[,1]>0,1], probs=seq(0, 1, 1/L), na.rm=TRUE), L)
            beta <- rep(1, L)	
            delta <- rep(0.5, L)
            nbm <- nbm_em(count, alpha, beta, delta)	
            alpha <- nbm$alpha	
            beta <- nbm$beta	
            Pi0 <- matrix(rep(nbm$wght, L), ncol=L, byrow=TRUE)
	    ## Order the parameters s.t. they are in increasing order of mean values
            myorder <- order(alpha/beta)	
            alpha <- matrix(alpha[myorder],ncol=length(myorder))
            beta <- matrix(beta[myorder],ncol=length(myorder))
            Pi0 <- matrix(Pi0[myorder, myorder], L)
            
            cat("    N.samples ",nInd,"\n    alpha0: ",alpha,"\n    beta0: ",beta,"\n",sep=" ")        
            ##start the HMM
            hmmRes <- nbHMM(count, alpha=alpha, beta=beta, TRANS=Pi0, delta=delta, genolike=genoMatrix, PLOIDYMAX=maxPloidy, keepStates=keptStates)
            hmmRes$'lociSNP' = sitesSNP;
        } ######end matrix cleaning if necessary

        if(!strcmp(useGeno,"yes")){
            ##find the CNVs according to the negative binomial HMM states
            if( length(hmmRes$states)>1 )
                hmmOut <- Viterbi(hmmRes)$y
            if( length(hmmRes$states)==1 )
                hmmOut <- rep( 1, length(hmmRes$dens[,1]) )
            ##print(hmmOut)
            muHmm = c()
            seenState = hmmOut[1] #unique does order things from command line Rscript. I do unique by hand.
            for(ll in hmmOut)
                if(!any(seenState == ll))
                    seenState <- c(seenState,ll)
            for(ll in seenState)
                muHmm <- c(muHmm, mean(count[hmmOut==ll,1]))
            
            ##muHmm <- sort(hmmRes$mu, index.return=TRUE)$x
            ##muIdx <- sort(hmmRes$mu, index.return=TRUE)$ix
            muIdx <- sort(muHmm, index.return=TRUE)$ix
            muHmm <- sort(muHmm, index.return=TRUE)$x
            ##print(muHmm)
            ##print(unique(hmmOut))
            ##print(genoResc)
            postProb1 = hmmRes$postprob
            postProb2 = postProb1
            dens1 = hmmRes$dens
            dens2 = dens1
            geno1 = hmmRes$geno
            geno2 = geno1
            
            nbStates <- unique(hmmOut)[muIdx]
            postProbRm <- c()
            if(length(nbStates)==1){
                hmmOut[1:length(hmmOut)] = which.max( apply(genoResc, 2, median) )#which.max( colSums( genoResc ) )
                uniquePloidy = which.max( apply(genoResc, 2, median) )
                isCNV=rep(FALSE,length(hmmOut))
            }
            if(length(nbStates)>1){         
                nbPloidy=rep(0, length(nbStates))
                truePloidy = c()
                for(nn in 1:length(nbStates)){
                    ##print(hmmOut==nbStates[nn])
                    genoBuffer = as.matrix(genoResc[hmmOut==nbStates[nn],])
                    ##print(apply(genoBuffer, 2, median))
                    ##print(colSums( genoBuffer ))
                    ##truePloidy[nn] = which.max( colSums( genoBuffer ) )
                    truePloidy[nn] = which.max( apply(genoBuffer, 2, median) )
                    ##postProb2[hmmOut==nbStates[nn],nn] = postProb1[hmmOut==nbStates[nn],truePloidy[nn]]
                }
                cat("trueploidy",truePloidy,"\n",sep=" ")
                isCNV=rep(FALSE,length(hmmOut))
                uniquePloidy = unique(truePloidy)
                ##print(uniquePloidy)
                if(length(uniquePloidy) == length(truePloidy)){
                    for(uu in 1:length(uniquePloidy)){
                        st = which(truePloidy == uniquePloidy[uu])
                        hmmOut[hmmOut==nbStates[st[1]]] = truePloidy[st[1]]
                    }
                }
                if(length(uniquePloidy) < length(truePloidy)){
                    for(uu in 1:length(uniquePloidy)){
                        st = which(truePloidy == uniquePloidy[uu])
                        ##print(length(st))
                        if(length(st)>1){
                            for(s in st[-1]){
                                ##print(hmmOut==nbStates[s])
                                ##isCNV[hmmOut==nbStates[s]] = TRUE
                                hmmOut[hmmOut==nbStates[s]] = truePloidy[st[1]]
                                postProb2[hmmOut==nbStates[s],st[1]] = postProb1[hmmOut==nbStates[nn],s]
                                dens2[hmmOut==nbStates[s],st[1]] = dens1[hmmOut==nbStates[nn],s]
                                geno2[hmmOut==nbStates[s],st[1]] = geno1[hmmOut==nbStates[nn],s]
                                postProbRm = c(postProbRm, s)
                            }
                        }
                    }
                }
            }
            print("CNV ploidy")
            ##print(hmmOut)
            alpha2 = hmmRes$alpha
            beta2 = hmmRes$beta
            delta2 = hmmRes$delta
            TRANS2 = hmmRes$TRANS
            mu2 = hmmRes$mu
            sigma2 = hmmRes$sigma
                       
            hmmRes$states = uniquePloidy
            if(!is.null(postProbRm)){
                postProb2 = as.matrix(postProb2[,-postProbRm])
                ##print(postProb2)
                for(nr in 1:dim(postProb2)[1]){
                    x = postProb2[nr,]
                    if(sum(x)>0)
                        postProb2[nr,] = x/sum(x)
                    if(sum(x)==0)
                        postProb2[nr,] = 1/rep(length(x),length(x))
                } 
                ##print(postProb2)
                hmmRes$postprob = postProb2
                hmmRes$alpha = alpha2[-postProbRm]
                hmmRes$beta = beta2[-postProbRm]
                delta2 = delta2[-postProbRm]
                hmmRes$TRANS = as.matrix(TRANS2[-postProbRm,-postProbRm])
                hmmRes$dens = as.matrix(dens2[,-postProbRm])
                hmmRes$geno = as.matrix(geno2[,-postProbRm])
                ##hmmRes$mu = mu2[-postProbRm]
                ##hmmRes$sigma = sigma2[-postProbRm]
            }
            if(sum(delta2)==0) delta2 = rep(1,length(delta2))/length(delta2)
            hmmRes$delta = delta2
            for(nr in 1:dim(TRANS2)[1]){
                x = TRANS2[nr,]
                if(sum(x)>0)
                    TRANS2[nr,] = x/sum(x)
                if(sum(x)==0)
                    TRANS2[nr,] = 1/rep(length(x),length(x))
            } 
        }
        ##print(str(hmmRes))
        ##print(length(hmmRes$states))
        ##if(strcmp(useGeno,"yes")){
        ##apply Viterbi and properly order ploidy labels
        if( length(hmmRes$states)==1 )
            V <- list(y=rep( hmmRes$states, length(hmmRes$dens[,1]) ), nu=matrix(1,length(hmmRes$dens[,1]),1) )
        if( length(hmmRes$states)>1 ){
            V <- Viterbi(hmmRes)
            for( ll in seq(length(sort(hmmRes$states, decreasing=TRUE)), 1, -1) ){
                idx <- (V$y==ll)
                V$y[idx] <- hmmRes$states[ll]
            }
        }
        hmmOut <- V$y
        ##}
   
        ## TO DO
        ##return something in a file
        cat("File: ",BASENAMEFILE[i],". Sample: ",inputNames[whichInd],"\n",sep="",file=outTxt[i],append=!(fileCounter==1)) #info
        
        cat(hmmRes$delta,"\n",file=outTxt[i],sep="\t",append=TRUE) #starting probabilities
        cat(hmmRes$TRANS,"\n",file=outTxt[i],sep="\t",append=TRUE) #transition probabilities
        cat(hmmRes$alpha,"\n",file=outTxt[i],sep="\t",append=TRUE) #alpha neg.bin. parameters
        cat(hmmRes$beta,"\n",file=outTxt[i],sep="\t",append=TRUE) #beta neg.bin. parameters
        cat(hmmRes$logl,"\n",file=outTxt[i],sep="\t",append=TRUE) #loglikelihood max(hmmRes$logl[hmmRes$logl<0])
        cat(hmmRes$states,"\n",file=outTxt[i],sep="\t",append=TRUE) #inferred states
        cat(hmmRes$postprob,"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE) #posterior probability matrix
        propStates <- colSums(hmmRes$postprob)/sum(hmmRes$postprob)
        ##cat(propStates,"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE) #posterior probability proportion of HMM states
        cat(windTable[,1],"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE) #starting loci for windows
        cat(windTable[,2],"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE) #ending loci for windows
        cat(hmmRes$count[,1],"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE)
        if(strcmp(useGeno,"yes"))
            cat(V$y,"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE) # inferred ploidy for the above defined window
        if(!strcmp(useGeno,"yes"))
            cat(hmmOut,"\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE)
        
        ##cat(sum( (V$y - truePl)!=0 )/length(V$y)  ,"\n",file=paste(outTxt[i],"Error",sep=""),sep="",append=TRUE) #error rate when comparing to simulations (only for article publication)
        ##print(hmmRes$geno)
        
        ##plot ploidy inference    
        stringPlot <- sprintf("\tInferred ploidies from %s\nSample: %s",BASENAMEFILE[i], inputNames[whichInd])
        hmmRes$'chrSNP' <- chrSNP
        if(strcmp(useGeno,"yes")){
            cat( "Inferred state sequence: ", V$y, "\n", sep=" ")
            hmmPlotting(hmmRes, V, truePl=NA, main=stringPlot, propStates=propStates, CNV=rep(FALSE, length(V$y)))
            ##print on screen    
            cat(sprintf("\tInferred ploidies from %s. Sample: %s\n", BASENAMEFILE[i], inputNames[whichInd]))
            cat("\t-----------------------------------------------------\n")
            fileCounter <- fileCounter + 1
        }
        if(!strcmp(useGeno,"yes")){
            cat("Inferred state sequence: ", hmmOut, "\n", sep=" ")
            hmmPlotting(hmmRes, V=list(y=hmmOut, nu=matrix(0,length(hmmOut),length(unique(hmmOut)))), truePl=NA, main=stringPlot, propStates=propStates, CNV=rep(FALSE, length(hmmOut))) #add loci from windows
            ##print on screen    
            cat(sprintf("\tInferred ploidies from %s. Sample: %s\n", BASENAMEFILE[i], inputNames[whichInd]))
            cat("\t-----------------------------------------------------\n")
            fileCounter <- fileCounter + 1
        }
    }
    ##close pdf plot connection
    dev.off()
}

