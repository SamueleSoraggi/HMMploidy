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
           fileList = NA, #list of basenames for GUNZIPPED .genolike, .mafs and .par files
           wind = 100, #size of window for depth and genotype likelihoods. we work on a chromosome-basis.
           #so I have not implemented a list of chromosomes/loci intervals to read. it might come if we decide it is smart.
           maxPloidy = 6, #maximum ploidy. Must change with choice of potential ploidies (e.g. haploid might be excluded a priori by users)
           minInd = 1, #min ind having reads #also in ANGSD, better to have it here as well
           chosenInd = NA, #which Individual to consider (one at the time for now)
           #Must implement option for: all ind together or one at the time.
           #isSims = FALSE, #data is simulated. if TRUE fileList will also refer to file(s) with true ploidies
           alpha = NA, #alpha parameters comma separated
           beta = NA, #beta parameters comma separated
           quantileTrim ="0,1", #quantiles for trimming
           eps = .0005 #effect of sequencing and mapping error
           )

#if no argument aree given prints the need arguments and the optional ones with default
des<-list(fileList="[string] list of .genolike and phat.mafs and eventual .windows files",
          wind="[integer] Size of window for depth and genotype likelihoods (NA)",
          minInd="[integer] min Nr of individuals per locus having data (1)",
          maxPloidy="[integer] Maximum ploidy allowed (6)", #have to implement case where ploidies are chosen
          chosenInd ="[integers] which Individual to consider. one at the time for now. (NA=all)",
          isSims="[bool] data is simulated using the simulation script (FALSE)",
          alpha ="[numerics] alpha parameters comma separated (NA, read from .par file)",
          beta="[numerics] beta parameters comma separated (NA, read from .par file)",
          quantileTrim="[integers] comma-separated quantiles for trimming (0,1)",
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

angsdVector <- c(); fileVector <- c(); outPdf <- c(); outTxt <- c();
BASENAMEFILE <- c();
for(i in 1:length(filez)){
    fileVector[i] <- paste(filez[i],".genolikes",sep="")
    angsdVector[i] <- paste(filez[i],".mafs",sep="")
    splittedName <- unlist(strsplit(filez[i],split="/"))
    BASENAMEFILE[i] <- splittedName[length(splittedName)]
    outPdf[i] <- paste(filez[i],".pdf",sep="")
    outTxt[i] <- paste(filez[i],".hiddenMarkovPloidy",sep="")
}

#read parameters from file if not provided as input directly
#params <- list()
directInputPar = FALSE
#if(is.na(alpha) | is.na(beta)){
#    for(i in 1:length(filez))
#        params[[i]] <- read.table(paste(filez[i],".par",sep=""), header=FALSE, as.is=T)
#}
if(!(is.na(alpha) | is.na(beta))){
    directInputPar = TRUE
#    for(i in 1:length(filez))
#        params[[i]] <- rbind(alpha,beta)
}



##numeric conversion of inputs
wind <- as.numeric(wind)
minInd <- as.numeric(minInd)
maxPloidy <- as.numeric(maxPloidy)
#flags
isNumericChosenInd <- all(!is.na(chosenInd)) #check for choice of individuals

##individuals chosen for analysis (one by one at the moment)
##if chosenInd id NA it will be assigned as all individuals later
if(isNumericChosenInd)
    chosenInd <- eval( parse( text=paste("c(",chosenInd,")",sep="") ) )
if(!(is.na(alpha) | is.na(beta))){
    alpha <- eval( parse( text=paste("c(",alpha,")",sep="") ) )
    beta <- eval( parse( text=paste("c(",beta,")",sep="") ) )
}

#print(alpha)
#print(beta)

##print Rscript input
cat("----------\nfileList: ", fileList, " wind: ", wind," minInd: ", minInd, " chosenInd: ", chosenInd ," maxPloidy: ", maxPloidy, " alpha: ", alpha, " beta: ", beta, " quantileTrim: ", quantileTrim, " eps: ", eps,  "\n-----------\n" )

############################
### Supporting functions ###
############################


hmmPlotting <- function(hmm, V, truePl=NULL, main="Inferred ploidies"){
    options(warn=-1)
    loci = hmm$lociSNP
    borderVal <- round( seq(min(loci),max(loci),length.out=min(20,length(V$y)) ) )
    xlabels=c()
    if(max(loci)>=1e+6){
        XLAB="Scaffold Position (Mb)"
        for(i in 1:(length(borderVal)))
            xlabels[i] <- sprintf("%.1f", borderVal[i]/(1e+6))
    }
    if(max(loci)<1e+6){
        XLAB="Scaffold Position (Kb)"
        for(i in 1:(length(borderVal)))
            xlabels[i] <- sprintf("%.1f", borderVal[i]/(1e+3))
    }
    #print(xlabels)
    
    layout(matrix(c(1,1,1,1,2,2), nrow = 3, ncol = 2, byrow = TRUE))
    
    plot( V$y, pch=15, lwd=.75, col="navyblue", main=main, xaxt="n", yaxt="n", ylab="Ploidy", xlab=XLAB, ylim=c(min(V$y,truePl)-.5, max(V$y,truePl)+1 ), cex=.5, cex.main=1.4, cex.lab=1.2)
    
    if(!is.null(truePl))
        points(truePl-.075 , pch=15, lwd=.75, col="coral", cex=.5)
   
    abline( h=seq( min(V$y,truePl), max(V$y,truePl) ), col="gray" )

    axis( side=1, at=seq(1,length(V$y),length.out=min(20,length(V$y))), labels=xlabels, las=2, cex=1 )
    axis( side=2, at=seq(min(V$y), max(V$y)) )
    

    postProb <- hmm$postprob
    counter=1
    for(yValue in intersect(hmm$states,seq(min(V$y,truePl), max(V$y,truePl)))){
        polygon( x=c( length(V$y),1, seq(1,length(V$y)), length(V$y) ), y= yValue + 0.025 + c( 0, 0, postProb[,counter], 0 )/2, col="deepskyblue1", border=NA)
        counter=counter+1
    }
    
    legendCol = c("coral","navyblue","deepskyblue1")
    legendPch = c(15,15,15)
    legendTxt = c( "True Ploidy", "Inferred Ploidy", "Posterior Prob." )
    if(is.null(truePl)){
        legendCol = c("navyblue","deepskyblue1")
        legendPch = c(15,15)
        legendTxt = c( "Inferred Ploidy", "Posterior Prob." )
    }
    
    legend(x=1, y = max(V$y,truePl)+1, legend=legendTxt, col = legendCol, lwd = rep(3,length(legendPch)), pch = legendPch, bty = "n", ncol = length(legendPch), cex=1.4)

    ##PLOT DEPTH
    plot(hmm$count[,1], type="p", lwd=2, col="deepskyblue1", xlab=XLAB, main="Window Mean Depth", xaxt="n", ylab="Window Mean Depth", bg=3, ylim=c(min(hmm$count[,1])-.1*min(hmm$count[,1]),max(hmm$count[,1])+.2*max(hmm$count[,1])), cex.lab=1.2, cex.main=1.3)
    axis( side=1, at=seq(1,length(V$y),length.out=min(20,length(V$y))), labels=xlabels, las=2, cex=1 )
    abline(h=hmm$mu, col="coral")
    legend(x=1, y = max(hmm$count[,1])+.25*max(hmm$count[,1]), legend=c("Mean Depth", "Distribution Mean"), col = c("deepskyblue1","coral"), lwd = rep(3,3), lty=c(NA,1), pch = c(20,NA), bty = "n", ncol = 2, cex=1.4)
options(warn=0)
}


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

# nbh_init  Initalize parameters for nbh_em
# 			Function nbm_em (NB mixture model) is used to find alpha,
# 			beta, and wght (mixprop); wght (1xN) is repeated N times row-wise
# 			to represent the initial TRANS for the subsequent nbh_em training
# 			Use: nbh0 <- nbh_init(count, K)
#			nbh0: list(TRANS, alpha, beta)

nbm_em <- function(count, alpha, beta, wght, NBM_NIT_MAX=250, NBM_TOL=1e-2){
	
	# Data length
	Total <- length(count)	
	#if(any(count < 0) || any(count != round(count))){
	#	stop("Data does not contain positive integers.")
	#}
	
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
    }
				
    llk <- log(sum(forwrd2[Total,])) + sum(log(scale2))

    #bckwrd[Total,] <- matrix(1, ncol=N)
    bckwrd2[Total,] <- matrix(1, ncol=N)
    
    for(t in (Total-1):1) {			
        #bckwrd[t,] <- (bckwrd[t+1,] * dens[t+1,]) %*% t(TRANS)
        bckwrd2[t,] <- (bckwrd2[t+1,] * dens2[t+1,]) %*% t(TRANS)
        #bckwrd[t,] <- bckwrd[t,] / scale[t]
        bckwrd2[t,] <- bckwrd2[t,] / scale2[t]
    }
	
    ni <- forwrd2 * bckwrd2
    ni <- ni / ( apply(ni, 1, sum) %*% matrix(1,ncol=N) )
    #bckwrd=bckwrd2
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
nbHMM <- function(count, delta, TRANS, alpha, beta, genolike=0, ws=1, PLOIDYMAX=maxPloidy, NBH_NIT_MAX=10000, NBH_TOL=1e-5, MAXALPHA=1e7, MAXBETA=1e7){

    stateVec <- 1:PLOIDYMAX
    geno <- genolike[,stateVec]
    bicIter <- 0
    Total <- dim(count)[1]	
    if(any(count < 0)) stop("Data does not contain positive integers.")
    
    N <- nrow(TRANS)
    alpha <- matrix(alpha, 1, N)
    beta <- matrix(beta, 1, N)	

    logl <- matrix(0, ncol=NBH_NIT_MAX)
    
    alpha0 <- alpha
    beta0 <- beta
    TRANS0 <- TRANS
    postprob0 <- matrix(1, Total, N)/N
    logl0 <- logl
    rate <- +Inf
    checkBIC <- TRUE
    
    ##Main loop of the EM algorithm
    for(nit in 1:NBH_NIT_MAX) {

        bicIter <- bicIter + 1

        ##check after some steps for some convergence and try to remove 1 state
        if(rate<=.0001 && bicIter>=30 && checkBIC){
            print(alpha)
            if(length(alpha)==2){ #if removing one state leaves only one ploidy
                alpha=matrix(alpha,ncol=length(alpha))
                beta=matrix(beta,ncol=length(beta))
                cat("\t==>reduction to 1 state start\n")
                compLL <- c()
                for(kk in 1:length(stateVec)){#try out one state at a time
                    diff <- +Inf
                    contDiff <- 0
                    gg <- geno[,kk]
                    aa <- alpha[,kk]; aa <- matrix(aa,ncol=length(aa))
                    bb <- beta[,kk]; bb <- matrix(bb,ncol=length(bb))
                    nanFlag <- FALSE
                    
                    while(diff >= 0.001 & contDiff < 150){
                        contDiff <- contDiff + 1
                        resSingle1 <- MStepSingle(count,aa,bb,gg)
                        llk1 <- resSingle1$llk
                        if(is.nan(E1$llk)){
                            nanFlag=TRUE
                            break
                            }
                        aa <- resSingle1$alpha
                        bb <- resSingle1$beta
                        contDiff <- contDiff + 1
                        resSingle2 <- MStepSingle(count,aa,bb,gg)
                        llk2 <- resSingle2$llk
                        if(is.nan(E1$llk)){
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
                while(diff >= 0.001 & contDiff < 150){
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

                    diff <- abs( (LLtest[1] - LLtest[2])/LLtest[1] )
                    if(is.nan(diff))
                        diff=0
                }

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
            print(alpha2)
            print(beta2)
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
            
            while(diff >= 0.001 & contDiff < 250){
                
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
        #message(sprintf('Iter %d:\tLLK=%.3f\tploidies=%d', (nit-1), logl[nit], dim(alpha)[2] ) )	
        if(is.nan(logl[nit])) {			
            warning("NaN logl data detected. Returning the previous training results")
            logl[nit] <- 0; TRANS=TRANS0; alpha=alpha0; beta=beta0;  logl=logl0; bckwrd=postprob0; dens=dens0
            break
        }
  
        M <- MStep(E,count,TRANS,alpha,beta)
        TRANS <- M$TRANS
        alpha <- M$alpha
        beta <- M$beta
        bckwrd <- E$ni
        dens <- E$dens2
        delta <- M$delta

        
        if(any(is.na(M$alpha))) {
            warning(sprintf("Updated alpha becomes NA probably %s","due to bad initial alpha or insuff. data"))
            TRANS=TRANS0; alpha=alpha0; beta=beta0; logl=logl0; bckwrd=postprob0; dens=dens0
            #break
        }
        if(any(is.na(M$beta))) {
            warning(sprintf("Updated beta becomes NA probably %s","due to bad initial alpha or insuff. data"))
            TRANS=TRANS0; alpha=alpha0; beta=beta0; logl=logl0; bckwrd=postprob0; dens=dens0
            #break
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

                                        
    list(count=count, delta=delta, TRANS=TRANS, alpha=alpha, beta=beta, logl=logl, postprob=bckwrd, dens=dens, mu = alpha/beta, sigma=alpha/beta+alpha/beta^2, states=stateVec)				
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
        return(list(delta=delta,TRANS=TRANS,alpha=alpha,beta=beta))
    
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
sumGeno <- function(dp,ws=1,loci,lociSNP=loci,avg=FALSE){   
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

##sum of logarithm on rows of a matrix
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
   

##Likelihood of f=data vector given genotype. gl=genotype likelihoods vector. h=inbreeding coefficient.
pGenoData <- function(f,gl,nInd=1,h=0){   
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


pGenoDataSingle <- function(f,gl,h=0){   
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
        stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])
    return(list(y=y,nu=nu))
}



wind <- as.numeric(wind) #window size
minInd <- as.numeric(minInd)
maxPloidy <- as.numeric(maxPloidy) 
#nInd=as.numeric(nInd) #number of individuals
#GLsingle #genolikes on ploidies (all and single individual)
#DPsingle #depths (all and single individual)
#hmmRes #resulting list containing results from files
#V;
#whichInd <- chosenInd[1] #I will make a loop over this afterwards

#################################################
### Begin file-by-file analysis #################
#################################################

for(i in 1:length(fileVector)){
    
    cat("==> Analyze ", filez[i], "\n",sep="")
    ##read in the data from .mafs and .genolikes files
    GL <- fread(input=fileVector[i],sep="\t",showProgress=TRUE,header=FALSE,data.table=FALSE)
    #FREQFILE <- fread(input=angsdVector[i],sep="\t",showProgress=TRUE,header=TRUE,data.table=FALSE)
    rowsGL <- dim(GL)[1]
    nInd <- length( unique( GL[,3] ) )
    sites <- unique( GL[ ,2] )    
    DP <- GL[ ,5]


    #calculate allele frequencies
    majorReads <- GL[,8]
    minorReads <- GL[,9]
    freqs <- alleleFrequencies(majorReads,minorReads,nInd,minInd,eps)
    
    GL <- GL[ ,-c(1:9)]
    
    if(!isNumericChosenInd)
        chosenInd <- 1:nInd

    ###############################
    ## begin of FOR loop to      ##
    ## read one genome at a time ##
    ###############################

    fileCounter = 1
    ##open pdf plot connection
    pdf(outPdf[i])

    for(whichInd in chosenInd){
        
        #if(directInputPar==FALSE){
        #    alpha=as.vector(as.numeric(params[[i]][2*whichInd - 1, 1:maxPloidy ]))
        #    beta=as.vector(as.numeric(params[[i]][2*whichInd, 1:maxPloidy]))
        #}
    
    ##select single individual depth and genolikes
        idxSingle <- seq(whichInd,rowsGL,nInd)
        DPsingle <- DP[idxSingle]; GLsingle <- GL[idxSingle, ] #individual depth/genolikes
    ##trim the depth at the chosen quantile
        quantiles <- eval( parse( text=paste("c(",quantileTrim,")",sep="") ) )
        q <- quantile( DPsingle, quantiles )  
        idx <- which( DPsingle<=as.numeric(q[2]) & DPsingle>=as.numeric(q[1]) )
        DPsingle <- DPsingle[idx] #individual filtered datad
        GLsingle <- GLsingle[idx, ] #...""
        sitesIndiv <- sites[idx] #......""
        freqsIndiv <- freqs[idx] #......""
        idxTot = as.vector( sapply(idx, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        GLfiltered <- GL[idxTot, ] #all data filtered
        DPfiltered <- DP[idxTot] #......""

        ##find SNPs with thresholds .1<f<.9
        findSNP <- which(freqsIndiv>.1 & freqsIndiv<.9)
        freqsSNP <- freqsIndiv[findSNP]
        sitesSNP <- sitesIndiv[findSNP]
        totSNP <- as.vector( sapply(findSNP, function(j) ((j-1)*nInd+1):(j*nInd) ) )
        ##DPSNP = DPsingle[totSNP] I think it is not needed

        ##THIS CALCULATION IS NEEDED ONLY THE FIRST TIME (fileCounter==1).    
        ##probability of data given genotype, ploidy and frequencies (per SNP)...
        if(fileCounter==1){
            geno2 <- matrix(0, nrow=maxPloidy, ncol=length(freqsSNP))
            for(pp in 1:maxPloidy) #change ploidy
                geno2[pp,] <- pGenoData( f=freqsSNP, gl=readGL( findSNP, pp, nInd=nInd, GLfiltered ), nInd=nInd )
            ##...and per window
            geno <- apply( geno2, 1, function(x) sumGeno(x,wind,sitesIndiv,sitesSNP) )
        }
        ##mean depth over each locus in a window (not only on SNPs)
        ##use sumGeno(DPsingle,wind,sitesIndiv,sitesSNP,avg=TRUE)
        ##to apply only the average on SNPs (very noisy result)
        DPmean <- sumGeno( DPsingle, wind, sitesIndiv, sitesIndiv, avg=TRUE )

        ##clean from NA, NaN or infinite values
        keepSites <- apply( geno, 1, function(x) sum(is.na(x) | is.nan(x) | is.infinite(x))==0 )   
        DPmean = DPmean[keepSites]
        geno = geno[keepSites, ]
        keepSites <- which( !is.na(DPmean) & !is.nan(DPmean) & !is.infinite(DPmean) )
        DPmean = DPmean[keepSites]
        geno = geno[keepSites,]

        ##rescale likelihood of the data (avoids underflow)
        genoResc <- t( apply( geno , 1, logRescale ) )

        ##some initial parameters
        delta=rep(1/maxPloidy,maxPloidy) #i think it is ok without prior info
        Pi0=matrix(1/maxPloidy,nrow=maxPloidy,ncol=maxPloidy) #tridiagonal makes more sense?
        count <- matrix(DPmean,ncol=1)
        if(directInputPar==FALSE){
            cat(sprintf("    Estimate parameters for %d states\n", maxPloidy))
            alpha <- tail(quantile(count[count[,1]>0,1], probs=seq(0, 1, 1/maxPloidy), na.rm=TRUE), maxPloidy)
            beta <- rep(1, maxPloidy)	
            delta <- rep(0.5, maxPloidy)
            nbm <- nbm_em(count, alpha, beta, delta)	
            alpha <- nbm$alpha	
            beta <- nbm$beta	
            Pi0 <- matrix(rep(nbm$wght, maxPloidy), ncol=maxPloidy, byrow=TRUE)
	    # Order the parameters s.t. they are in increasing order of mean values
            myorder <- order(alpha/beta)	
            alpha <- matrix(alpha[myorder],ncol=length(myorder))
            beta <- matrix(beta[myorder],ncol=length(myorder))
            Pi0 <- Pi0[myorder, myorder]		
        }
        
        cat("    N.samples ",nInd," alpha0: ",alpha," beta0: ",beta,"\n",sep=" ")        
        ##start the HMM
        hmmRes <- nbHMM(count, alpha=alpha, beta=beta, TRANS=Pi0, delta=delta, genolike=genoResc, PLOIDYMAX=maxPloidy)
        hmmRes$'lociSNP' = sitesSNP; hmmRes$'geno' = geno;

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

        ## TO DO
        ##return something in a file
        cat("File: ",BASENAMEFILE[i]," individual  ",whichInd," out of ",nInd,"\n",sep="",file=outTxt[i],append=!(fileCounter==1))
        cat(hmmRes$delta,"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(hmmRes$TRANS,"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(hmmRes$alpha,"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(hmmRes$beta,"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(max(hmmRes$logl[hmmRes$logl<0]),"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(hmmRes$states,"\n",file=outTxt[i],sep="\t",append=TRUE)
        cat(hmmRes$postprob,"\n\n",file=outTxt[i],sep="\t",append=TRUE,fill=FALSE)
        
        
    ##plot ploidy inference    
        stringPlot <- sprintf("\tInferred ploidies from %s\nindividual %d", BASENAMEFILE[i], whichInd)
        hmmPlotting(hmmRes, V, truePl=NULL, main=stringPlot)
    
    ##print on screen    
        cat(sprintf("\tInferred ploidies from %s individual %d\n", BASENAMEFILE[i], whichInd))
        cat("\t-----------------------------------------------------\n")
        fileCounter <- fileCounter + 1
    }

    ##close pdf plot connection
    dev.off()
}
