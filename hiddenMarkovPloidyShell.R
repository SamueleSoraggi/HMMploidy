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
           wind = NA, #size of window for depth and genotype likelihoods. we work on a chromosome-basis.
           #so I have not implemented a list of chromosomes/loci intervals to read. it might come if we decide it is smart.
           maxPloidy = 6, #maximum ploidy. Must change with choice of potential ploidies (e.g. haploid might be excluded a priori by users)
           minInd = 1,#min ind having reads #also in ANGSD, better to have it here as well
           chosenInd = NA,#which Individual to consider (one at the time for now)
           #Must implement option for: all ind together or one at the time.
           isSims = FALSE, #data is simulated. if TRUE fileList will also refer to file(s) with true ploidies
           alpha = NULL, #alpha parameters comma separated
           beta = NULL #beta parameters comma separated
           )

#if no argument aree given prints the need arguments and the optional ones with default
des<-list(fileList="list of .genolike and phat.mafs and eventual .windows files",
          wind="Size of window for depth and genotype likelihoods (NA)",
          minInd="min Nr of individuals per locus having data (1)",
          maxPloidy="Maximum ploidy allowed (6)", #have to implement a case where ploidies are chosen
          chosenInd = "which Individual to consider. one at the time for now. (NA=all)",
          isSims="data is simulated using the simulation script (FALSE)",
          alpha ="alpha parameters comma separated",
          beta="beta parameters comma separated"
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


##numeric conversion of inputs
wind <- as.numeric(wind)
minInd <- as.numeric(maxPloidy)
maxPloidy <- as.numeric(maxPloidy)
##individuals chosen for analysis (one by one at the moment)
##if chosenInd id NA it will be assigned as all individuals later
if(all(!is.na(chosenInd)))
    chosenInd <- eval( parse( text=paste("c(",chosenInd,")",sep="") ) )

alpha <- eval( parse( text=paste("c(",alpha,")",sep="") ) )
beta <- eval( parse( text=paste("c(",beta,")",sep="") ) )

##print Rscript input
cat("----------\nfileList: ", fileList," wind: ", wind," minInd: ", minInd, " chosenInd: ", chosenInd ," maxPloidy: ", maxPloidy, " isSims: ", isSims,  " alpha: ", alpha, " beta: ", beta, "\n-----------\n" )

############################
### Supporting functions ###
############################


hmmPlotting <- function(hmm, V, truePl=NULL, main="Inferred ploidies"){
    border = list()
    loci = hmm$lociSNP
    for(i in 1:length(loci))
        border <- round( seq(1,length(loci),length.out=10) )
    idxBorder = c()
    for(i in 1:length(border))

   
    layout(matrix(c(1,1,1,1,2,2), nrow = 3, ncol = 2, byrow = TRUE))
 
    plot( V$y, pch=15, lwd=.75, col="navyblue", main=main, xaxt="n", yaxt="n", ylab="Ploidy", xlab="Scaffold Position (Kb)", ylim=c(min(V$y,truePl)-.5, max(V$y,truePl)+.75 ), cex=.5, cex.main=1.5, cex.lab=1.2)

    if(!is.null(truePl))
        points(truePl-.075 , pch=15, lwd=.75, col="coral", cex=.5)

    xlabels=c()
    for(i in 1:length(border))
        xlabels <- c(xlabels,  sprintf("%.1f",loci[ border[i] ]/1e+3) )
    
    abline( h=seq( min(V$y,truePl), max(V$y,truePl) ), col="gray" )

    axis( side=1, at=seq(1,length(V$y),length.out=10), labels=xlabels, las=0, cex=1.2 )
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
  
    legend(x=1, y = max(V$y,truePl)+.85, legend=legendTxt, col = legendCol, lwd = rep(3,length(legendPch)), pch = legendPch, bty = "n", ncol = length(legendPch), cex=1.5)

    ##PLOT DEPTH
    plot(hmm$count[,1], type="p", lwd=2, col="deepskyblue1", xlab="Scaffold Position (Mb)", main="Mean Depth", xaxt="n", ylab="Mean Depth", bg=3, ylim=c(min(hmm$count[,1])-.1*min(hmm$count[,1]),max(hmm$count[,1])+.125*max(hmm$count[,1])), cex.lab=1.2)
    axis( side=1, at=seq(1,length(V$y),length.out=10), labels=xlabels, las=0, cex=1.2 )
    abline(h=hmm$mu, col="coral")
    legend(x=1, y = max(hmm$count[,1])+.1*max(hmm$count[,1]), legend=c("Mean Depth", "Neg.Bin. Mean"), col = c("deepskyblue1","coral"), lwd = rep(3,3), lty=c(NA,1), pch = c(20,NA), bty = "n", ncol = 2, cex=1.5)

cat("==> Plot Done :)")

}


##fast function to match loci between ANGSD .maf file (freq) and .genolikes file (refer)
cppFunction('NumericVector matchSites(NumericVector freq, NumericVector refer) {
  int n = refer.size();
  NumericVector out(n);
  int found = 0;
  int j = 0;

  for (int i = 0; i < n; i++) {
    found = 0;
    while (found == 0) {
      if (refer[i] == freq[j]){
         out[i] = j;
         found = 1;
       }
      if (refer[i] < freq[j]){
         found = 1;
       }
    j++;
    } 
  }
  return out;
}')


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
    
    while (any(tmp <= 0)){
        warning(sprintf("Alpha (%.4f)<0 ! Try smaller (10%s) Newton step ...\n", tmp_step,"%"))
        tmp_step <- tmp_step/10
        tmp <- alpha + tmp_step
    }
		
    alpha <- tmp

    return(list(delta=delta,TRANS=TRANS,alpha=alpha,beta=beta,remStates=remStates))
}

##Mstep for one-state conditional EM optimization
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

    
    densLog <-  matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha)) - count %*% log(1+beta) + lgamma(count + alpha) - dnorm  + geno

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


##EM algorithm for HMM optimization
nbHMM <- function(count, delta, TRANS, alpha, beta, genolike=0, ws=1, PLOIDYMAX=6, NBH_NIT_MAX=10000, NBH_TOL=1e-5, MAXALPHA=1e7, MAXBETA=1e7){

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
            if(dim(alpha)[2]==2){ #if removing one state leaves only one ploidy
                compLL <- c()
                for(kk in 1:length(stateVec)){#try out one state at a time
                    diff <- +Inf
                    contDiff <- 0
                    gg <- geno[,kk]
                    aa <- alpha[,kk]
                    bb <- beta[,kk]
                    
                    while(diff >= 0.001 & contDiff < 150){
                        contDiff <- contDiff + 1
                        resSingle1 <- MStepSingle(count,aa,bb,gg)
                        llk1 <- resSingle1$llk
                        aa <- resSingle1$alpha
                        bb <- resSingle1$beta
                        contDiff <- contDiff + 1
                        resSingle2 <- MStepSingle(count,aa,bb,gg)
                        llk2 <- resSingle2$llk
                        aa <- resSingle2$alpha
                        bb <- resSingle2$beta
                        diff <- abs( (llk1 - llk2)/llk1 )
                    }

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

                cat("Best newBIC ", newBIC," oldBIC", oldBIC, "reduce states ", (newBIC>oldBIC), "\n") 

                if(oldBIC > newBIC){
                    checkBIC <- FALSE
                    break
            }
            else{
                geno <- geno[,whichBIC ]
                alpha <- alpha[,whichBIC]
                beta <-  beta[,whichBIC]
                stateVec <- stateVec[ whichBIC ]
                bckwrd <- matrix(1,nrow(count),1)
            }
                cat("!-!-!- States Relation : ", stateVec,"\n")
                cat("!-!-!- alpha : ", alpha, " beta : ", beta,"\n")
                break
        }

            
            K <- dim(alpha)[2]
            
            ##if remoiving one state leaves at least other two, then what follows will happen
            cat("check reduction to ",K-1,"states start\n\talpha: ", as.vector(alpha),"\n\tbeta: ",as.vector(beta),"\n",sep=" ")
            bicIter <- 0
            viterbiIter <- 0
            combMat <- combs( 1:K, K-1 )
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
            cat("Best newBIC ", newBIC," oldBIC", oldBIC, "reduce states ", (newBIC>oldBIC), "\n") 

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
            cat("!-!-!- States Relation : ", stateVec,"\n")
            cat("!-!-!- alpha : ", alpha, "\nbeta : ", beta, "\nmu : ", alpha/beta, "\n")
        }

        N <- nrow(TRANS)        
        E <- EStep(count,delta,TRANS,alpha,beta,geno)

        logl[nit] <- E$llk
        message(sprintf('Iter %d:\tLLK=%.3f\tploidies=%d', (nit-1), logl[nit], dim(alpha)[2] ) )	
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
minInd <- as.numeric(maxPloidy)
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
    FREQFILE <- fread(input=angsdVector[i],sep="\t",showProgress=TRUE,header=TRUE,data.table=FALSE)
    rowsGL <- dim(GL)[1]
    nInd <- length( unique( GL[,3] ) )
    sites <- unique( GL[ ,2] )
    DP <- GL[ ,5]
    TRUEREF <- GL[seq(1,rowsGL,nInd),6]
    GL <- GL[ ,-c(1:7)]
    FQInd <- FREQFILE[ ,7] #per-loci Nr individuals in frequency estimate
    FF <- FREQFILE[ ,2] #sites to be kept...
    FF <- FF[FQInd>=minInd] #...filtered using FQInd
    sameSites <- matchSites(FF,sites) + 1 #Format matching with .mafs output (+1=adapting indexing RC++)
    FQ <- FREQFILE[sameSites,6] 
    FQREF <- FREQFILE[sameSites,3] #reference in the .mafs file
    sameREF <- (TRUEREF==FQREF) #reference matching with the one in .genolikes
    freqs <- as.numeric(FQ)
    freqs[!sameREF] <- 1-freqs[!sameREF] #for freq f of non-matching reference allele, do 1-f
    if(any(is.na(chosenInd)))
        chosenInd <- 1:nInd
    
    ###############################
    ## begin of FOR loop to      ##
    ## read one genome at a time ##
    ###############################

    fileCounter = 1
    ##open pdf plot connection
    pdf(outPdf[i])

    for(whichInd in chosenInd){
        ##select single individual depth and genolikes
        idxSingle <- seq(whichInd,rowsGL,nInd)
        DPsingle <- DP[idxSingle]; GLsingle <- GL[idxSingle, ]
    
        ##trim the depth at the .1 and .9 quantile
        TRIM=TRUE
        if(TRIM){
            q <- quantile( DPsingle, c(.1,.9) )  
            idx <- which( DPsingle<as.numeric(q[2]) & DPsingle>as.numeric(q[1]) )
            DPsingle <- DPsingle[idx]
            GLsingle <- GLsingle[idx, ]
            FQ <- FQ[idx]
            sitesIndiv <- sites[idx] #individual filtered sites
            freqsIndiv <- freqs[idx] #individual filtered frequencies
            idxTot = as.vector( sapply(idx, function(j) ((j-1)*nInd+1):(j*nInd) ) )
            GLfiltered <- GL[idxTot, ] #remember the filtering of GL
            DPfiltered <- DP[idxTot]
        }

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
            ##geno2[pp,] <- pGenoData(f=freqsSNP, gl=readGL( findSNP, pp, nInd=1, GLsingle[[i]]),nInd=1)      
            ##...and per window
            geno <- t( apply( geno2, 1, function(x) sumGeno(x,wind,sitesIndiv,sitesSNP)) )
            geno <- t( geno )
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
        delta=rep(1/6,6) #i think it is ok without prior info
        Pi0=matrix(1/6,nrow=6,ncol=6) #tridiagonal makes more sense?
        count <- matrix(DPmean,ncol=1)
        ##start the HMM
        hmmRes <- nbHMM(count, alpha=alpha, beta=beta, TRANS=Pi0, delta=delta, genolike=genoResc)
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
    

    ##plot ploidy inference    
    stringPlot <- sprintf("Inferred ploidies from\n%s\nfor individual %d", BASENAMEFILE[i], whichInd)
    hmmPlotting(hmmRes, V, truePl=NULL, main=stringPlot)
    
    ##print on screen    
    cat(sprintf("Inferred ploidies from\n%s\nfor individual %d\n", BASENAMEFILE[i], whichInd))
    fileCounter <- fileCounter + 1
    }

    ##close pdf plot connection
    dev.off()
}
