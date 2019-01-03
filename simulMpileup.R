
# simulate mpileup file with different ploidy 1-5, assuming only 2 alleles and genotype frequencies defined by K and Ne and F

# assume, population allele frequencies drawn from an exponential distribution

#seq1 272 T 24 ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&

library("getopt")

# http://www.inside-r.org/packages/cran/getopt/docs/getopt.package
spec=matrix(c(
	      'out',	'o', 2, "character", "output files for real data (and log if verbose), (mpileup is in stdout)",
	      'copy', 	'c', 1, "character", "ploidy per sample, e.g. 2x3,4 is 2,2,2,4",
	      'sites', 	's', 2, "integer", "number of sites [default 1,000]",
	      'depth', 	'd', 2, "double", "mean haploid depth per sample [default 20.0]",
	      'lendepth', 'l', 2, "integer", "mean length of sites with increasing/decreasing depth [default 0, disabled]",
	      'errdepth', 'e', 2, "double", "error rate in mean depth [default 0.05]",
	      'qual', 'q', 2, "integer", "mean base quality in phred score [default 20]",
	      'pvar', 'r', 2, "double", "probability that site is variabile in the population [1.0]",
	      'ksfs', 'k', 2, "double", "coeff. for shape of SFS default [1.0]",
	      'panc', 'a', 2, "double", "probability that ancestor state is correct [1.0]",
	      'ne', 'n', 2, "integer", "effective population size [default 10,000]",
	      'pool', 'p', 0, "logical", "enable pool data",
	      'help', 'h', 0, "logical", "print help message",
	      'verbose', 'v', 0, "logical", "verbose creates log file",
              'offset', 'f', 0, "integer", "offset value for genomic position",
              'seed','u', 2, "integer", "random seed for simulations reproducibility [default 180218]"
	      ), byrow=TRUE, ncol=5)
opt <- getopt(spec)

# help
if ( !is.null(opt$help) ) {
	write.table(spec, sep="\t", quote=F, col.names=F, row.names=F)
	q(status=1)
}

# default values
if (is.null(opt$sites)) opt$sites <- 1e3
if (is.null(opt$depth)) opt$depth <- 20.0
if (is.null(opt$qual)) opt$qual <- 20
if (is.null(opt$ksfs)) opt$ksfs <- 1.0
if (is.null(opt$panc)) opt$panc <- 1.0
if (is.null(opt$pvar)) opt$pvar <- 1.0
if (is.null(opt$ne)) opt$ne <- 1e4
if (is.null(opt$verbose)) opt$verbose <- FALSE
if (is.null(opt$pool)) opt$pool <- FALSE
if (is.null(opt$offset)) opt$offset <- 0
if (is.null(opt$lendepth)) opt$lendepth <- 0
if (is.null(opt$errdepth)) opt$errdepth <- 0.05
if (is.null(opt$seed)) opt$seed <- 180218

# set seed 
set.seed(opt$seed)

# switch panc to 1-panc for old consistency to previous version
opt$panc <- 1-opt$panc

# assign to old variables (then in the future change this)
fout_log <- paste(opt$out, ".log", sep="", collapse="")
fout_real <- opt$out
nsites <- opt$sites
mbqual <- opt$qual
K <- opt$ksfs
Ne <- opt$ne

# parse ncopy
ncopy <- c()
tmp <- strsplit(opt$copy, split=",")[[1]]
for (i in 1:length(tmp)) {
	tmp2 <- strsplit(tmp[i], split="x")[[1]]
	if (length(tmp2)>1) {
		ncopy <- c(ncopy, rep(tmp2[1], times=tmp2[2]))
	} else {
		ncopy <- c(ncopy, tmp2[1])
	}
}
rm(tmp); rm(tmp2)
ncopy=as.numeric(ncopy)
if (max(ncopy)>6) {
	cat("Max ploidy is 6.\n")
	q(status=1)
}

# init files
if (!is.null(opt$out)) cat("", file=fout_real)

# how many samples
nsams <- length(ncopy)

# write to log file
if (opt$verbose & !is.null(opt$out)) {
	cat("", file=fout_log)
	write.table(spec, sep="\t", quote=F, col.names=F, row.names=F, file=fout_log, append=T)
	cat(unlist(opt),"\nNr of samples:", nsams, "\nPloidies:", ncopy, "\n",file=fout_log, append=T)
	cat("Files:", fout_real, ",",fout_log,"\n", append=T, file=fout_log)
}

# sample depths and qualities, the latter are centered around phred score = 10

rangeLams <- rep(0,nsites)
sampledLams <- rep(0,nsites)
depth <- matrix( 0, nrow=nsams, ncol=nsites )
for(samIdx in 1:nsams){
    rangeLams <- opt$depth * ncopy[samIdx] + (opt$depth * ncopy[samIdx] * opt$errdepth * c(-1,1))
    sampledLams <- runif(nsites, min=rangeLams[1], max=rangeLams[2])
    depth[samIdx,] <- rpois(nsites,sampledLams)
    }

rm(sampledLams)
rm(rangeLams)

if (opt$lendepth>0) {

	conDepth <- matrix(NA, nrow=nrow(depth), ncol=ncol(depth))
	conDepth[,1] <- depth[,1]

	for (j in 1:nsams) {

            indexes <- c(1)
            toBeTaken <- 2:nsites

            i <- 2
            lengths <- rpois(nsites, opt$lendepth)
            increasing <- sample(x=c(0,1),size=nsites,prob=c(.5,.5),replace=TRUE)
            
		while (i <= nsites) {

			#lenSeg <- rpois(1, opt$lendepth)
                        lenSeg <- lengths[i]
			#increasing <- sample(c(0,1),1)

			ind <- c()
			if (increasing[i]) {
				ind <- toBeTaken[which(depth[j,toBeTaken]>=depth[j,i])[1:lenSeg]]
			} else {
				ind <- toBeTaken[which(depth[j,toBeTaken]<=depth[j,i])[1:lenSeg]]
			}
			ind <- ind[which(!is.na(ind))]
		
			if (length(ind)>0) {
				if (increasing) {
					ind <- ind[sort(depth[j,ind], dec=F, index.ret=T)$ix]
				} else {
					ind <- ind[sort(depth[j,ind], dec=T, index.ret=T)$ix]
				}

				indexes <- c(indexes, ind)
				toBeTaken <- setdiff(toBeTaken, ind)
				i <- length(indexes) + 1 
			}
		}
		conDepth[j,] <- depth[j,indexes]
	}
depth <- conDepth 
rm(conDepth)
}

#write.table(depth, sep="\t", quote=F, col.names=F, row.names=F)

# ascii, already starting at +33
pscores <- '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^_`abcdefghijklmnopqrstuvwxyz{|}~'

# assuming:
#ref=sample(c("A","C","G","T"),1)
#nonref=sample(setdiff(c("A","C","G","T"),ref),1)
# assume that ref is ancestral and nonref is derived, so at the population level
#major=ref
#minor=nonref
ref <- "A"
nonref <- "C"

# population allele frequencies
# if polymorphic
ee <- (1/(1:(Ne-1))^(1/K))
ee <-ee/sum(ee)
# the sum must be = pvar
ee <- ee*(opt$pvar)
# add that some sites might not be polymorphic in the populationder <- c(0, ee, 0)
pder <-c( (1-opt$pvar)*(9/10), ee, (1-opt$pvar)*(1/10) )

# this is the expected p
#Ep=weighted.mean(seq(0,Ne,1), pder)/Ne
#Ep

# this is prob of Major being the ancestral
# sum(pder[1:(floor(Ne/2)+1)]); 1- sum(pder[1:(floor(Ne/2)+1)])

# ascii phred score
# http://www.omixon.com/bioinformatics-for-beginners-file-formats-part-2-short-reads/

# if depth is 0, replace with only 1 read with very low quality

#for (i in valid) { # cycle across sites

# sample derived allele frequency
qqVector <- sample(x=seq(0,Ne,1),size=opt$sites,prob=pder,replace=TRUE)/Ne 
# probability of incorrectly assigning the ancestral state
pAncErr <- sample(x=c(0,1),size=opt$sites,prob=c(1-opt$panc,opt$panc),repl=TRUE)
qqVector[which(pAncErr==1)] <- 1-qqVector[which(pAncErr==1)]


for (i in 1:opt$sites) {

	# if pool, initialise
	pool_alls <- pool_bqs <- c()

	# first elements of line: chrom, pos, reference
	linea <- c(paste("copy_",opt$copy,sep="",collapse=""), i+opt$offset, ref)
	# for real data output
	linea_real <- c(linea, nonref)

	# count derived alleles and print on file
	daf <- 0

	# sample derived allele frequency
	#qq <- sample(x=seq(0,Ne,1),size=1,prob=pder)/Ne
        qq <- qqVector[i] 

        # probability of incorrectly assigning the ancestral state
	#if (sample(x=c(0,1),size=1,prob=c(1-opt$panc,opt$panc),repl=F)) qq <- 1-qq

	pp <- 1-qq
	linea_real <- c(linea_real, qq)

	for (n in 1:nsams) { # cycle across samples

		alls <- bqs <- c() # init bases and qualities

		# haploid case
		if (ncopy[n]==1) {

			# genotype probs assuming HWE
			priors <- c(pp, qq);

			# sample genotypes according to previously calculated probs
			genos <- c("A", "C")
			geno <- genos[sample(1:length(genos),1, prob=priors)]

			# daf
			if (geno=="C") daf <- daf+1

			if (depth[n,i]>0) { # if data

              		  	# sample base qualities for all reads around the mean
                		ibq <- round(rnorm(mean=mbqual,sd=2,n=depth[n,i]))
				ibq[which(ibq<0)] <- 0

				# for each read
				for (j in 1:depth[n,i]) {

					# base quality in ASCII character
					bqs <- c(bqs, substring(pscores, ibq[j]-1, ibq[j]-1))

					# from base quality calculate base probability
		      			ps <- ibq[j]
					p <- 10^(-(ps/10)) # probability

					# probabilities of sampling bases depending of base qualities
					if (geno=="A") probs <- c( 1-p, p/3, p/3, p/3);
					if (geno=="C") probs <- c( p/3, 1-p, p/3, p/3)

					# sample reads and concatenate
					alls=c(alls, sample(x=c(".","C","G","T"), size=1, prob=probs))

				} # end each read

			} else { # end if data
			
				bqs=c(bqs,"!")
				alls=c(alls,",")
				depth[n,i]=1

			} # end if not data

		} # end if haploid

		# diploid case
		if (ncopy[n]==2) {
	
			# priors: AA AC CC
			priors=c(pp^2,2*pp*qq,qq^2)

			genos=c("AA", "AC", "CC")
			geno=genos[sample(1:length(genos),1, prob=priors)]

			# daf
			if (geno=="AC") daf=daf+1
			if (geno=="CC") daf=daf+2

			if (depth[n,i]>0) { # if data

        	        	# sample qualities
                		ibq=round(rnorm(mean=mbqual,sd=2,n=depth[n,i])); bqs=c(); alls=c(); ibq[which(ibq<0)]=0
			
				# for each read
				for (j in 1:depth[n,i]) {
					bqs=c(bqs, substring(pscores, ibq[j]-1, ibq[j]-1))
					ps=ibq[j]; p=10^(-(ps/10)) # probability
					# probabilities of sampling bases depending of base qualities
					if (geno=="AA") probs=c( 1-p, p/3, p/3, p/3)
					if (geno=="AC") probs=c( ((1/2)*(1-p))+((1/2)*(p/3)), ((1/2)*(1-p))+((1/2)*(p/3)), ((0/2)*(1-p))+((2/2)*(p/3)), ((0/2)*(1-p))+((2/2)*(p/3)) )
					if (geno=="CC") probs=c( p/3, 1-p, p/3, p/3)
					alls=c(alls, sample(x=c(".","C","G","T"), size=1, prob=probs))
				} # end for each read

			} else { # end if data
				
				bqs=c(bqs,"!")
				alls=c(alls,",")
				depth[n,i]=1
							
			} # end if not data# end if data

		} # end if diploid

		# triploid case
		if (ncopy[n]==3) {
	
			# priors: AAA AAC ACC CCC
			priors=c(pp^3,3*pp^2*qq,3*pp*qq^2,qq^3)
			genos=c("AAA", "AAC", "ACC", "CCC")
			geno=genos[sample(1:length(genos),1, prob=priors)]

			# daf
			if (geno=="AAC") daf=daf+1
			if (geno=="ACC") daf=daf+2
			if (geno=="CCC") daf=daf+3

			if (depth[n,i]>0) { # if data

				# sample qualities
       		  		ibq=round(rnorm(mean=mbqual,sd=2,n=depth[n,i])); bqs=c(); alls=c(); ibq[which(ibq<0)]=0
		
				# for each read
				for (j in 1:depth[n,i]) {
					bqs=c(bqs, substring(pscores, ibq[j]-1, ibq[j]-1))
		        		ps=ibq[j]; p=10^(-(ps/10)) # probability
					# probabilities of sampling bases depending of base qualities
					if (geno=="AAA") probs=c( (1-p), p/3, p/3, p/3)
					if (geno=="CCC") probs=c( p/3, 1-p, p/3, p/3)
					if (geno=="AAC") probs=c( ((2/3)*(1-p))+((1/3)*(p/3)), ((1/3)*(1-p))+((2/3)*(p/3)), ((0/3)*(1-p))+((3/3)*(p/3)), ((0/3)*(1-p))+((3/3)*(p/3)) )
					if (geno=="ACC") probs=c( ((1/3)*(1-p))+((2/3)*(p/3)), ((2/3)*(1-p))+((1/3)*(p/3)), ((0/3)*(1-p))+((3/3)*(p/3)), ((0/3)*(1-p))+((3/3)*(p/3)) )
					alls=c(alls, sample(x=c(".","C","G","T"), size=1, prob=probs))
				}

			} else { # end if data
				
				bqs=c(bqs,"!")
				alls=c(alls,",")
				depth[n,i]=1
							
			} # end if not data# end if data

		} # end if triploid

		# tetraploid case
		if (ncopy[n]==4) {
	
			# priors: AAAA AAAC AACC ACCC CCCC
			priors=c(pp^4,4*pp^3*qq,6*pp^2*qq^2,4*pp*qq^3,qq^4)
			genos=c("AAAA", "AAAC", "AACC", "ACCC", "CCCC")
			geno=genos[sample(1:length(genos),1, prob=priors)]
		
			# daf
			if (geno=="AAAC") daf=daf+1
			if (geno=="AACC") daf=daf+2
			if (geno=="ACCC") daf=daf+3
			if (geno=="CCCC") daf=daf+4

			if (depth[n,i]>0) { # if data

				# sample qualities
				ibq=round(rnorm(mean=mbqual,sd=2,n=depth[n,i])); bqs=c(); alls=c(); ibq[which(ibq<0)]=0
			
				# for each read
				for (j in 1:depth[n,i]) {
					bqs=c(bqs, substring(pscores, ibq[j]-1, ibq[j]-1))
					ps=ibq[j]; p=10^(-(ps/10)) # probability
					# probabilities of sampling bases depending of base qualities
					if (geno=="AAAA") probs=c( (1-p), p/3, p/3, p/3)
					if (geno=="CCCC") probs=c( p/3, 1-p, p/3, p/3)
					if (geno=="AAAC") probs=c( ((3/4)*(1-p))+((1/4)*(p/3)), ((1/4)*(1-p))+((3/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)) )
					if (geno=="AACC") probs=c( ((2/4)*(1-p))+((2/4)*(p/3)), ((2/4)*(1-p))+((2/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)))
					if (geno=="ACCC") probs=c( ((1/4)*(1-p))+((3/4)*(p/3)), ((3/4)*(1-p))+((1/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)), ((0/4)*(1-p))+((4/4)*(p/3)))
					alls=c(alls, sample(x=c(".","C","G","T"), size=1, prob=probs))
				}

			} else { # end if data
				
				bqs=c(bqs,"!")
				alls=c(alls,",")
				depth[n,i]=1
							
			} # end if not data# end if data

		} # end if tetraploid

		# pentaploid case
		if (ncopy[n]==5) {

			# priors: AAAAA AAAAC AAACC AACCC ACCCC CCCCC
			priors=c(pp^5,5*pp^4*qq,10*pp^3*qq^2,10*pp^2*qq^3,5*pp*qq^4,qq^5)
			genos=c("AAAAA", "AAAAC", "AAACC", "AACCC", "ACCCC", "CCCCC")
			geno=genos[sample(1:length(genos),1, prob=priors)]
                
			# daf
			if (geno=="AAAAC") daf=daf+1
			if (geno=="AAACC") daf=daf+2
			if (geno=="AACCC") daf=daf+3
			if (geno=="ACCCC") daf=daf+4
			if (geno=="CCCCC") daf=daf+5

			if (depth[n,i]>0) { # if data

				# sample qualities
                		ibq=round(rnorm(mean=mbqual,sd=2,n=depth[n,i])); bqs=c(); alls=c(); ibq[which(ibq<0)]=0
		
				# for each read
				for (j in 1:depth[n,i]) {
					bqs=c(bqs, substring(pscores, ibq[j]-1, ibq[j]-1))
			       	 	ps=ibq[j]; p=10^(-(ps/10)) # probability
					# probabilities of sampling bases depending of base qualities
					if (geno=="AAAAA") probs=c( (1-p), p/3, p/3, p/3)
					if (geno=="CCCCC") probs=c( p/3, 1-p, p/3, p/3)
					if (geno=="AAAAC") probs=c( ((4/5)*(1-p))+((1/5)*(p/3)), ((1/5)*(1-p))+((4/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)) )
					if (geno=="AAACC") probs=c( ((3/5)*(1-p))+((2/5)*(p/3)), ((2/5)*(1-p))+((3/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)))
					if (geno=="AACCC") probs=c( ((2/5)*(1-p))+((3/5)*(p/3)), ((3/5)*(1-p))+((2/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)))
					if (geno=="ACCCC") probs=c( ((1/5)*(1-p))+((4/5)*(p/3)), ((4/5)*(1-p))+((1/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)), ((0/5)*(1-p))+((5/5)*(p/3)))
					alls=c(alls, sample(x=c(".","C","G","T"), size=1, prob=probs))
				}

			} else { # end if data
				
				bqs=c(bqs,"!")
				alls=c(alls,",")
				depth[n,i]=1
							
			} # end if not data# end if data

		} # end if pentaploid

		# line for mpileup
		if (opt$pool==FALSE) {
			linea=c(linea, depth[n,i], paste(alls, sep="", collapse=""), paste(bqs,sep="",collapse=""))
		} else {
			pool_alls=c(pool_alls, alls)
			pool_bqs=c(pool_bqs, bqs)
		}

		# line for real data
		linea_real=c(linea_real, geno)

	} # end for samples	
	
	if (opt$pool==FALSE) {
		cat(linea, sep="\t")
		cat("\n")
	} else {
		linea=c(linea, sum(depth[,i]), paste(pool_alls, sep="", collapse=""), paste(pool_bqs,sep="",collapse=""))
		if (sum(depth[,i])!=length(pool_alls)) cat("ERROR!!!, depth is different than nr of reads!!!! QUIT!")
		cat(linea, sep="\t")
		cat("\n")
	}

	# write real	
	if (!is.null(opt$out)) {
		# daf
		nchroms=sum(as.numeric(ncopy))
	        linea_real=c(linea_real, (daf/nchroms))
		# write
		cat(linea_real, sep="\t", file=fout_real, append=T)
		cat("\n", file=fout_real, append=T)
	}

} # end for sites






