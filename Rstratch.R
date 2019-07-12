pascalGet <- function(h) {
  for(i in 0:(h-1)) {
    s <- 0
    for(k in 0:(h-i)) s <- paste(s, "  ", sep="")
    for(j in 0:i) {
      s <- paste(s, sprintf("%3d ", choose(i, j)), sep="")
    }
    print(s)
  }
}


binomialCoefficents <- function(ploidy) {
  coef = c()   # Initial vector of binomial coefficents
  for(i in 0:ploidy) { coef = c(coef, choose(ploidy, i)) }
  return(coef)
}

binomialCoefficents(3)






pPowers <- function(ploidy) {
  ppow = c()   # Initial vector of pp powers
  for (i in ploidy:0) { ppow = c(ppow, i) }
  qpow = rev(ppow)
  return(pow)
}

qPowers <- function(ploidy) {
  qpow = rev(ppow)
  return(pow)
}



pqPowers(2)

hwe <- function(ploidy) { # Hardy Weinberg Equilibrium genotype probs
  pq = c(pp, qq)
  pqPowers(ploidy)
  for (i in 1:ploidy) {
  priors <- c(pq[1] ,pq[2])
  }
}

hwe(ploidy)
priors

ploidy = 1

pp = 1
qq = 2

pq <- c(pp,qq)

a = character(0)
x = character(length=ploidy)

for (i1 in 2:(ploidy+1)){
  for (i2 in pq){
    a = paste(a, i2)
  }
} 

a
x

eval(parse(text = paste(paste("priors <- c(", paste(pq, collapse=",") ),")")))

a[1][1]
pp = 1
qq = 0

hwe(3)

cList
a <- c(1,2,3)
rev(a)

for (i in 1:10) { x <- c(0, x) + c(x, 0); print(x) }

choose(3,3)

x = 0.4
y= 0.6
polym(c(x,y), degree = 2, raw = TRUE)

# R package 'mpoly' required
mp()
mpoly()


