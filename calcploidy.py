import math

def combinations_with_rep(iterable, r): # combinations with replacements, edited function from itertools
    pool = list(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield list(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield list(pool[i] for i in indices)

def calcGenoLogLikeN_MajorMinor(N,read,site,major,minor):
    alleles = ['A','C','G','T']
    log_likes=[0.0]*(ploidy+1)
    it = -1
    phredScale=33
    mm = [major,minor] 
    mmList = list(combinations_with_rep(mm,ploidy) # List of major minor combinations
    nList = list(range(1,ploidy+1))
    jList = ["j"+ str(i) for i in nList]
    # cycle across all possible genotypes
    readLen = len(read.base)
    for subList in mmList:
        it += 1
        for item in subList:
            for i in range(readLen):
                bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
                sublike = 0.0
                if alleles[item] == read.base[i]:
                    sublike += (1-bP)/ploidy
                else:
                    sublike += (bP/3)/ploidy
                log_likes[it] += math.log(sublike)
    return(log_likes)


