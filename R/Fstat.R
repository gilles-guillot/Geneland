#' Fstat
#' @description Computes F statistics according to Weir and Cockerham's 
#'   estimators.
#'   Missing values are allowed and accounting for in the computation of
#'   Fst.
#'   The presence of missing values involves a downward bias in the
#'   computation of Fis.
#'   This function should not be used on haploid data.
#'    @param genotypes Diploid codominant genotype data. A matrix with one line per
#'     individual and 2 columns per locus
#'   @param npop total number of popluation present in the dataset
#'   @param pop.mbrship Vector of integers giving the population membership
#'     for each individual
#'    @param ploidy Integer: 1 or 2 (default is 2) under development. Do
#'      not use for haploid data.
#' 
#' @return A list with components 
#'   Fis  (vector of estimations of within-population Fis)
#'   and Fst  (matrix of estimations the pairwise population
#'     Fst's) 
#' @export 
Fstat <- function(genotypes,npop,pop.mbrship,ploidy=2)
  { 
    #print("debut fonction R Fstat ")
    #print(genotypes)
    #print(allele.numbers)
    #print(npop)
    #print(pop.mbrship)

    
##     if(ploidy == 1)
##       {
##         data.tmp <- matrix(nrow=nrow(genotypes),
##                            ncol=ncol(genotypes)*2)
##         data.tmp[,seq(1,ncol(genotypes)*2-1,2)] <- genotypes
##         data.tmp[,seq(2,ncol(genotypes)*2,2)] <- genotypes
##         genotypes <- data.tmp
##       }
     if(ploidy == 1) stop("Fstat not implemented for haploid data")
     
    format <- FormatGenotypes(genotypes,ploidy=ploidy)
    genotypes <- format$genotypes
    allele.numbers <- format$allele.numbers
    
    if(sum(is.na(genotypes)) != 0)
      {
        warning("Genotypes contain missing values which might bias computations")
        genotypes[is.na(genotypes)] <- -999
      }
                                        # Pairwise computations
    Fis=rep(-999,npop)
    Fst=matrix(nrow=npop,ncol=npop,-999)

    if(npop > 1)
      {
        for(iclass1 in 1:(npop-1))
          {
            for(iclass2 in (iclass1+1):npop)
              {
                sub1 <- pop.mbrship==iclass1
                sub2 <- pop.mbrship==iclass2
                if((sum(sub1)!=0)  & (sum(sub2)!=0))
                  {
                    ztmp <- genotypes[sub1 | sub2,]
                    nindivtmp <- nrow(ztmp)
                    pop.mbrshiptmp <- pop.mbrship[sub1 | sub2]
                    pop.mbrshiptmp[pop.mbrshiptmp==iclass1] <- 1
                    pop.mbrshiptmp[pop.mbrshiptmp==iclass2] <- 2
                    tabindiv <- matrix(nrow=nindivtmp,ncol=2,data=-999)
                    kk <- numeric(2)
                    effcl <- table(pop.mbrshiptmp)
                    nloc <- length(allele.numbers)
                    nloc2 <- 2*nloc
                    Fistmp <- Fsttmp <- Fittmp <- -999
                    ## print("Computing Fst")
                    res<- .Fortran("ggfst",
                                   PACKAGE="Geneland",
                                   as.integer(nindivtmp),
                                   as.integer(nloc),
                                   as.integer(nloc2),
                                   as.integer(allele.numbers),
                                   as.integer(2),
                                   as.integer(effcl),
                                   as.integer(ztmp),
                                   as.integer(pop.mbrshiptmp),
                                   as.integer(tabindiv),
                                   as.integer(kk),
                                   as.double(Fistmp),
                                   as.double(Fsttmp),
                                   as.double(Fittmp))
                    Fst[iclass1,iclass2] <- res[[12]][1]
                  }
              }
          }
      }
    for(iclass1 in 1:npop)
      {
        sub1 <- pop.mbrship==iclass1
        if((sum(sub1)!=0))
          {

            ztmp <-  rbind(genotypes[sub1,],genotypes[sub1,])
            nindivtmp <- nrow(ztmp)
            pop.mbrshiptmp <- c(rep(1,sum(sub1)),rep(2,sum(sub1)))
            tabindiv <- matrix(nrow=nindivtmp,ncol=2,data=-999)
            kk <- numeric(2)
            effcl <- table(pop.mbrshiptmp)
            nloc <- length(allele.numbers)
            nloc2 <- 2*nloc
            Fistmp <- Fsttmp <- Fittmp <- -999
            ## print("Computing Fis")
            res<- .Fortran("ggfst",
                           PACKAGE="Geneland",
                           as.integer(nindivtmp),
                           as.integer(nloc),
                           as.integer(nloc2),
                           as.integer(allele.numbers),
                           as.integer(2),
                           as.integer(effcl),
                           as.integer(ztmp),
                           as.integer(pop.mbrshiptmp),
                           as.integer(tabindiv),
                           as.integer(kk),
                           as.double(Fistmp),
                           as.double(Fsttmp),
                           as.double(Fittmp))
            Fis[iclass1] <- res[[11]][1]
          }
      }
    Fst[lower.tri(Fst,diag=TRUE)] <- 0
    Fst <- Fst + t(Fst)
    Fis[Fis==-999] <- NA
    Fst[Fst==-999] <- NA     
    list(Fis=Fis,
         Fst=Fst)
  }
 



##################
##################

