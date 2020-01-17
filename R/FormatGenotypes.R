#' FormatGenotypes

#' @description Takes genotype data as a matrix with
#'   one line per individual and two columns per locus,
#'   with alleles coded by integers (number of replications for
#'                                   micro-satellites data).
#'   Build a new matrix with alleles codes as consecutive integers.
#'   If a locus  has 7 alleles they will be coded as 1,2,...7.
#'   Since version 1.0.1, this function does not have to be called by
#'   users. It is called through MCMC.
#'  @param genotypes A matrix with
#'     one line per individual and two columns per locus,
#'     with alleles coded by integers
#'  @param ploidy 1 or 2
#' @return 
#'   A list with elements: 
#'     genotypes (matrix with one line per individual and
#'       one or two columns per locus with alleles coded by integers) and 
#'  allele.numbers (a vector giving the number of possible alleles per locus
#'  @export

FormatGenotypes <- function(genotypes,ploidy)
  {
    formatted <- genotypes
    ## vector whose entries
    ## will be the numbers
    ## of allele per locus
    nall <- numeric(ncol(genotypes)/2)

    if(ploidy==2)
      {
        for(ic in 1:(ncol(genotypes)/2))
          {
            ens <- sort(unique(c(genotypes[,2*ic-1],genotypes[,2*ic])))
            max.ens <- length(ens)
            for(il in 1:(dim(genotypes)[1]))
              {
                formatted[il,2*ic-1] <- ifelse(is.na(genotypes[il,2*ic-1]),
                                               NA,
                                               (1:max.ens)[genotypes[il,2*ic-1]==ens])
                
                formatted[il,2*ic]   <- ifelse(is.na(genotypes[il,2*ic]),
                                               NA,
                                               (1:max.ens)[genotypes[il,2*ic]==ens])
              }
            nall[ic] <- max.ens
          }
      }
    if(ploidy==1)
      {
         for(ic in 1:(ncol(genotypes)))
           {
             ens <- sort(unique(genotypes[,ic]))
             max.ens <- length(ens)
             for(il in 1:(dim(genotypes)[1]))
               {
                 formatted[il,ic] <- ifelse(is.na(genotypes[il,ic]),
                                                NA,
                                                (1:max.ens)[genotypes[il,ic]==ens])
               }
             nall[ic] <- max.ens
           }
      }
    formatted <- as.matrix(formatted)
    ## keep NA code as NA at this step
    list(genotypes=formatted,allele.numbers=nall)
  }
