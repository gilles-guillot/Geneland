#' @title nullify
#' @description Simulates null alleles
#' @usage nullify(genotypes, nall.null = 1, nloc.null)
#'  @param genotypes  a matrix of genotypes as produced by simFmodel and simIBD
#'  @param nall.null number of null alleles on each locus
#'  @param nloc.null number of loci carrying null alleles
#' @return A list with component: 
#'    @param genotypes  The new genotypes after alteration.
#'     @export
nullify <- function(genotypes,nall.null=1,nloc.null)
  {
    ## simulate a null allele for each allele in nall.null
    ## return the empirical frequency of this allele
    z <- genotypes
    for(iloc in 1:nloc.null)
      {
        ## ina <- sample(1:length(unique(as.vector(z[,c(2*iloc-1,2*iloc)]))),nall.null)
        for(iall in nall.null)
          {
            zz <- z[,c(2*iloc-1,2*iloc)]
            #print(iall)
            for(iindiv in 1:nrow(genotypes))
              {
                if(sum(is.na(z[iindiv,c(2*iloc-1,2*iloc)]))==0)
                  {
                    if(sum(z[iindiv,c(2*iloc-1,2*iloc)]==iall)==2)
                      {
                        zz[iindiv,] <- c(NA,NA)
                      }else{
                        if(z[iindiv,2*iloc-1]==iall) zz[iindiv,1] <- z[iindiv,2*iloc]
                        if(z[iindiv,2*iloc]==iall) zz[iindiv,2] <- z[iindiv,2*iloc-1]
                      }
                  }else{
                    if(sum(is.na(z[iindiv,c(2*iloc-1,2*iloc)]))==2)
                      {
                        zz[iindiv,] <- z[iindiv,c(2*iloc-1,2*iloc)]
                      }else{
                        if(is.na(z[iindiv,2*iloc]))
                          {
                            if(z[iindiv,2*iloc-1]==iall) zz[iindiv,1] <- NA
                          }
                        if(is.na(z[iindiv,2*iloc-1]))
                          {                         
                            if(z[iindiv,2*iloc]==iall) zz[iindiv,2] <- NA
                          }
                      }
                  }
              }
                                        # print(zz)
            z[,c(2*iloc-1,2*iloc)] <- zz
                                        # print(z)
          }
      }
    ## list(genotypes=z)
    return(z)
  }









