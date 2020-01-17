#
# Spatialised populations via Poisson-Voronoi tiling
# + genotypes according to F-model or D-model
#
#' @title    Simulation of multi-locus genetic data from the spatial F-model
#'
#' @param nindiv Simulation of multi-locus genetic data from the spatial F-model 
#' @param coordinates Matrix (2 rows, nindiv columns) of spatial coordinates of individuals
#' @param coord.lim Vector of limits of spatial domain to be considered (x min, x max, y min, y max)
#' @param number.nuclei Integer: number of nuclei in the Voronoi tessellation
#' @param coord.nuclei Coordinates of nuclei of Voronoi tessellation
#' @param color.nuclei Population labels of the nuclei (vector of integers of size number.nuclei)
#' @param nall Vector of integers giving number of alleles at each locus
#' @param npop Number of populations
#' @param freq.model model for frequencies:"Correlated" or "Uncorrelated"
#' @param drift Vector (of size \code{npop}) of drift factors between 0 and 1 (only for the Correlated model)
#' @param dominance A character string "Codominant" or "Dominant". If
#' "Dominant" is chosen, the first allele is treated as a recessive
#' allele and all heterozigous are converted into homozigous for the
#' dominant allele. The presence of the "dominant" allele is coded as 1,
#' the absence of the "dominant" allele is coded as 0.
#' @param plots Logical: if TRUE, spatial coordinates are ploted
#' @param ploth Logical: if TRUE, barplots for allele frequencies are ploted
#'
#' @return A list of variables involved in the simulation. The elements of
#' this list are: 
#'   coordinates, genotypes, allele.numbers, number.nuclei, coord.nuclei, color.nuclei, 
#' frequencies,  ancestral.frequencies, drifts, index.nearest.nucleus
#' @export
simFmodel <- function(nindiv,
                      coordinates,
                      coord.lim=c(0,1,0,1),
                      number.nuclei,
                      coord.nuclei,
                      color.nuclei,
                      nall,
                      npop,
                      freq.model="Uncorrelated",
                      drift,
                      dominance="Codominant",
                      plots=FALSE,
                      ploth=FALSE)
  {
    if(dominance=="Dominant" & sum(nall !=rep(2,length(nall)))>0)
      {stop("Dominant option only for bi-allelic loci")}
                                        # spatial coord of indviduals
    if(missing(coordinates))
      {
        coordinates <- rbind(runif(min=coord.lim[1],max=coord.lim[2],nindiv),
                   runif(min=coord.lim[3],max=coord.lim[4],nindiv))
      }
    
                                        # centers and colors of tiles
    if(missing(coord.nuclei))
      {
        coord.nuclei <-  rbind(runif(min=coord.lim[1],max=coord.lim[2],number.nuclei),
                               runif(min=coord.lim[3],max=coord.lim[4],number.nuclei))

        color.nuclei <- sample(x=1:npop,size=number.nuclei,replace=TRUE)
        
      } 
    ppvois <- numeric(nindiv)
    for(iindiv in 1:nindiv)
      {
        k <- 1
        dd <- 1e+9
        for(ipp in 1:number.nuclei)
          {
            ddnew <- (coord.nuclei[1,ipp]-coordinates[1,iindiv])^2+
              (coord.nuclei[2,ipp]-coordinates[2,iindiv])^2
            if(ddnew < dd)
              {
                k <- ipp
                dd <- ddnew
              }
          } 
        ppvois[iindiv] <- k
      }

    nloc <- length(nall)


    ## drift and freq in ancestral pop
    if(freq.model == "Correlated")
      {
        # drift parameters
        if(missing(drift))
          {
            drift <- rbeta(shape1=2,shape2=20,n=npop)
          }
        # alleles frequencies in ancestral population
        fa <- matrix(nrow=nloc,ncol=max(nall),data=-999)
        for(iloc in 1:nloc)
          {
            fa[iloc,1:nall[iloc]] <- rexp(n=nall[iloc])
            fa[iloc,1:nall[iloc]] <- fa[iloc,1:nall[iloc]] /
              sum(fa[iloc,1:nall[iloc]])
          }
      }
    
    ## alleles frequencies in present time population
    freq <- array(dim=c(npop,nloc,max(nall)),data=-999)
    for(iclass in 1:npop)
      {
        for(iloc in 1:nloc)
          {
            if(freq.model == "Uncorrelated")
              {
                freq[iclass,iloc,1:nall[iloc]] <- rgamma(n=nall[iloc],
                                                         scale=1,
                                                         shape=1)
                freq[iclass,iloc,1:nall[iloc]] <- freq[iclass,iloc,1:nall[iloc]] /
                  sum(freq[iclass,iloc,1:nall[iloc]])
              }
            if(freq.model == "Correlated")
              {
                freq[iclass,iloc,1:nall[iloc]] <- rgamma(n=nall[iloc],
                                                         scale=1,
                                                         shape=fa[iloc,(1:nall[iloc])]*
                                                         (1-drift[iclass])/drift[iclass])
                freq[iclass,iloc,1:nall[iloc]] <- freq[iclass,iloc,1:nall[iloc]] /
                  sum(freq[iclass,iloc,1:nall[iloc]])
              }
          }
      }
    
    ##                                     # impose condition on  freqs :  f111 > f211
    ## if(npop > 1)
    ##   {
    ##     if(freq[1,1,1] < freq[2,1,1])
    ##       {
    ##         tmp <- freq[2,,]
    ##         freq[2,,] <- freq[1,,]
    ##         freq[1,,] <- tmp
    ##       }
    ##   }
          
    
    ## genotypes     
    z <- matrix(nrow=nindiv,ncol=nloc*2)
    for(iclass in 1:npop)
      {
                                        #for(iindiv in (1:nindiv)[color.nuclei[ppvois]==iclass] )
                                        #{
        subclass <- (1:nindiv)[color.nuclei[ppvois]==iclass]
        for(iloc in 1:nloc)
          {
            z[subclass,c(2*iloc-1,2*iloc)] <- sample(x=1:nall[iloc],
                                                     size=2*length(subclass),
                                                     prob=freq[iclass,iloc,1:nall[iloc]],
                                                     replace=TRUE)
                #z[iindiv,2*iloc-1] <- rdiscr(freq[iclass,iloc,1:nall[iloc]])
                #z[iindiv,2*iloc]   <- rdiscr(freq[iclass,iloc,1:nall[iloc]])
                                        #              }
          }
      }
    
    if(plots==TRUE)
      { 
        dev.new()
        plot(coordinates[1,],coordinates[2,],
             xlab="x coordinates",ylab="y coordinates")
        points(coord.nuclei[1,],coord.nuclei[2,],col=color.nuclei,cex=2,pch=16)
        text(coord.nuclei[1,],coord.nuclei[2,],1:number.nuclei,pos=1,col=2,pch=2,cex=1.2)
        text(coordinates[1,],coordinates[2,],ppvois,pos=1)
      }
    if(ploth==TRUE)
      {
        if(freq.model=="Correlated")
          {
           dev.new()
            par(mfrow=c(floor(sqrt(nloc)+1),
                  floor(sqrt(nloc))))
                                        #par(mfrow=c(1,nloc))
            for(iloc in 1:nloc)
              {
                plot(1:nall[iloc],fa[iloc,1:nall[iloc]],
                     type="h",col=2,xlab="",axes=FALSE,
                     sub=paste("Locus",iloc),ylim=c(0,1),
                     main="Frequencies in ancestral population")
              }
          }
        
       dev.new()
        #par(mfrow=c(npop,nloc))
        for(iclass in 1:npop)
          {
           dev.new()
            par(mfrow=c(floor(sqrt(nloc)+1),
                  floor(sqrt(nloc))))
            #par(mfrow=c(1,nloc))
            for(iloc in 1:nloc)
              {
                hist(c(z[color.nuclei[ppvois]==iclass,2*(iloc-1)+1],
                       z[color.nuclei[ppvois]==iclass,2*(iloc-1)+2]),
                     breaks=seq(.5,nall[iloc]+.5,1),
                     prob=TRUE,main=paste("Histogram of freq. in pop.",iclass,
                              ", locus",iloc),
                     xlab="",ylim=c(0,1),axes=FALSE)
                points(1:nall[iloc],freq[iclass,iloc,1:nall[iloc]],
                       type="h",col=2)
              }
          }
      }
    
    true.codom.genotypes <- NULL
    if(dominance == "Dominant")
      {
        true.codom.genotypes <- z
        for(iloc in 1:nloc)
          {
            for(iindiv in 1:nindiv)
              {
                if(sum(z[iindiv,c(2*iloc-1,2*iloc)]) ==3) z[iindiv,c(2*iloc-1,2*iloc)] <- c(2,2)
              }
          }
        z <- z[,seq(1,2*nloc-1,2)]
        z <- z-1
      }
    
 ##    if(freq.model=="Uncorrelated")
##       {
##         res <- list(coordinates=t(coordinates),
##                     genotypes=z,
##                     allele.numbers=nall,
##                     number.nuclei=number.nuclei,
##                     coord.nuclei=t(coord.nuclei),
##                     color.nuclei=color.nuclei,
##                     frequencies=freq,
##                     index.nearest.nucleus=ppvois)
##         return(res)
##       }
##     if(freq.model=="Correlated")
##       {
##         res <- list(coordinates=t(coordinates),
##                     genotypes=z,
##                     true.codom.genotypes=true.codom.genotypes,
##                     allele.numbers=nall,
##                     number.nuclei=number.nuclei,
##                     coord.nuclei=t(coord.nuclei),
##                     color.nuclei=color.nuclei,
##                     frequencies=freq,
##                     ancestral.frequencies=fa,
##                     drifts=drift,
##                     index.nearest.nucleus=ppvois)
##         return(res)
##       }

    res <- list(coordinates=t(coordinates),
                genotypes=z,
                allele.numbers=nall,
                number.nuclei=number.nuclei,
                coord.nuclei=t(coord.nuclei),
                color.nuclei=color.nuclei,
                frequencies=freq,
                index.nearest.nucleus=ppvois,
                dominance=dominance)
    

    if(freq.model=="Correlated")
      {
        res <- c(res,
                 list(ancestral.frequencies=fa,
                      drifts=drift))
      }
    if(dominance=="Dominant")
      {
        res <- c(res,
                 list(true.codom.genotypes=true.codom.genotypes))
      }

    return(res)
}

 
#############
##
## debugging
##



# simdata <- simFmodel(nindiv=100,
#                           coordinates = NULL,
#                           coord.lim=c(0,1,0,1),
#                           number.nuclei=1,
#                           #coord.nuclei = u,
#                           #color.nuclei = c,
#                           nall=rep(2,2),
#                           npop=1,
#                           #freq.model="Correlated",
#                           drift=rbeta(shape1=2,shape2=20,K[iset]),
#                           #seed=iset,
#                           plots =FALSE,
#                           ploth = TRUE)
    
