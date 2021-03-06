#' 
#' @title Spatialised populations via  colored Poisson-Voronoi tiling
#' genotypes according to RFDM
#' quantitative variable: Normal,Gamma prior Normal likelihood
#' @param  nindiv Number of indivuals 
#' 
#' 
#' @param  coord.indiv Coordinates of the individuals 
#' 
#' @param  coord.lim Limits of the geographical domain. The domain is
#'   supposed to be rectangular and the limits are given as  (abs min, abs
#'                                                            max, ord min, ord max) 
#' 
#' @param  npop Number of Populations 
#' 
#' @param  rate Rate of the Poisson process governing the hidden tessellation 
#' 
#' @param  number.nuclei Number of nuclei in the tessellation (if given, 
#'                                                            then \code{rate}  is ignored)
#' 
#' @param  coord.nuclei Coordinates of the nuclei (the number of
#'                                                coordinates of the nuclei
#'                                                given here as a matrix has to comply with \code{number.nuclei} ) 
#' 
#' @param  color.nuclei Population membeship of the nuclei: a vector of
#'   integer of length \code{number of nuclei} whose values are between 1
#'   and \code{npop} 
#' 
#' @param  sim.gen Logical to say whether genetic data should be simulated  
#' 
#' @param  allele.numbers A vector giving the number of alleles observed at
#'   each locus 
#' 
#' @param  IBD Logical. If TRUE, then the allele frequencies are simulated
#'   according to an IBD model. If FALSE, panmixia is assumed. 
#' 
#' @param  model Model of spatial covariance function used for the
#'   underlying Gaussian fields (see documentation of package
#'                               \code{RandomFields}
#'                               for details) 
#' 
#' @param  alpha Parameter of the spatial Dirichlet vector field
#'   of frequencies (a positive real)  
#' 
#' @param  beta Scale parameter of the spatial covariance function used for the
#'   underlying Gaussian fields. A positive real number (see documentation of package
#'                                                       \code{RandomFields}
#'                                                       for details) 
#' 
#' @param  gamma Smoothing parameter of spatial covariance function used for the
#'   underlying Gaussian fields. (see documentation of package
#'                                \code{RandomFields}
#'                                for details) 
#' 
#' @param  sim.quanti Logical to say whether quantitative data should be
#'   simulated  
#' 
#' @param  nquanti.var Number of quantitative variables to be simulated 
#' 
#' 
#' @param  mean.quanti Mean of the quantitative variables in the various
#'   groups. A matrix with \code{npop} lines and \code{nvar.quant} columns  
#' 
#' @param  sd.quanti Standard deviation  of the quantitative variables in
#'   the various groups.
#'   A matrix with \code{npop} lines and \code{nvar.quant} columns 
#' 
#' 
#' 
#' @param  seed.coord Random seed to initialise the simulation of the
#'   coordinates
#'   (mostly for debugging) 
#' 
#' @param  seed.tess Random seed to initialise the simulation of the
#'   tessellation  (mostly for debugging) 
#' 
#' @param  seed.freq  Random seed to initialise the simulation of the
#'   frequencies  (mostly for debugging) 
#' 
#' @param  give.freq.grid Logical to tell whether frequencies  on a grid are
#'   also returned 
#' 
#' @param  give.tess.grid Logical to tell whether population memberships of pixels  on a grid are also returned 
#' 
#' @param  npix A vector of two integers telling how many horizontal and
#'   vertical pixel should contain the grid for the graphical
#'   representations 
#' 
#' @param  comp.Fst Logical to tell whether Fst, Fis and Fit should be
#'   computed 
#' 
#' @param  comp.Dsigma2 Logical to tell whether IBD index Dsgma2 should be
#'   computed 
#' 
#' @param  comp.diff Logical to tell whether the local differentiation across the barriers
#'   should be computed 
#' 
#' @param  width Real number specifying the width around the barrier in
#'   the computation of its local differentiation 
#' 
#' @param  plot.pairs.borders Logical to tell whether the pairs of
#'   individuals coming into the computation of the differentiation of the barriers
#'   should be plotted 
#' @export
#'
#' @examples \dontrun{
#' dataset <- simdata(nindiv=100,
#'                    sim.gen=TRUE,
#'                    number.nuclei=10,
#'                    allele.numbers=rep(5,3),
#'                    model="stable",
#'                    IBD=TRUE,
#'                    alpha=1,
#'                    beta=1,
#'                    gamma=1,
#'                    npop=3,
#'                    sim.gen=TRUE,
#'                    give.tess.grid=TRUE,
#'                    give.freq.grid=TRUE,
#'                    npix=c(10,10),
#'                    comp.Fst=TRUE,
#'                    comp.Dsigma2=TRUE,
#'                    comp.diff=TRUE,
#'                    width=0.1,
#'                    plot.pairs.borders=TRUE)
#' }
#' @export
simdata <- function(nindiv,
                    coord.indiv,
                    coord.lim=c(0,1,0,1),
                    npop,
                    rate,                        
                    number.nuclei,
                    coord.nuclei,
                    color.nuclei,
                    allele.numbers,
                    sim.gen=FALSE,
                    IBD=TRUE,
                    model="stable",
                    alpha=1,
                    beta=1,
                    gamma=1.8,
                    sim.quanti=FALSE,
                    nquanti.var,
                    mean.quanti,
                    sd.quanti,
                    seed.coord,
                    seed.tess, 
                    seed.freq,
                    give.tess.grid=FALSE,
                    give.freq.grid=FALSE,
                    npix,
                    comp.Fst=FALSE,
                    comp.Dsigma2=FALSE,
                    comp.diff=FALSE,
                    width,
                    plot.pairs.borders=FALSE)

  {
    if(missing(rate)) rate <- nindiv/2
    if(missing(number.nuclei)) number.nuclei <- rpois(1,rate)

    
    ### spatial coord of indviduals
    if(missing(coord.indiv))
      {
        if(!missing(seed.coord)) set.seed(seed.coord)
        coord.indiv <- rbind(runif(min=coord.lim[1],max=coord.lim[2],nindiv),
                             runif(min=coord.lim[3],max=coord.lim[4],nindiv))
      }
    if(give.tess.grid  | give.freq.grid)
      {
        coord.grid <- matrix(nrow=2,ncol=npix[1]*npix[2],NA)
        coord.grid[1,] <- rep(seq(from=coord.lim[1],to=coord.lim[2],
                                  length=npix[1]),
                              npix[2])
        coord.grid[2,] <- as.vector(matrix(nrow=npix[1],ncol=npix[2],byrow=TRUE,
                                           rep(seq(from=coord.lim[3],to=coord.lim[4],
                                                   length=npix[2]),
                                                   npix[1])))
      }
    if(give.freq.grid){coord.all <- cbind(coord.indiv,coord.grid)}else
    {coord.all <- coord.indiv}
      
    
    ## centers and colors of tiles
    if(!missing(seed.tess)) set.seed(seed.tess)
    if(missing(coord.nuclei))
      {
        coord.nuclei <-  rbind(runif(min=coord.lim[1],max=coord.lim[2],number.nuclei),
                               runif(min=coord.lim[3],max=coord.lim[4],number.nuclei))
      }
    if(missing(color.nuclei))
      {
        color.nuclei <- sample(1:npop,number.nuclei,replace=TRUE)
      }


    if(give.tess.grid  | give.freq.grid)
      {
        coord.all.tmp <- cbind(coord.indiv,coord.grid)
      }else

    {
       coord.all.tmp <-coord.indiv
     }
    nearest.nucleus <- numeric(ncol(coord.all.tmp))
    for(iindiv in 1:ncol(coord.all.tmp))
      {
        k <- 1
        dd <- 1e+9
        for(ipp in 1:number.nuclei)
          {
            ddnew <- (coord.nuclei[1,ipp]-coord.all.tmp[1,iindiv])^2+
              (coord.nuclei[2,ipp]-coord.all.tmp[2,iindiv])^2
            if(ddnew < dd)
              {
                k <- ipp
                dd <- ddnew
              }
          }
        nearest.nucleus[iindiv] <- k
      }
  
    nearest.nucleus.indiv <- nearest.nucleus[1:nindiv]
    if(give.freq.grid | give.tess.grid)
      {
        nearest.nucleus.grid <- nearest.nucleus[-(1:nindiv)]
      }

   ## simulating genotypes 
   if(sim.gen)
     {
       nloc <- length(allele.numbers)
       
       ## Gaussian random fields
       ## param=c(NA,variance,nugget,scale,...)
       if(!missing(seed.freq)) set.seed(seed.freq)
       e <- array(NA,dim=c(npop,ncol(coord.all),nloc,max(allele.numbers)))
       ## Gaussian field stored in an
       ## intermediate array to comply with
       ## RandomFileds package format
       if(IBD)
         {
           print("Calling GaussRF")
           ff <- GaussRF(x=coord.all[1,],
                         y=coord.all[2,],
                         grid=FALSE,
                         model=model,
                         param=c(0,1,0,beta,gamma),##param=c(mean, variance, nugget, scale, ...). Here ‘...’ stands for additional parameters such as ν in the whittle model.
                         n=sum(allele.numbers)*npop)
           print("End of GaussRF")
         }else
       {
         ff <- matrix(nrow=ncol(coord.all),
                      ncol=sum(allele.numbers)*npop,
                      data=rnorm(sum(allele.numbers)*npop),
                      byrow=TRUE)
       }
       
       ## arrange in a better format 
       count <- 1
       for(ipop in 1:npop)
         {
           for(iloc in 1:nloc)
             {
               for(iall in 1:allele.numbers[iloc])
                 {
                   e[ipop,,iloc,iall] <- ff[,count]#10*(iall-1.5)*(1.5-ipop)#ff[,count]#
                   count <- count + 1
                 }
             }
         }
       
       
       ## tranformation into exponential-like
       ## values
       ## if(length(alpha)==1) CORRECTED ON 17/06/2008
       if(alpha==1)
         {
           y <- qexp(pnorm(e),rate=1)
         }
       else
         {
           y <- qgamma(pnorm(e),shape=alpha)
         }

    
       ## allele frequencies
       freq <- array(NA,dim=c(npop,ncol(coord.all),nloc,max(allele.numbers)))
       for(ipop in 1:npop)
         {
           for(iloc in 1:nloc)
             {
               for(iall in 1:allele.numbers[iloc])
                 {
                   freq[ipop,,iloc,iall] <- y[ipop,,iloc,iall]/
                     apply(y[ipop,,iloc,1:allele.numbers[iloc]],1,sum)
                 }
             }
         }
       print("freq. stored")
                    
       ## genotypes     
       z <- matrix(nrow=nindiv,ncol=nloc*2)
       for(iindiv in 1:nindiv)
         {
           ipop <- color.nuclei[nearest.nucleus[iindiv]]
           for(iloc in 1:nloc)
             {
               ff <- freq[ipop,iindiv,iloc,1:allele.numbers[iloc]]
               ##z[iindiv,2*iloc-1] <- rdiscr(ff)
               ##z[iindiv,2*iloc]   <- rdiscr(ff)
               z[iindiv,c(2*iloc-1,2*iloc)] <- sample(x=1:allele.numbers[iloc],
                                                      size=2,
                                                      replace=TRUE,
                                                      prob=ff)
             }
         }
       print("genotypes have been simulated")
       
       freq.indiv <- array(data=freq[,1:nindiv,,],
                           dim=c(npop,nindiv,nloc,max(allele.numbers)))
       
       gf.indiv <- array(data=e[,1:nindiv,,],
                         dim=c(npop,nindiv,nloc,max(allele.numbers)))
       
       
       if(give.freq.grid)
         {
           freq.grid <- array(data=freq[,-(1:nindiv),,],
                              dim=c(npop,prod(npix),
                                nloc,max(allele.numbers)))
           gf.grid <- array(data=e[,-(1:nindiv),,],
                            dim=c(npop,prod(npix),
                              nloc,max(allele.numbers)))
         }

    
       ## Fstat
       Fst <- Fis <- NA
       if(comp.Fst)
         {
           print("Computing Fstat")
           Fstat.res <- Fstat(npop=npop,
                              genotypes=z,
                              pop.mbrship=color.nuclei[nearest.nucleus.indiv])
           Fst <-    Fstat.res$Fst
           Fis <- Fstat.res$Fis
           print("End of call to Fstat")
         }
       
       ## IBD index Dsigma^2
       Dsigma2 <- rep(NA,npop)
       if(comp.Dsigma2)
         {
           print("Computing IBD index Dsigma^2")
           for(ipop in 1:npop)
             {
               sub.pop <- color.nuclei[nearest.nucleus.indiv]==ipop
               size.pop <- sum(sub.pop)
               a <- matrix(nrow=size.pop,ncol=size.pop,data=-999)
               if(size.pop > 0)
                 { 
                   res <- .Fortran(##name="areq4",
                                   "areq4",
                                   PACKAGE="Geneland",
                                   as.integer(size.pop),
                                   as.integer(nloc),
                                   as.integer(nloc*2),
                                   as.integer(allele.numbers),
                                   as.integer(max(allele.numbers)),
                                   as.integer(z[sub.pop,]),
                                   as.double(a))
                   a <- matrix(nrow=size.pop,ncol=size.pop,res[[7]])
                   sub.upt <- upper.tri(a)
                   aa <- a[sub.upt]
                   r <- as.matrix(dist(t(coord.indiv[,sub.pop]),upper=TRUE))
                   rr <- r[sub.upt]
                   coeff <- lm(aa~log(rr))$coefficients
                   Dsigma2[ipop] <- 1/(4*pi*coeff[2])
                 }
             }
           print("End of computation of IBD index Dsigma^2")
         }
       
       
       ## Local diff accross the barrier between pop
       diff.B <- matrix(nrow=npop,ncol=npop,NA)
       if(comp.diff)
         {
           pop.mbrshp=color.nuclei[nearest.nucleus.indiv]
           if(plot.pairs.borders)
             {
             dev.new()
             plot(t(coord.indiv),type="n")
             }
           for(ipop1 in 1:(npop-1))
             {
               sub1 <- (1:nindiv)[pop.mbrshp == ipop1]
               if(length(sub1) > 0)
                 {
                   nindiv1 <- length(sub1)
                   for(ipop2 in (ipop1+1):npop)
                     {
                       sub2 <- (1:nindiv)[pop.mbrshp == ipop2]
                       if(length(sub2) > 0)
                         {
                           nindiv2 <- length(sub2)
                           dd <- as.matrix(dist(t(coord.indiv[,c(sub1,sub2)])))
                           diag(dd) <- Inf
                           npairs <- sum(dd[1:nindiv1,(nindiv1+1):(nindiv1+nindiv2)] < width)
                           ##print(npairs)
                           ind.pairs <- matrix(nrow=npairs,ncol=2)
                           k <- 1
                           for(iindiv1 in 1:nindiv1)
                             {
                               for(iindiv2 in 1:nindiv2)
                                 {
                                   if(dd[iindiv1,nindiv1+iindiv2] < width)
                                     {
                                       ind.pairs[k,] <- t(c(sub1[iindiv1],sub2[iindiv2]))
                                       k <- k+1
                                     }
                                 }
                             }
                           diff.B[ipop1,ipop2] <- mean(abs(freq[ipop1,ind.pairs[,1],,]-
                                                           freq[ipop2,ind.pairs[,2],,]))
                           if(plot.pairs.borders)
                             {
                               
                               points(t(coord.indiv[,ind.pairs[,1]]),col=ipop1)
                               points(t(coord.indiv[,ind.pairs[,2]]),col=ipop2)
                             }
                         }
                     }
                 }
             }
         }
       
       ## Local differentiation within pop
       diff.W <- rep(NA,npop)
       if(comp.diff)
         {
           pop.mbrshp=color.nuclei[nearest.nucleus.indiv]
           for(ipop in 1:npop)
             {
               sub1 <- (1:nindiv)[pop.mbrshp == ipop]
               if(length(sub1) > 0)
                 {
                   nindiv1 <- length(sub1)
                   dd <- as.matrix(dist(t(coord.indiv[,sub1])))
                                        #print(dd)
                   diag(dd) <- Inf
                   npairs <- sum(dd < width)/2
                                        #print(npairs)
                   ind.pairs <- matrix(nrow=npairs,ncol=2)
                   k <- 1
                   for(iindiv1 in 1:(nindiv1-1))
                     {
                       for(iindiv2 in (iindiv1+1):nindiv1)
                         {
                           if(dd[iindiv1,iindiv2] < width)
                             {
                               ind.pairs[k,] <- t(c(sub1[iindiv1],sub1[iindiv2]))
                               k <- k+1
                             }
                         }
                     }
                   diff.W[ipop] <- mean(abs(freq[ipop,ind.pairs[,1],,]-
                                            freq[ipop,ind.pairs[,2],,]))
                 }
             }
         }
       
       
       ## mean heterozygoty (average over indiv and loci)
       {
         iall1 <- seq(1,2*nloc-1,2)
         iall2 <- seq(2,2*nloc,2)
         heter <- mean(z[,iall1]!=z[,iall2])
       }
     }

    if(sim.quanti)
      {
        quanti.var <- matrix(nrow=nindiv,ncol=nquanti.var)
        for(iindiv in 1:nindiv)
         {
           ipop <- color.nuclei[nearest.nucleus[iindiv]]
           quanti.var[iindiv,] <- rnorm(n=nquanti.var,
                                        mean=mean.quanti[ipop,],
                                        sd=sd.quanti[ipop,])
           
         }
      }

    res <- list(coord.lim=coord.lim,
                coord.indiv=t(coord.indiv),
                npop=npop,
                nearest.nucleus.indiv=nearest.nucleus.indiv,
                number.nuclei=number.nuclei,
                coord.nuclei=t(coord.nuclei),
                color.nuclei=color.nuclei)

    if(sim.gen)
      {
        res <- c(res,
                list(allele.numbers=allele.numbers,
                     model=model,
                     alpha=alpha,
                     beta=beta,
                     gamma=gamma,
                     genotypes=z,
                     allele.numbers=allele.numbers,
                     freq.indiv=freq.indiv,
                     gf.indiv=gf.indiv,
                     Fst=Fst,
                     Fis=Fis,
                     Dsigma2=Dsigma2,
                     heter=heter,
                     diff.B=diff.B,
                     diff.W=diff.W)
                 )
      }
    if(sim.quanti)
      {
        res <- c(res,
                 list(mean.quanti=mean.quanti,
                      sd.quanti= sd.quanti,
                      quanti.var=quanti.var)
                 )
      }
    
    if(give.tess.grid | give.freq.grid )
      {
        res <- c(res,
                 list(coord.grid=t(coord.grid),
                      nearest.nucleus.grid=nearest.nucleus.grid,
                      npix=npix))
      }
    if(give.freq.grid)
      {
        res <- c(res,
                 list(freq.grid=freq.grid,         
                      gf.grid=gf.grid)) 
      }
    return(res)
  }

