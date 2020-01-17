#'
#' @title MCMC inference in Geneland
#' @description  Markov Chain Monte-Carlo inference of clusters from genotpype data
#' @param coordinates Spatial coordinates of individuals. A matrix with 2 columns and one line per individual.
#'
#' @param geno.dip.codom Genotypes  for diploid data with codominant markers.
#' A matrix with one line per individual and two  columns per locus.
#' Note that the object has to be of type matrix not table. This can be
#' forced by function \code{as.matrix}.
#' @param geno.dip.dom  Genotypes for diploid data with dominant
#' markers. A matrix with one line per individual and one  column per
#' locus. Presence/absence of a band should be
#' coded as 0/1 (0 for absence / 1 for presence). Dominant and
#' codominant markers can be analyzed jointly by passing variables to arguments 
#' geno.dip.codom and geno.dip.dom.
#' Haploid data and diploid dominant data can not be analyzed jointly in
#' the current version.
#' Note that the object has to be of type matrix not table. This can be
#' forced by function \code{as.matrix}.
#' @param geno.hap Genotypes of haploid data.
#' A matrix with one line per individual and one  column per
#' locus.     Dominant diploid data and haploid data 
#' can be analyzed jointly (e.g. to analyse microsatelite data or SNP
#'                          data together with mtDNA.
#'                          Haploid data and diploid dominant data can not be analyzed jointly in
#'                          the current version.
#'                          Note that the object has to be of type matrix not table. This can be
#'                          forced by function \code{as.matrix}.
#' @param qtc A matrix of continuous quantitative  phenotypic variables. One line per individual and one  column per
#' phenotypic variable. 
#' Note that the object has to be of type matrix not table. This can be
#' forced by function \code{as.matrix}.
#' @param qtd A matrix of discrete quantitative phenotypic variables. NOT IMPLEMENTED YET
#' @param ql A matrix of categorical phenotypic variables. NOT IMPLEMENTED YET
#' @param path.mcmc Path to output files directory. It seems that the
#' path should be given in the Unix style even under Windows (use \/
#'                                                              instead of \\).
#' This path *has to*  end with a slash (\/)
#' (e.g. path.mcmc="/home/me/Geneland-stuffs/")
#' @param rate.max Maximum rate of Poisson process  (real number >0).
#' Setting \code{rate.max} equal to the number of individuals in the
#' dataset has proved to be efficient in many cases.
#' @param delta.coord Parameter prescribing the amount of uncertainty attached
#' to spatial coordinates. If \code{delta.coord}=0 spatial coordinates are
#' considered as true coordinates, if \code{delta.coord}>0 it is assumed that observed
#' coordinates are true coordinates blurred by an additive noise uniform
#' on a square of side \code{delta.coord} centered on 0.
#' @param shape1 First parameter in the Beta(shape1,shape2) prior
#' distribution of the drift parameters in the Correlated model.
#' @param shape2 Second parameter in the Beta(shape1,shape2) prior
#' distribution of the drift parameters in the Correlated model.
#' @param npopmin Minimum number of populations (integer >=1) 
#' @param npopinit Initial number of populations
#' ( integer sucht that
#'   \code{npopmin} =< \code{npopinit} =< \code{npopmax})
#' @param npopmax Maximum number of populations (integer >=
#' \code{npopinit}).
#' There is no obvious rule to select \code{npopmax},
#' it should be set to a value larger than any value that
#' you can reasonably expect for your data.
#' @param nb.nuclei.max Integer: Maximum number of nuclei in the
#' Poisson-Voronoi tessellation. A good guess consists in setting this
#' value equal to \code{3*rate.max}. Lower values
#' can also be used in order to speed up computations. The relevance of
#' the value set can be
#' checked by inspection of the MCMC run. The number of tiles should not
#' go too close to \code{nb.nuclei.max}. If it does, you should re-run your
#' chain  with a larger value for \code{nb.nuclei.max}. In case of use
#' of the option \code{SPATIAL=FALSE}, \code{nb.nuclei.max} should be
#' set equal to the number of individuals.
#' @param nit Number of MCMC iterations
#' @param thinning Number of MCMC iterations between two writing steps (if \code{thinning}=1, all
#' states are saved whereas if e.g. \code{thinning}=10 only each 10 iteration is saved)
#' @param freq.model Character: "Correlated" or "Uncorrelated" (model for
#' frequencies). 
#' See also details in detail section of \code{\link{Geneland}} help page.
#' @param varnpop Logical: if TRUE the number of class is treated as
#' unknown and will vary along the MCMC inference, if FALSE it will be
#' fixed to the initial value \code{npopinit}. 
#' @param spatial Logical: if TRUE the colored Poisson-Voronoi
#' tessellation is used as a prior for the spatial organisation of
#' populations. If FALSE, all clustering receive equal prior
#' probability. In this case spatial information (i.e coordinates)
#' are not used  and the locations of  the nuclei are initialized and
#' kept fixed at the locations of individuals.
#' @param jcf Logical: if true update of c and f are performed jointly
#' @param filter.null.alleles Logical: if TRUE, tries to filter out null
#' alleles. An extra fictive null allele is created at each locus coding
#' for all putative null allele. Its frequency  is estimated and can be
#' viewed with function \code{PlotFreq}. This option is available only
#' with \code{freq.model="Uncorrelated"}.
#' @param prop.update.cell Integer between 0 and 1. Proportion of cell updated. For
#' debugging only.
#' @param write.rate.Poisson.process Logical: if TRUE (default) write rate
#' of Poisson process simulated by MCMC
#' @param write.number.nuclei Logical: if TRUE (default) write number of nuclei simulated by MCMC
#' @param write.number.pop Logical: if TRUE (default) write number of populations simulated by MCMC 
#' @param write.coord.nuclei Logical: if TRUE (default) write coordinates
#' of nuclei simulated by MCMC 
#' @param write.color.nuclei Logical: if TRUE (default) write color of
#' nuclei simulated by MCMC
#' @param write.freq Logical: if TRUE (default is FALSE) write allele
#' frequencies simulated by MCMC
#' @param write.ancestral.freq Logical: if TRUE (default is FALSE) write
#' ancestral allele frequencies simulated by MCMC
#' @param write.drifts Logical: if TRUE (default is FALSE) write drifts simulated by MCMC 
#' @param write.logposterior Logical: if TRUE (default is FALSE) write
#' logposterior simulated by MCMC 
#' @param write.loglikelihood Logical: if TRUE (default is FALSE) write
#' loglikelihood simulated by MCMC
#' @param write.true.coord Logical: if TRUE (default is FALSE) write true
#' spatial coordinates simulated by MCMC 
#' @param write.size.pop Logical: if TRUE (default is FALSE) write size of
#' populations simulated by MCMC
#' @param write.mean.quanti Logical: if TRUE (default is FALSE) write
#' means of quantitatives variables in the various groups simulated by MCMC
#' @param write.sd.quanti Logical: if TRUE (default is FALSE) write
#' standard deviations of quantitatives variables in the various groups simulated by MCMC
#' @param write.betaqtc Logical: if TRUE (default is FALSE) write
#' hyper-parameter beta of distribution of quantitatives variables simulated by MCMC
#' @param miss.loc A matrix with \code{nindiv} lines and \code{nloc}
#' columns of 0 or 1. For each individual, at each locus it says if the
#' locus is genuinely missing (no attempt to measure it). This info is
#' used under the option \code{filterNA=TRUE} do decide how a double
#' missing value should be treated (genuine missing data
#'                                  or double null allele).
#' @return              Successive states of all blocks of parameters are written in files
#' contained in \code{path.mcmc} and named after the type of parameters they contain.
#' All parameters processed by function \code{\link{MCMC}} are
#' written in the directory  specified by \file{path.mcmc} as follows:
#'   \itemize{
#'     \item File \file{population.numbers.txt} contains values of the number of
#'     populations (\code{nit} lines, one line per iteration of the MCMC
#'                   algorithm).
#'     
#'     
#'     \item File \file{population.numbers.txt} contains values of the number of
#'     populations (\code{nit} lines, one line per iteration of the MCMC algorithm).
#'     
#'     \item  File \file{nuclei.numbers.txt} contains the number of points in the Poisson
#'     point process generating the Voronoi tessellation.
#'     
#'     \item  File \file{color.nuclei.txt} contains vectors of integers of
#'     length \code{nb.nuclei.max} coding the class membership of each Voronoi tile.
#'     Vectors of class membership for successive states of the chain are
#'     concatenated in one column. Some entries of the vector containing
#'     clas membership for a current state may have missing values as the
#'     actual number of polygon may be smaller that the maximum number allowed
#'     \code{nb.nuclei.max}. This file has \code{nb.nuclei.max*chain/thinning} lines.
#'     
#'     \item  File \file{coord.nuclei.txt} contains coordinates of points in the Poisson
#'     point process generating the Voronoi tessellation. It has
#'     \code{nb.nuclei.max*chain/thinning} lines
#'     and two columns (hor. and vert. coordinates).
#'     
#'     \item  File \file{drifts.txt} contains the drift factors for each
#'     population, (one column per population).
#'     
#'     \item  File \file{ancestral.frequencies.txt} contains allele frequencies in ancestral
#'     population. Each line contains all frequencies of the current state.
#'     The file has \code{nit} lines.
#'     In each line, values of allele frequencies are stored by increasing
#'     allele index and and locus index (allele index varying first).
#'     
#'     \item  File \file{frequencies.txt}contains allele frequencies of present time
#'     populations. Column xx contains frequencies of population numer xx.
#'     In each column values of allele frequencies are stored by increasing
#'     allele index and and locus index (allele index varying first), and
#'     values of successive iterations are pasted.
#'     The file has \code{nallmax*nloc*nit/thinning} lines where \code{nallmax} 
#'     is the maximum number of alleles over all loci.
#'     
#'     \item  File \file{Poisson.process.rate.txt} contains rates of Poisson
#'     process.
#'     
#'     \item  File \file{hidden.coord.txt} contains the coordinates of each
#'     individual as updated along the chain if those given as input are not
#'     considered as exact coordinates (which is specified by 
#'                                      \code{delta.coord} to a non zero value).
#'     
#'     \item  File \file{log.likelihood.txt} contains log-likelihood of data
#'     for the current state of parameters of the Markov chain.
#'     
#'     \item  File \file{log.posterior.density.txt} contains log of posterior probability
#'     (up to marginal density of data) of the
#'     current state of parameters in the Markov chain.
#'     
#'   }
#'  @references 
#'         \itemize{
#' \item G. Guillot, Estoup, A., Mortier, F. Cosson, J.F. A spatial statistical
#' model for landscape genetics. Genetics, 170, 1261-1280, 2005.
#' 
#' \item  G. Guillot, Mortier, F., Estoup, A. Geneland : A program for landscape
#' genetics. Molecular Ecology Notes, 5, 712-715, 2005.
#' 
#' 
#' \item   Gilles Guillot, Filipe Santos  and Arnaud Estoup,
#' Analysing georeferenced population genetics data with
#' Geneland: a new algorithm to deal with null alleles and a friendly graphical user interface
#' Bioinformatics 2008 24(11):1406-1407.
#' 
#' \item G. Guillot.  Inference of structure in subdivided populations at
#' low levels of genetic differentiation. The correlated allele
#' frequencies 
#' model revisited. Bioinformatics, 24:2222-2228, 2008
#' 
#' \item G. Guillot and F. Santos  A computer program to simulate
#' multilocus 
#' genotype data with spatially auto-correlated allele frequencies.
#' Molecular Ecology Resources, 2009
#' 
#' \item G. Guillot, R. Leblois, A. Coulon, A. Frantz  Statistical
#' methods in spatial genetics, Molecular Ecology, 2009.
#' 
#' \item B. Guedj and G. Guillot. Estimating the location and shape of hybrid
#' zones. Molecular Ecology Resources, 11(6) 1119-1123, 2011
#' 
#' \item G. Guillot, S. Renaud, R. Ledevin, J. Michaux and J. Claude. A
#' Unifying Model for the Analysis of Phenotypic, Genetic and Geographic
#' Data. Systematic Biology, to appear, 2012.}
#' @export
MCMC <- function(
                 ## input data
                 coordinates=NULL, # spatial coordinates
                 geno.dip.codom=NULL, # diploid codominant markers
                                      # one line per indiv.
                                      # two column per marker
                 geno.dip.dom=NULL, # diploid dominant markers
                                    # one line per indiv.
                                    # one column per marker
                 geno.hap=NULL, # haploid
                                # one line per indiv.
                                # one column per marker
                 qtc=NULL, # quantitative continuous variables
                 qtd=NULL, # quantitative discrete variables
                 ql=NULL, # qualitative variables
                 ## path to output directory
                 path.mcmc,
                 ## hyper-prior parameters
                 rate.max,delta.coord=0,shape1=2,shape2=20,
                 npopmin=1,npopinit,npopmax,
                 ## dimensions
                 nb.nuclei.max,
                 ## mcmc computations options
                 nit,
                 thinning=1,
                 freq.model="Uncorrelated",
                 varnpop=TRUE,
                 spatial=TRUE,
                 jcf=TRUE,
                 filter.null.alleles=TRUE,
                 prop.update.cell=0.1,
                 ## writing mcmc output files options 
                 write.rate.Poisson.process=FALSE,
                 write.number.nuclei=TRUE,
                 write.number.pop=TRUE,
                 write.coord.nuclei=TRUE,
                 write.color.nuclei=TRUE,
                 write.freq=TRUE,
                 write.ancestral.freq=TRUE, 
                 write.drifts=TRUE,
                 write.logposterior=TRUE,
                 write.loglikelihood=TRUE,
                 write.true.coord=TRUE,
                 write.size.pop=FALSE,
                 write.mean.quanti=TRUE,
                 write.sd.quanti=TRUE,
                 write.betaqtc=FALSE,
                 miss.loc=NULL)
  {

    ploidy <- 2 ## default init
    if(!is.null(geno.dip.dom) & !is.null(geno.hap))
      {
        stop("It is currently not possible to analyze jointly diploid dominant genotypes and haploid genotypes")
      }
    if(is.null(geno.dip.dom) & is.null(geno.hap))  ploidy <- 2
    if(!is.null(geno.dip.dom)) ploidy <- 2
    if(!is.null(geno.hap)) ploidy <- 1
    ## define geno1 and geno2
    ## geno1 has L columns
    ## geno2 has 2L columns
    geno2 <- geno.dip.codom
    ## print(paste("is.null(geno.dip.codom)",is.null(geno.dip.codom)))
    ## print(geno2)
    if(ploidy==2)
      {
        geno1 <- geno.dip.dom
      }
    if(ploidy==1)
      {
        geno1 <- geno.hap
      }
 
    
    
    ######################################
    ## checking most common errors on parameters
    ## path.mcmc ends with a /
    if(substring(path.mcmc,
                 first=nchar(path.mcmc),
                 last=nchar(path.mcmc)) != "/")
      {
        path.mcmc <- paste(path.mcmc,"/",sep="")
      }
    short.path <- substring(path.mcmc,first=1,last=nchar(path.mcmc)-1)
    if(!file_test("-d",short.path)) stop(paste("Directory ", path.mcmc,
                                          "does not exist."))

    if((nit %% thinning) != 0)stop('nit/thinning is not an integer')
    if(missing(npopmax)) stop('Argument npopmax is missing with no default')
    if(missing(npopinit)) npopinit <- npopmax
    if(npopinit > npopmax) stop('npopinit > npopmax')
    if(npopinit < npopmin) stop('npopinit < npopmin')
    if((freq.model != "Correlated") & (freq.model != "Uncorrelated"))
      {
        stop(paste('Error:',freq.model, 'is not a frequency model. Check spelling (case sensitive) '))
      }
    if((ploidy != 1) & (ploidy != 2))
      {
        ##print(paste("ploidy = ",ploidy," is not a valid value."))
        stop(paste("ploidy = ",ploidy," is not a valid value."))
      }
    if(!is.null(geno1) & is.null(geno2))
      {
        if(filter.null.alleles)
          {
            if(ploidy == 1)
              {
                stop("Algorithm for filtering null alleles not compatible with haploid data")
              }
            if(ploidy == 2)
              {
                stop("Algorithm for filtering null alleles not compatible with dominant markers")
              }
          }
      }
    if(!is.null(geno1) & (ploidy==2))
      ##recode codominant markers data passed as 0/1 into 1/2
       {
         geno1 <- geno1 + 1
       }

    ## read dimensions of data matrices
    if(is.null(geno1))
      {
        nloc.geno1 <- 1
        nindiv.geno1 <- 0
      }else
    {
      nloc.geno1 <- ncol(geno1)
      nindiv.geno1 <- nrow(geno1)
      print(c("in MCMC.R nindiv.geno1=",nindiv.geno1))
    }
    if(is.null(geno2))
      {
        nloc.geno2 <- 1
        nindiv.geno2 <- 0
        print(c("in MCMC.R nindiv.geno2=",nindiv.geno2))
      }else
    {
      nloc.geno2 <- ncol(geno2)/2
      nindiv.geno2 <- nrow(geno2)
    }
    if(is.null(qtc))
      {
        nqtc <- 1
        nindivqtc <- 0
      }else
    {
      nqtc <- ncol(qtc)
      nindivqtc <- nrow(qtc)
    }
    if(is.null(qtd))
      {
        nqtd <- 1
        nindivqtd <- 0
      }else
    {
      nqtd <- ncol(qtd)
      nindivqtd <- nrow(qtd)
    }
    if(is.null(ql))
      {
        nql <- 1
        nindivql <- 0
      }else
    {
      nql <- ncol(ql)
      nindivql <- nrow(ql)
    }
    
    ## check consistency of dimensions
    nnn <- c(nindiv.geno1,nindiv.geno2,nindivqtc,nindivqtd,nindivql)
    sub <- nnn>0
    if(length(unique(nnn[sub])) > 1)
      {
        print(paste('nindiv.geno1 = ',nindiv.geno1))
        print(paste('nindiv.geno2 = ',nindiv.geno2))
        print(paste('nindivqtc = ',nindivqtc))
        print(paste('nindivqtd = ',nindivqtd))
        print(paste('nindivql = ',nindivql))
        stop('Number of rows of data matrices do not match')
      }else
    { nindiv <- nnn[sub][1] }
    ## say if arrays should be used or just treated as dummy
    use.geno1 <- nindiv.geno1 > 0
    use.geno2 <- nindiv.geno2 > 0
    use.qtc <- nindivqtc > 0
    use.qtd <- nindivqtd > 0
    use.ql <- nindivql > 0


    
    ## define dummy arrays
    print('defining dummy data arrays')
    if(nindiv.geno1==0)
      { geno1 <- matrix(nrow=nindiv,ncol=1,data=NA) }
    if(nindiv.geno2==0)
      { geno2 <- matrix(nrow=nindiv,ncol=2,data=NA) }
    if(nindivqtc==0)
      { qtc <- matrix(nrow=nindiv,ncol=nqtc,data=NA) }
    if(nindivqtd==0)
      { qtd <- matrix(nrow=nindiv,ncol=nqtd,data=NA) }
    if(nindivql==0)
      { ql <- matrix(nrow=nindiv,ncol=nql,data=NA) }
    
    
    ## define dummy coordinates if coord are missing
    print('defining dummy coordinates if coord are missing')
    if(is.null(coordinates))
      {
        if(spatial)
          {
            stop('Please give spatial coordinates of individuals or set argument spatial to FALSE')
          }else
        {
          n.int <- ceiling(sqrt(nindiv))
          x <- rep(seq(from=0,to=1,length=n.int),n.int)
          y <- rep(seq(from=0,to=1,length=n.int),n.int)
          y <- as.vector(t(matrix(nrow=n.int,ncol=n.int,y,byrow=FALSE)))
          coordinates <- cbind(x,y)[1:nindiv,]
        }
      }else
    {
      if(ncol(coordinates) != 2)stop('matrix of coordinates does not have 2 columns')
      if(nrow(coordinates) != nindiv)
        {
          print(paste('number of individuals in data matrix =',nindiv))
          print(paste('number of individuals in coordinate matrix =',nrow(coordinates)))
          stop('Number of rows in coordinate matrix and data matrices do not match ')
        }     
    }
    ## reformatting coord data as a matrix
    if(!is.matrix(coordinates)) coordinates <- as.matrix(coordinates)
     
    ## define dummy matrix indicating genuinely missing data in geno2
    print('defining dummy matrix indicating genuinely missing data in geno2')
    if(is.null(miss.loc))
      {
        miss.loc <- matrix(nrow=nindiv,
                           ncol=ncol(geno2)/2,
                           data=0)
      }
     
 
    
    ## #############################################
    ##
    ## Initializing variables
    ##
    ## ############################################

    
    ## default values for rate.max and nb.nuclei.max
    if(missing(rate.max)) rate.max <- nindiv
    if(missing(nb.nuclei.max))
      {
        nb.nuclei.max <-  ifelse(spatial, 2*nindiv,nindiv)
      }
    if(nb.nuclei.max < nindiv)
      {
        stop('nb.nuclei.max should be at least equal to the number of individuals')
      }
    if(spatial & (nb.nuclei.max < 2*rate.max))
      {
        stop('nb.nuclei.max is too small as compared to rate.max')
      }

    ## reformat alleles and qualit. variables as consecutive integers
    if(nindiv.geno1 >0)
      { res <- FormatGenotypes(as.matrix(geno1),ploidy=1)
        geno1 <- res$genotypes
        allele.numbers.geno1 <- res$allele.numbers
        print(c("In R function MCMC, allele.numbers.geno1=",allele.numbers.geno1))
        if(sum(allele.numbers.geno1==1) > 0){stop('Some of the markers do not display any polymorphism')}        
        ##         if(ploidy == 1)
        ##           {
        ##             ## dealing with  apparently non polymorphic loci
        ##             print("Incrementing allele.numbers.geno1")
        ##             allele.numbers.geno1 <- allele.numbers.geno1 + 1 
        ##           }
        ##         if(ploidy == 2)
        ##           {
        ##             ## dealing with  apparently non polymorphic loci
        ##             print("Incrementing allele.numbers.geno1 ")
        ##             allele.numbers.geno1[allele.numbers.geno1==1] <- 2
        ##           }
        print(paste('Number of missing data in matrix geno1: ',sum(is.na(geno1))))
      }else
    {allele.numbers.geno1 <- -999}
    
    if(nindiv.geno2 > 0)
      { res <- FormatGenotypes(as.matrix(geno2),ploidy=2)
        geno2 <- res$genotypes
        allele.numbers.geno2 <- res$allele.numbers
        print(c("In R function MCMC, allele.numbers.geno2=",allele.numbers.geno2))
        if(sum(allele.numbers.geno2==1) > 0){stop('Some of the markers do not display any polymorphism')}
        ##         ## dealing with  apparently non polymorphic loci
        ##         print("Incrementing allele.numbers.geno2")
        ##         allele.numbers.geno2[allele.numbers.geno2==1] <- 2
        print(paste('Number of missing data in matrix geno2: ',sum(is.na(geno2))))
      }else
    { allele.numbers.geno2 <- -999}

    if(nindivql > 0)
      { res <- FormatGenotypes(as.matrix(ql),ploidy=1)
        ql <- res$genotypes
        allele.numbers.ql <- res$allele.numbers
        print(paste('Number of missing obs. for qualitative variables: ',sum(is.na(ql))))
      }else
    { allele.numbers.ql <- -999}

    if(nindivqtc > 0)
      {
        print(paste('Number of missing obs. for quantitative continuous variables: ',sum(is.na(qtc))))
      }

    if(nindivqtd > 0)
      {
        print(paste('Number of missing obs. for quantitative discrete variables: ',sum(is.na(qtd))))
      }
    nchar.path <- nchar(path.mcmc)
    
    
    ## Add a fictive allele (coding for several potential null alleles)
    if(filter.null.alleles)
      {allele.numbers.geno2 <-  allele.numbers.geno2 + 1}
    nalt<- c(allele.numbers.geno2,allele.numbers.geno1,allele.numbers.ql)
    nalmax <- max(1,max(nalt))

    ## define working arrays for parameters of spatial model
    print("Defining working arrays for parameters of spatial model...")
    u <- matrix(nrow=2,ncol=nb.nuclei.max,data=-999)
    utemp <- matrix(nrow=2,ncol=nb.nuclei.max,data=-999)
    c <- rep(times=nb.nuclei.max,-999)
    ctemp <- rep(times=nb.nuclei.max,-999)
    t <-matrix(nrow=2,ncol=nindiv,data=-999)
    ttemp <-matrix(nrow=2,ncol=nindiv,data=-999)
    indcell <- rep(times=nindiv,-999)
    indcelltemp <- rep(times=nindiv,-999)
    distcell <- rep(times=nindiv,-999)
    distcelltemp <- rep(times=nindiv,-999)
    xlim <- ylim <- rep(-999,times=2)
    
    ## define working arrays for parameters of genetic model(s)
    print("Defining working arrays for parameters of genetic model...")
    ncolt <- nloc.geno2 + nloc.geno1 + nql
    f <- array(dim=c(npopmax,ncolt,nalmax),data=-999)
    ftemp <- array(dim=c(npopmax,ncolt,nalmax),data=-999)
    fa <- array(dim=c(ncolt,nalmax),data=-999)
    drift <- rep(-999,npopmax)
    drifttemp <- rep(-999,npopmax)
    n <-  array(dim=c(npopmax,ncolt,nalmax),data=-999)
    ntemp <-  array(dim=c(npopmax,ncolt,nalmax),data=-999)
    a <- rep(times=nalmax,-999)
    ptemp <- rep(times=nalmax,-999)
    cellclass <- rep(times=nb.nuclei.max,-999)
    listcell <- rep(times=nb.nuclei.max,-999)
    fmodel <- ifelse(freq.model=="Correlated",1,0) # Correlated or Uncorrelated model
    kfix <- 1-as.integer(varnpop)
    full.cond.y <- matrix(nrow=nalmax,ncol=2,0)

    ## define working arrays for parameters of quantitative variables
    print("Defining working arrays for parameters of quantitative variables...")
    meanqtc <- sdqtc <- meanqtctmp <- sdqtctmp <- matrix(nrow=npopmax,ncol=nqtc,data=-999)
    nnqtc <- sqtc <- ssqtc <- matrix(nrow=npopmax,ncol=nqtc,data=-999)
    ## define and init working arrays for hyper-parameters of quantitative variables
    ksiqtc <- kappaqtc <- alphaqtc <- betaqtc <- gbeta <- hbeta <- rep(-999,nqtc)
    if(use.qtc)
      {
        ksiqtc <- apply(qtc,2,mean,na.rm=TRUE)
        kappaqtc <- hbeta <- 2/(apply(qtc,2,max,na.rm=TRUE) - apply(qtc,2,min,na.rm=TRUE))^2
        alphaqtc <- rep(2,nqtc)
        gbeta <- rep(.5,nqtc)
        betaqtc <- rgamma(n=nqtc,shape=gbeta,rate=hbeta)
      }
    print(paste("ksiqtc",ksiqtc))
    print(paste("kappaqtc",kappaqtc))
    print(paste("alphaqtc",alphaqtc))
    print(paste("hbeta",hbeta))
    print(paste("gbeta",gbeta))
    print(paste("betaqtc",betaqtc))
    
    ## change NA code for Fortran
    geno2.999 <- geno2
    geno2.999[is.na(geno2)] <- -999
    geno1.999 <- geno1
    geno1.999[is.na(geno1.999)] <- -999
    qtc.999 <- qtc
    qtc.999[is.na(qtc.999)] <- -999
    qtd.999 <- qtd
    qtd.999[is.na(qtd.999)] <- -999
    ql.999 <- ql
    ql.999[is.na(ql.999)] <- -999
    
    ## init true genotypes
    ## as given data for codominant markers part...
    true.geno <- cbind(geno2.999,
                        matrix(nrow=nindiv,ncol=nloc.geno1*2,data=-999))
    ## ... and as "diploidised homozigous" data for dominant markers part
    sub <- seq(nloc.geno2*2+1,nloc.geno2*2+nloc.geno1*2-1,2)
    true.geno[,sub] <- geno1.999
    sub <- seq(nloc.geno2*2+2,nloc.geno2*2+nloc.geno1*2,2)
    true.geno[,sub] <- geno1.999


    ## put together scalar parameters 
    integer.par <- c(write.rate.Poisson.process,
                     write.number.nuclei,
                     write.number.pop,
                     write.coord.nuclei,
                     write.color.nuclei,
                     write.freq,
                     write.ancestral.freq,
                     write.drifts,
                     write.logposterior,
                     write.loglikelihood,
                     write.true.coord,
                     write.size.pop,
                     write.mean.quanti,
                     write.sd.quanti,
                     write.betaqtc,
                     fmodel,
                     kfix,
                     spatial,
                     jcf,
                     filter.null.alleles,
                     ploidy,
                     nchar.path,
                     nit,
                     thinning,
                     use.geno1,
                     use.geno2,
                     use.qtc,
                     use.qtd,
                     use.ql)
    integer.par <- c(integer.par,rep(-999,100-length(integer.par)))
    double.par <-  c(rate.max,
                     delta.coord,
                     shape1,
                     shape2,
                     prop.update.cell)
    double.par <- c(double.par,rep(-999,100-length(double.par)))     

    ## define arrays for output objects
    ## written on disk from R after completion of MCMC run
    nitsaved <- nit/thinning
    
    ## out1[,1] : rate Poisson process
    ## out1[,2] : nb tiles
    ## out1[,3] : nb clusters
    ## out1[,4] : log likelhood
    ## out1[,5] : log posterior
    out1 <- array(dim=c(nitsaved,5),data=-999)

    ## outspace
    ## outspace[,1:2,] : u
    ## outspace[,3,] : c
    ## outspace[,4:5,] : t
    outspace <- array(dim=c(nitsaved,5,nb.nuclei.max),data=-999)
    
    ## outfreq
    ## outfreq[,1,,,] : f
    ## outfreq[,2,1,,] : fa
    ## outfreq[,3,,1,1] : d
    outfreq <- array(dim=c(nitsaved,3,npopmax,ncolt,nalmax),data=-999)

    ## outqtc
    ## outqtc[,1,,] : meanqtc
    ## outqtc[,2,,] : sdqtc
    outqtc <- array(dim=c(nitsaved,2,npopmax,nqtc),data=-999)

    
    ## print(c("in MCMC.R nalt=",nalt))
    
    ## Starting MCMC computations
    out.res<- .Fortran("mcmcgld",
                       PACKAGE="Geneland",
                       #
                       as.double(t(coordinates)),
                       as.integer(geno2.999),
                       as.integer(miss.loc),
                       #
                       as.integer(geno1.999),
                       #
                       as.integer(ql.999),
                       as.integer(nql),
                       #
                       as.double(qtc.999),
                       as.integer(nqtc),
                       #
                       as.character(path.mcmc),
                       as.integer(integer.par),
                       as.double(double.par),
                       as.integer(nindiv),
                       as.integer(nloc.geno2),
                       as.integer(nloc.geno1),
                       as.integer(ncolt),
                       as.integer(nalt),
                       as.integer(nalmax),
                       as.integer(nb.nuclei.max),
                       as.integer(npopinit),
                       as.integer(npopmin),
                       as.integer(npopmax),
                       as.double(xlim),
                       as.double(ylim),
                       as.integer(indcell),
                       as.integer(indcelltemp),
                       as.double(distcell),
                       as.double(distcelltemp),                      
                       as.double(t),
                       as.double(ttemp),
                       as.double(u),
                       as.double(utemp),
                       as.integer(c),
                       as.integer(ctemp),
                       #
                       as.double(f),
                       as.double(ftemp),
                       as.double(fa),
                       as.double(drift),
                       as.double(drifttemp),
                       as.integer(n),
                       as.integer(ntemp),
                       as.double(a),
                       as.double(ptemp),
                       #
                       as.double(meanqtc),
                       as.double(sdqtc),
                       as.double(meanqtctmp),
                       as.double(sdqtctmp),                       
                       as.integer(nnqtc),
                       as.double(sqtc),
                       as.double(ssqtc),
                       as.double(ksiqtc),
                       as.double(kappaqtc),
                       as.double(alphaqtc),
                       as.double(betaqtc),
                       as.double(gbeta),
                       as.double(hbeta),
                       as.integer(cellclass),
                       as.integer(listcell),
                       as.integer(true.geno),
                       as.double(full.cond.y),
                       #
                       as.integer(nitsaved),
                       as.double(out1),
                       as.double(outspace),
                       as.double(outfreq),
                       as.double(outqtc)
                       )

    ## write info on the present run in an ascii file
    param <- c(paste("nindiv :",nindiv),
               paste("rate.max :",rate.max),
               paste("nb.nuclei.max :",nb.nuclei.max),
               paste("nit :",nit),
               paste("thinning :",thinning),
               paste("varnpop :",varnpop),
               paste("npopmin :",npopmin),
               paste("npopinit :",npopinit),
               paste("npopmax :",npopmax),
               paste("spatial :",spatial),
               paste("delta.coord :",delta.coord),
               paste("use.geno1 :",use.geno1),
               paste("use.geno2 :",use.geno2),
               paste("use.qtc :",use.qtc),
               paste("use.qtd :",use.qtd),
               paste("use.ql :",use.ql)         
               )
    if(use.geno1)
      {
        param <- c(param,
                   paste("nloc.geno1 :",nloc.geno1),
                   paste("filter.null.alleles :",filter.null.alleles))
      }
    if(use.geno2)
      {
        param <- c(param,
                   paste("nloc.geno2 :",nloc.geno2),
                   paste("filter.null.alleles :",filter.null.alleles)
                   )
      }
    if((use.geno1 | use.geno2) | use.ql)
      {
        param <- c(param,
                   paste("nalmax :",nalmax),
                   paste("freq.model :",freq.model))
      }
     if(use.geno1)
      {
        param <- c(param,
                   paste("ploidy :",ploidy))
      }
    if(use.qtc)
      {
        param <- c(param,
                   paste("nqtc :",nqtc))
      }
     if(use.qtd)
      {
        param <- c(param,
                   paste("nqtd :",nqtd))
      }
     if(use.ql)
      {
        param <- c(param,
                   paste("nql :",nql))
      }       
    write.table(param,file=paste(path.mcmc,"parameters.txt",sep=""),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(allele.numbers.geno1,file=paste(path.mcmc,"allele.numbers.geno1.txt",sep=""),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(allele.numbers.geno2,file=paste(path.mcmc,"allele.numbers.geno2.txt",sep=""),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(allele.numbers.ql,file=paste(path.mcmc,"number.levels.ql.txt",sep=""),
                quote=FALSE,row.names=FALSE,col.names=FALSE)

    ##  reception mcmc outputs and write them in external text files
    print('Writing MCMC outputs in external text files')
    out1 <-  array(dim=c(nitsaved,5),data=out.res[[61]])
    write.table(file=paste(path.mcmc,"Poisson.process.rate.txt",sep=""),
                out1[,1],row.names=FALSE,col.names=FALSE)
    write.table(file=paste(path.mcmc,"nuclei.numbers.txt",sep=""),
                out1[,2],row.names=FALSE,col.names=FALSE)
    write.table(file=paste(path.mcmc,"populations.numbers.txt",sep=""),
                out1[,3],row.names=FALSE,col.names=FALSE)
    write.table(file=paste(path.mcmc,"log.likelihood.txt",sep=""),
                out1[,4],row.names=FALSE,col.names=FALSE)
    write.table(file=paste(path.mcmc,"log.posterior.density.txt",sep=""),
                out1[,5],row.names=FALSE,col.names=FALSE)
    ##
    outspace <- array(dim=c(nitsaved,5,nb.nuclei.max),data=out.res[[62]])
    for(iitstor in 1:(nit/thinning))
      {
        append <- ifelse(iitstor>1,TRUE,FALSE)
        write.table(file=paste(path.mcmc,"coord.nuclei.txt",sep=""),
                    t(outspace[iitstor,1:2,]),
                    row.names=FALSE,col.names=FALSE,append=append)
        write.table(file=paste(path.mcmc,"color.nuclei.txt",sep=""),
                    outspace[iitstor,3,],row.names=FALSE,col.names=FALSE,
                    append=append)
        write.table(file=paste(path.mcmc,"hidden.coord.txt",sep=""),
                    t(outspace[iitstor,4:5,]),row.names=FALSE,col.names=FALSE)
      }
    ##
    outfreq <- array(dim=c(nitsaved,3,npopmax,ncolt,nalmax),data=out.res[[63]])
    for(iitstor in 1:(nit/thinning))
      {
        append <- ifelse(iitstor>1,TRUE,FALSE)
        write.table(file=paste(path.mcmc,"ancestral.frequencies.txt",sep=""),
                    outfreq[iitstor,2,1,,],
                    row.names=FALSE,col.names=FALSE,append=append)
        write.table(file=paste(path.mcmc,"drifts.txt",sep=""),
                    t(outfreq[iitstor,3,,1,1]),
                    row.names=FALSE,col.names=FALSE,append=append)
        for(iloc in 1:ncolt)
          {
            append <- ifelse(iitstor>1 | iloc>1,TRUE,FALSE)
            write.table(file=paste(path.mcmc,"frequencies.txt",sep=""),
                        t(outfreq[iitstor,1,,iloc,]),
                        row.names=FALSE,col.names=FALSE,append=append)
          }
      }
    ##
    outqtc <- array(dim=c(nitsaved,2,npopmax,nqtc),data=out.res[[64]])
     for(iitstor in 1:(nit/thinning))
      {
        append <- ifelse(iitstor>1,TRUE,FALSE)
        write.table(file=paste(path.mcmc,"mean.qtc.txt",sep=""),
                    outqtc[iitstor,1,,],
                    row.names=FALSE,col.names=FALSE,append=append)
        write.table(file=paste(path.mcmc,"sd.qtc.txt",sep=""),
                    outqtc[iitstor,2,,],row.names=FALSE,col.names=FALSE,
                    append=append)
      }
  }
 
