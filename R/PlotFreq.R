#' @title PlotFreq
#' @description  Plot the trace of allele frequencies in present time population
#'   and optionnaly print it.
#'   for allele \code{iall}  of locus \code{iloc}
#' 
#' @usage  PlotFreq(path.mcmc,ipop,iloc,iall,printit=FALSE,path)
#
#'  @param path.mcmc  Path to output files directory 
#'  @param ipop Integer number : index of population
#'  @param iloc Integer number : index of locus
#'  @param iall Integer number : index of allele. If \code{MCMC}
#'     was launched with option \code{filter.null.alleles=TRUE}, an extra
#'     fictive allele standing for putative null alleles is created. It
#'     estimated frequency can be also plotted. If there is say, 5 alleles
#'     at a locus, the estimated frequency of null alleles can be seen
#'     invoking \code{PlotFreq} with \code{iall=6}.
#'  @param printit Logical : if TRUE, figures are also printed
#'  @param path Character : Path to directory where figures
#'     should be printed
#' @export
PlotFreq <- function(path.mcmc,ipop,iloc,iall,printit=FALSE,path)
  {
                                        # get informations about the MCMC run 
    fileparam <- paste(path.mcmc,"parameters.txt",sep="")
    param <- as.matrix(read.table(fileparam))
    nit <- as.numeric(param[param[,1]=="nit",3])
    thinning <-  as.numeric(param[param[,1]=="thinning",3])
    filter.null.alleles <- as.logical(param[param[,1]=="filter.null.alleles",3])
    
##     if(ploidy == 1)
##       {
##         ## diploidize the data
##         data.tmp <- matrix(nrow=nrow(genotypes),
##                            ncol=ncol(genotypes)*2)
##         data.tmp[,seq(1,ncol(genotypes)*2-1,2)] <- genotypes
##         data.tmp[,seq(2,ncol(genotypes)*2,2)] <- genotypes
##         genotypes <- data.tmp
##        }
##     allele.numbers <- FormatGenotypes(genotypes)$allele.numbers
##     if(filter.null.alleles){allele.numbers <-  allele.numbers + 1}
    allele.numbers <- c(scan(paste(path.mcmc,"allele.numbers.geno2.txt",sep="")),
                        scan(paste(path.mcmc,"allele.numbers.geno1.txt",sep="")),
                        scan(paste(path.mcmc,"number.levels.ql.txt",sep="")))
    nalmax <- max(allele.numbers)
       
    filef <-  paste(path.mcmc,"frequencies.txt",sep="")
    f <- as.matrix(read.table(filef))
    npopmax <- ncol(f)
    nloc <- length(allele.numbers)
                                        # extract frequencies
                                        # from messy matrix f
    sub1 <- rep(FALSE,times=iall-1)
    sub2 <- TRUE
    sub3 <- rep(FALSE,times=nalmax-iall)
    sub <- c(sub1,sub2,sub3)
    sub1 <- rep(FALSE,(iloc-1)*nalmax)
    sub2 <- sub
    sub3 <- rep(FALSE,(nloc-iloc)*nalmax)
    sub <- rep(c(sub1,sub2,sub3),times=nit/thinning)
    plot(f[sub,ipop],
         xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
         ylab=paste("Frequency of allele",
           iall,"at locus",iloc),type="l",ylim=c(0,1))
    if(iall < allele.numbers[iloc])
      {
        title(main=paste("Allele frequencies of allele",iall,
                         "for locus",iloc,
                         "in population",ipop))
      }
    if(filter.null.alleles & (iall == allele.numbers[iloc]))
      {
        title(main=paste("Frequencies of null alleles",
                         "for locus",iloc,
                         "in population",ipop))
      }     
    if(printit==TRUE)
      {
        postscript(file=paste(path,"freq.pop",ipop,".loc",iloc,".ps",sep=""))
        par(mfrow=c(ceiling(sqrt(allele.numbers[iloc])),ceiling(sqrt(allele.numbers[iloc]))))
                                        # extract frequencies
                                        # in this messy tabular
        sub1 <- rep(FALSE,times=iall-1)
        sub2 <- TRUE
        sub3 <- rep(FALSE,times=nalmax-iall)
        sub <- c(sub1,sub2,sub3)
        sub1 <- rep(FALSE,(iloc-1)*nalmax)
        sub2 <- sub
        sub3 <- rep(FALSE,(nloc-iloc)*nalmax)
        sub <- rep(c(sub1,sub2,sub3),times=nit/thinning)
        plot(f[sub,ipop],
             xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
             ylab=paste("Frequency of allele",
               iall,"at locus",iloc),type="l",ylim=c(0,1))
        if(iall < allele.numbers[iloc])
          {
            title(main=paste("Allele frequencies of allele",iall,
                    "for locus",iloc,
                    "in population",ipop))
          }
        if(filter.null.alleles & (iall == allele.numbers[iloc]))
          {
            title(main=paste("Frequencies of null alleles",
                    "for locus",iloc,
                    "in population",ipop))
          }     
        dev.off()
      }
  }
