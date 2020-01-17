#' @title PlotFreqA
#' @description Plot frequency of  allele \code{iall}  of locus \code{iloc}
#'   in ancestral population
#' 
#' @usage PlotFreqA(path.mcmc,iloc,iall,printit=FALSE,path)
#'  @param path.mcmc  Path to output files directory 
#'  @param iloc Integer number : index of locus
#'  @param iall Integer number : index of allele
#'  @param printit Logical : if TRUE, figures are also printed
#'  @param path Character : Path to directory where figures
#'     should be printed
#' @export
PlotFreqA <- function(path.mcmc,iloc,iall,printit=FALSE,path)
  {
   
    
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])
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
    
    nloc <- length(allele.numbers)
    filefa <-  paste(path.mcmc,"ancestral.frequencies.txt",sep="")
    fa <- as.matrix(read.table(filefa))
    dev.new()
    plot(fa[,(iloc-1)*nalmax+iall],
         xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
         ylab=paste("Frequency of allele",
           iall,"at locus",iloc),type="l",ylim=c(0,1))
    title(main=ifelse(iall==1,
            paste("Allele freq. in ancestral pop. at locus",iloc),
            ""))
    if(printit==TRUE)
      {
        postscript(file=paste(path,"freq.ancestral.pop.loc",iloc,".ps",sep=""))
        plot(fa[,(iloc-1)*nalmax+iall],
             xlab=paste("Index of MCMC iteration"," (x ",thinning,")",sep=""),
             ylab=paste("Frequency of allele",
               iall,"at locus",iloc),type="l",ylim=c(0,1))
        title(main=ifelse(iall==1,
                paste("Allele freq. in ancestral pop. at locus",iloc),
                ""))
        dev.off()
      }

  }


