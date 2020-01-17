#' Fstat.output
#' @description Computes F statistics according to Weir and Cockerham's 
#'   estimators.    Missing values are allowed but as of today, the NA code is only treated as an extra
#'   allele which might bias the result. This function should not be used on haploid data.
#'   @param coordinates Matrix with one line per individual and two columns
#'    @param genotypes Genotypes of individuals. A matrix with one line per
#'      individual and 2 columns per locus
#'    @param burnin Integer: number of saved iterations to discard.
#'    @param ploidy Integer: 1 or 2 (default is 2).
#'      Do not use for haploid data.
#'   @param path.mcmc  Path to output files directory 
#' @return A list with components 
#'  Pairwise.Fis (matrix of real numbers estimating the pairwise Fis) and 
#'  Pairwise.Fst (matrix of real numbers estimating the pairwise Fs
#' @export 

Fstat.output <- function(coordinates=NULL,genotypes,ploidy=2,burnin=NULL,path.mcmc)
  {
    param <- as.matrix(read.table(file=paste(path.mcmc,"parameters.txt",sep="")))
    post.param <- as.matrix(read.table(file=paste(path.mcmc,"postprocess.parameters.txt",sep="")))
    
    pmpi <- read.table(file=paste(path.mcmc,"modal.pop.indiv.txt",sep=""))[,3]
    npop <- scan(paste(path.mcmc,"populations.numbers.txt",sep=""))
    burnin <- as.numeric(post.param[post.param[,1]=="burnin",3])
    npop.max <- as.numeric(param[param[,1]=="npopmax",3])
    npop.est <- order(hist(x=npop[-(1:burnin)],
                           breaks=seq(.5,npop.max+.5,1),plot=FALSE)$counts,
                      decreasing=TRUE)[1]
    Fstat(genotypes,npop.est,pmpi,ploidy=ploidy)
  }
