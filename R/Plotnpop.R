#' @title  Plotnpop
#' @description Gives a plot of the number of populations along the MCMC run
#' @usage  Plotnpop(path.mcmc,burnin,printit=FALSE,file,format="pdf")
#'   @param path.mcmc  Path to output files directory 
#'   @param printit Logical : if TRUE, figures are also printed
#'   @param file Character : Path to file where figures
#'     should be printed
#'   @param format format of the output file, should be either \code{"ps"}
#'     or \code{"pdf"}
#'   @param burnin An integer: number of saved iterations to discard for
#'     the representation of the histogram of the chain
#' @export
#' 


Plotnpop <- function(path.mcmc,burnin,printit=FALSE,file,
                     format="pdf")
  {
    fileparam <- paste(path.mcmc,"parameters.txt",sep="/")
    param <- as.matrix(read.table(fileparam))
    thinning <- as.numeric(param[param[,1]=="thinning",3])
    
    filenpop <- paste(path.mcmc,"populations.numbers.txt",sep="")
    npop <- scan(filenpop)
    if(burnin >0)
      {
        sub <- -(1:burnin)
      }else
    {
      sub <- 1:length(npop)
    }
    #layout(mat=matrix(nr=1,nc=2,data=1:2),width=c(4,2))
    if(printit==TRUE)
      {
        if(format=="ps"){postscript(file)}
        if(format=="pdf"){pdf(file)}
        {
          par(mfrow=c(1,2))
          plot((1:length(npop))*thinning,npop,type="l",ylab="Number of clusters",
               xlab=paste("Index of MCMC iteration","\n Whole chain",sep=""),
               ylim=c(1,max(npop)+0.5))
          hist(npop[sub],plot=TRUE,prob=TRUE,breaks=seq(.5,max(npop)+0.5,1),
               xlab=paste("Nb. of clusters along the chain \n(after a burnin of ",
                 burnin,"x",thinning,"i t.)",sep=""),
               main="Number of clusters\n along the chain \n after burnin")
          dev.off()
        }
      }
    par(mfrow=c(1,2))
    plot((1:length(npop))*thinning,npop,type="l",ylab="Number of classes",
         xlab=paste("Index of MCMC iteration","\n Whole chain",sep=""),
         ylim=c(1,max(npop)+0.5))
    hist(npop[sub],plot=TRUE,prob=TRUE,breaks=seq(.5,max(npop)+0.5,1),
         xlab=paste("Nb. of clusters along the chain \n(after a burnin of ",
           burnin,"x",thinning,"i t.)",sep=""),
         main="Number of clusters\n along the chain \n after burnin")
  }






