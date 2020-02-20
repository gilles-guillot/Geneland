#' gl2gp
#' @description Takes the matrices  of  coordinates and genotypes in the Geneland format
#'   and writes it as an ascci file in the Genepop format
#'   \@param coordinates   Matrix of coordinates
#'   \@param genotypes  MAtrix of genotypes
#'   \@param file   Character string giving the path to the file
#'     where the data in Genepop format should be written
#' @return An ascci file is written. This file may need to be
#'   process in order to convert eol character from unix style to dos style. 
#' @export
gl2gp <- function(coordinates,genotypes,file)
  {
    ## recode missing values as 0
    genotypes[is.na(genotypes)] <- 0
    
    nindiv <- nrow(genotypes)
    nloc <- ncol(genotypes)/2
    write.table("Header line", file, quote = FALSE, col.names = FALSE, 
                row.names = FALSE)
    ## open file and write header line and list of loci
    for (iloc in 1:nloc)
      {
        write.table(paste("loc", iloc), file, quote = FALSE, 
                    col.names = FALSE, row.names = FALSE, append = TRUE)
      }

    ## write genotypes following genepop format
    dum.coord <- c(-999,-999)
    numb.pop <- 0
    for (iindiv in 1:nindiv) {
      if (sum(coordinates[iindiv,] != dum.coord)>0)
        {
          write.table("Pop", file, quote = FALSE, col.names = FALSE, 
                      row.names = FALSE, append = TRUE)
          dum.coord <- coordinates[iindiv,]
          numb.pop <- numb.pop + 1
        }
      string <- "sample"
      string <- paste(string,numb.pop,sep="")
      x <- c(round(coordinates[iindiv,], digits = 3))
      stringbis <- paste(x[1],x[2],sep=" ")
      stringbis <- paste("(",stringbis,")",sep="")
      string <- paste(string,stringbis,sep=" ")
      string <- c(string,",")
      for (iloc in 1:nloc)
        {
          all1 <- genotypes[iindiv, 2 * iloc - 1]
          all2 <- genotypes[iindiv, 2 * iloc]
          char1 <- "000"
          if (all1 < 10) {
            substring(char1, 3, 3) <- as.character(all1)
          }
          if ((all1 > 9) & (all1 < 100)) {
            substring(char1, 2, 3) <- as.character(all1)
          }
          if ((all1 > 99) & (all1 < 1000)) {
            substring(char1, 1, 3) <- as.character(all1)
          }
          char2 <- "000"
          if (all2 < 10) {
            substring(char2, 3,3) <- as.character(all2)
          }
          if ((all2 > 9) & (all2 < 100)) {
            substring(char2, 2, 3) <- as.character(all2)
          }
          if ((all2 > 99) & (all2 < 1000)) {
            substring(char2, 1, 3) <- as.character(all2)
          }
          char <- paste(char1, char2, sep = "")
          string <- c(string, char)
        }
      if (iindiv < nindiv) {
        write.table(x = t(string), file = file, quote = FALSE, 
                    col.names = FALSE, row.names = FALSE, append = TRUE)
      }
      if (iindiv == nindiv) {
        write.table(x = t(string), file = file, eol = "", 
                    quote = FALSE, col.names = FALSE, row.names = FALSE, 
                    append = TRUE)
      }
    }
  }

##
## geno <- as.matrix(read.table("~/Desktop/donn?es/Cod_skagrerak/4th sending/genonew2.txt"))
## coord <- as.matrix(read.table("~/Desktop/donn?es/Cod_skagrerak/4th sending/coord.new.txt"))

## gl2gp(coord,geno,"/home/tanguyr/Desktop/genepop.txt")





