GenelandStartupMessage <- function()
{
  msg <- c(paste0(
    "Geneland   version ", 
    packageVersion("Geneland")),
    "\nType 'citation(\"Geneland\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- GenelandStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'Geneland' version", packageVersion("Geneland"))
  packageStartupMessage(msg)      
  invisible()
}