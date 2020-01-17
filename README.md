# Geneland package installation instructions

To install this package from github, you need to install first the 
R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

Windows users need also to install the program  [Rtools](https://cran.r-project.org/bin/windows/Rtools)

For installation under R: 

`if( ! 'devtools' %in% installed.packages() )  { install.packages('devtools') }`

`devtools::install_github('gilles-guillot/Geneland', build_vignettes = TRUE)`

Load the package with: `library(Geneland)`

And explore it with: `vignette('Geneland')`

