# Geneland package installation instructions

R package devtools https://cran.r-project.org/web/packages/devtools/index.html

For Windows users: Rtools  https://cran.r-project.org/bin/windows/Rtools/

For installation under R: 

`if( ! 'devtools' %in% installed.packages() )  { install.packages('devtools') }`

`devtools::install_github('gilles-guillot/Geneland', build_vignettes = TRUE)`

Load the package with: `library(Geneland)`

And explore it with: `vignette('Geneland')`

