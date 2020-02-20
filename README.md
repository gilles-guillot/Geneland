# Geneland package installation instructions

## Requirements

### R package devtools

To install this package from github, you need to install first the 
R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html)

### Windows users

* Windows users need also to install the program  [Rtools](https://cran.r-project.org/bin/windows/Rtools)
* To build the documentation (vignette) it is necessary to have pdlatex installed through e.g 
[MiKTeX](https://miktex.org/download)



## Installation of Geneland

`if( ! 'devtools' %in% installed.packages() )  { install.packages('devtools') }`

`devtools::install_github('gilles-guillot/Geneland', build_vignettes = TRUE)`

or if you do not have pdflatex installed, simply: 
`devtools::install_github('gilles-guillot/Geneland')`

## Getting started

Load the package with: `library(Geneland)`

And explore it with: `vignette('Geneland')`

