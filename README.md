# Geneland package installation instructions

R package devtools https://cran.r-project.org/web/packages/devtools/index.html
13
​
14
For Windows users: Rtools  https://cran.r-project.org/bin/windows/Rtools/
15
​
16
For installation under R: 
17
​
18
`if( ! 'devtools' %in% installed.packages() )
19
  { install.packages('devtools') }`
20
`devtools::install_github('gilles-guillot/Geneland', build_vignettes = TRUE)`
21
​
22
Load the package with: 
23
`library(Geneland)`
24
​
25
And explore its feature with:
26
`vignette('Geneland')`
27
