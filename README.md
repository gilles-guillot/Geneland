# Geneland

Bayesian models and MCMC algorithms for inferring population genetic structure
from georeferenced individual multi-locus genotypes and quantitative phenotypes.
Individuals are clustered into populations whose spatial extent is modelled by
a coloured Poisson-Voronoi tessellation.

**Version 5.0.0** adds a Shiny interface and ggplot2 graphics. The model and
the sampler are unchanged: for a given seed, `MCMC()` reproduces the output of
4.9.2 exactly. See [NEWS.md](NEWS.md) for the full list of changes, including
several long-standing bugs.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("gilles-guillot/Geneland", build_vignettes = TRUE)
```

or, without the vignette if you have no LaTeX installation:

```r
remotes::install_github("gilles-guillot/Geneland")
```

Windows users need [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/).
Building the vignette needs `pdflatex`, from e.g.
[MiKTeX](https://miktex.org/download) or the R package
[tinytex](https://CRAN.R-project.org/package=tinytex).

### Optional packages

| Package | Needed for |
|---|---|
| `shiny`, `bslib` | the graphical interface |
| `fields` | frequency-surface plots in `show.simdata()` |
| `RandomFields` | `simdata(IBD = TRUE)` only — see below |

`RandomFields` was archived from CRAN in 2022 and is distributed via GitHub,
so it cannot be installed with `install.packages()`:

```r
remotes::install_github("cran/RandomFields")
```

Every other function works without it.

## Getting started

### Graphical interface

```r
library(Geneland)
Geneland.GUI()          # or run_geneland_app()
```

Six tabs — Data, Model, Run, Diagnostics, Maps, Results — covering the standard
workflow. Everything is written to a temporary directory; nothing touches your
filespace unless you use a download button.

The previous Tk interface is unchanged and still available as
`Geneland.GUI.tcltk()`. It covers some features the Shiny app does not yet
have: parallel multi-chain runs, F-statistics, data simulation and the
hybrid-zone tools.

### From the console

```r
library(Geneland)

coord <- as.matrix(read.table(system.file("extdata", "coordinates.txt",
                                          package = "Geneland")))
geno  <- as.matrix(read.table(system.file("extdata", "genotypes.txt",
                                          package = "Geneland")))

path <- file.path(tempdir(), "run1/")
dir.create(path)

MCMC(coordinates = coord, geno.dip.codom = geno, path.mcmc = path,
     rate.max = nrow(coord), nb.nuclei.max = 3 * nrow(coord),
     npopmin = 1, npopinit = 5, npopmax = 8,
     nit = 100000, thinning = 100,
     freq.model = "Uncorrelated", varnpop = TRUE, spatial = TRUE)

PostProcessChain(coordinates = coord, path.mcmc = path,
                 nxdom = 100, nydom = 100, burnin = 200)
```

### Figures

```r
run <- read_geneland(path)
run

autoplot(run, "npop")                          # trace of the number of clusters
autoplot(run, "npop_post")                     # its posterior distribution
autoplot(run, "map",   coordinates = coord)    # modal cluster membership
autoplot(run, "proba", coordinates = coord)    # per-cluster probability surfaces

ggplot2::ggsave("map.pdf", autoplot(run, "map"), width = 8, height = 6)
```

`autoplot()` returns a ggplot object, so it can be modified like any other:

```r
autoplot(run, "npop") + ggplot2::labs(title = "My dataset")
```

Available displays: `npop`, `npop_post`, `ntile`, `ntile_post`, `drift`,
`loglik`, `logpost`, `rate`, `map`, `proba`.

The base-graphics functions (`Plotnpop()`, `PosteriorMode()`, `PlotDrift()`
and the rest) are unchanged and still work.

For the full worked example:

```r
vignette("Geneland")
```

## Citing

```r
citation("Geneland")
```

## Building the source package

The vignette PDF must be compacted, otherwise `R CMD check --as-cran` raises a
WARNING about PDF size:

```sh
R CMD build --compact-vignettes=both Geneland
R CMD check --as-cran Geneland_5.0.0.tar.gz
```
