# Geneland 5.0.0

A packaging and interface release. The statistical model and the MCMC sampler
are **unchanged**: for a given seed and the same arguments, `MCMC()` produces
output identical to that of version 4.9.2, verified file by file.

## Breaking changes

* `Geneland.GUI()` now launches a Shiny application instead of the Tk
  interface. The Tk interface is unchanged and remains available as
  `Geneland.GUI.tcltk()`.

* The vignette no longer ships the maintainer's full BibTeX collection
  (`inst/biblio.bib`, `inst/gilles.bib`, ~1400 entries, 353 KB). Only the 31
  references actually cited are kept, in `vignettes/Geneland.bib`.

## New: a Shiny interface

* `run_geneland_app()` starts a six-tab application — Data, Model, Run,
  Diagnostics, Maps, Results — covering the standard workflow: load data,
  choose a model, run the sampler, inspect the chain, view and export the
  maps.

* Everything is written under `tempdir()`; nothing touches your filespace
  unless you use a download button. Results can be exported as a zip archive
  and reopened later with `read_geneland()`.

* `shiny` and `bslib` are optional (`Suggests`), so the package remains
  installable without them.

* The sampler runs in the process serving the app, so the interface is
  unresponsive while a chain runs. Long chains are still better launched from
  the console with `MCMC()`.

* The application does not yet cover every feature of the Tk interface;
  parallel multi-chain runs, F-statistics, data simulation and the hybrid-zone
  tools remain available through `Geneland.GUI.tcltk()` and the console API.

## New: ggplot2 graphics

The existing `Plot*` functions are unchanged. Alongside them there is now a
ggplot2 layer that returns plot objects rather than drawing to a device, so
they can be modified, arranged and saved by the caller.

* `read_geneland()` reads a run directory into a `geneland_run` object, with
  `print()` and `plot()` methods.

* `autoplot()` on that object produces ten displays: `"npop"`, `"npop_post"`,
  `"ntile"`, `"ntile_post"`, `"drift"`, `"loglik"`, `"logpost"`, `"rate"`,
  `"map"` and `"proba"`.

* `theme_geneland()`, `geneland_pal()`, `scale_colour_geneland_d()`,
  `scale_fill_geneland_d()` and `scale_fill_geneland_c()`.

Improvements over the base-graphics equivalents:

* The colours were validated for colour-vision deficiency. `terrain.colors()`
  is not colourblind-safe and implies an ordering that clusters do not have.
  Because any two clusters can be adjacent on a map, colour alone cannot carry
  cluster identity, so the membership map also draws the cluster number.

* Posterior probability surfaces use a single-hue sequential ramp and one
  panel per cluster, replacing a loop that opened one device per cluster.

* Traces *shade* the burn-in instead of silently discarding it, so a burn-in
  that is too short stays visible.

* The posterior of the number of clusters is drawn over the whole prior
  support, so a degenerate posterior reads as "all the mass on one value"
  rather than as a single bar.

* No `dev.new()` and no changes to `par()`; the caller's graphics device is
  left alone.

## Bug fixes

* `Plotntile()` was documented but **never exported**: its `@export` tag was
  indented and roxygen treated it as prose. It is now exported.

* `gl2gp()` had its `@param` tags escaped as `\@param`, so none of its three
  arguments were documented.

* `MCMC()` passed the output path to Fortran as a character string. Passing
  character to `.Fortran` is deprecated, and the receiving code had not used
  the argument since output writing moved to R. Passing `nchar(path.mcmc)`
  unchecked into a `character*255` buffer was also an out-of-bounds read for
  paths longer than 255 characters. The argument has been removed.

* `MCMC()` extracted its results from the `.Fortran` return list by hard-coded
  positions (`[[61]]`…`[[64]]`), which silently break whenever an argument is
  added or removed. The four outputs are now named.

* `simdata(IBD = TRUE)` called `GaussRF()` from `RandomFields` with no check.
  `RandomFields` was archived from CRAN in 2022, so this failed with
  "could not find function". It now reports what is missing and how to install
  it. See "Optional dependencies" below.

* `show.simdata()` called `as.image()` and `image.plot()` from `fields`
  without declaring or checking for the package.

* Removed a leftover debug statement in `HZ()` that printed `"coucou avant
  .Fortran"` on every call.

## Optional dependencies

* `RandomFields` is needed only by `simdata(IBD = TRUE)`. It was archived from
  CRAN and is distributed via GitHub, so it cannot be declared in `Suggests`
  (there is no CRAN-style repository to install it from). Install it with
  `remotes::install_github("cran/RandomFields")`. Every other function works
  without it.

* `mapproj`, `maps` and `PBSmapping` were listed in `Suggests` but never used
  anywhere in the package; they have been removed.

## Packaging

Version 4.9.2 could not be checked at all — `R CMD check` aborted with an
ERROR. `R CMD check --as-cran` now reports 2 NOTEs on Linux.

* `DESCRIPTION`: `Description` is now a proper sentence; `tcltk` moved from
  `Depends` to `Imports`; `parallel` declared (it was used by the Tk interface
  but never declared); `knitr` declared, which was the cause of the ERROR.

* `NAMESPACE`: all base-package functions are now imported. Previously there
  were 961 undefined global references, and the namespace could not be loaded
  with its stated dependencies.

* Hand-written `@usage` blocks removed throughout. `@param` tags indented
  underneath them were being read as usage continuation lines, which produced
  every one of the malformed-`\usage` and undocumented-argument warnings.

* The vignette declares `\VignetteIndexEntry` and `\VignetteEngine` and now
  builds. Note that the source package must be built with
  `R CMD build --compact-vignettes=both`, otherwise the vignette PDF triggers
  a size WARNING.

* `inst/CITATION` uses `bibentry()` instead of the deprecated `citEntry()`.

* Compiled objects, editor backups and LaTeX intermediates are no longer
  tracked or shipped.

## Testing

The package now has a `testthat` suite: **193 assertions across 10 files**,
running in about 16 seconds.

* Numerical regression on the sampler, pinning `MCMC()` output to recorded
  values so that the pending Fortran rewrite cannot silently change the
  numbers. These are skipped on CRAN: the chain's path depends on floating
  point accept/reject decisions, so another compiler or optimisation level may
  legitimately produce a different, equally valid chain. The portable
  invariants (bounds, shapes, finiteness, parameter round-trip) are not
  skipped.
* Post-processing invariants: membership probabilities lie in [0, 1] and sum
  to one; the modal cluster is the argmax of those probabilities.
* `Fstat()` properties: the pairwise matrix is symmetric with a zero diagonal,
  and a random split of a panmictic sample gives near-zero Fst.
* The plotting layer: every display builds, the palette refuses a ninth hue,
  the membership map labels every cluster it draws, and plotting leaves `par()`
  and the device list untouched.
* The Shiny interface is driven end to end through `shiny::testServer()`:
  load, configure, sample, post-process, and render every plot.

Writing the suite turned up four defects, now fixed:

* Ten of the `write.*` arguments of `MCMC()` documented a default opposite to
  the one in the code (for example `write.drifts` is `TRUE`, documented as
  `FALSE`).
* Drift factors for clusters that are absent from the model at a given
  iteration are written as the sentinel `-999`. `autoplot(run, "drift")` was
  plotting it, so the trace plunged off the panel; the sentinel is now treated
  as missing and the line simply stops.
* The Shiny data loader assigned uploaded files to the active dataset *before*
  validating them, so a rejected upload (mismatched row counts, say) could
  still be handed to the sampler.
* `geom_label(label.size = )` is deprecated in ggplot2 3.5.0; the minimum
  ggplot2 version is now 3.5.0.

## Known limitations

* `PostProcessChain()` and `HZ()` still do their own file I/O from Fortran,
  which is the remaining `R CMD check` NOTE. Because those reads carry no
  `iostat=` guard, a truncated or missing chain file aborts the R process
  rather than raising a catchable R error. Moving the reads into R is planned.

* `Fstat()` and the plotting functions carry the most test coverage; the
  hybrid-zone model (`HZ()`, `show.estimate.hz()`) and the simulation
  functions are exercised only lightly.


# Geneland 4.9.2 and earlier

See the package vignette and
<https://github.com/gilles-guillot/Geneland> for the earlier history.
