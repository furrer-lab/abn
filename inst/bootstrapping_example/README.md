# Overview

This directory illustrates how to perform a bootstrap analysis. An older version of the vignette referred to this directory quite extensively. The vignette is now in the directory `inst/old_vignette`.

Note that most of the R or script files need manual adaption depending on the system or setup. 
The files are for `abn` version >= 2.2.

To start, launch `bootstrapping.R` in Rstudio. To launch many batches, launch the bash script `bootstrapping.bash`.

## Files in the directory

* `README.md` this file
* `pigs_post_params.Rds` contains a output of a posterior analysis.
* `bootstrapping.R`, `bootstrapping.bash` R script for the main bootstrapping loop and bash to launch several of these loops (to make use of several cores, for example).
* `pigs_model.bug` JAGS model

## The following files are constructed and then deleted

* `pigs_post_params.R` temporary data for JAGS algorithm
* `script_*`, `ini_*`, `out_*chain1.txt`, `out_*index.txt`, 

## Note

Older versions of the `abn` pakage contain a few more files in the directory. Please see the CRAN abn-archive, e.g.,
https://cran.r-project.org/src/contrib/Archive/abn/abn_1.0.tar.gz for a more detailed exposition.





