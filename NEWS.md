
# abn 3.1.9
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.9)

> 2025-06-26


# abn 3.1.8
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.8)

> 2025-06-24


# abn 3.1.7
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.7)

> 2025-06-06


# abn 3.1.6
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.6)

> 2025-06-05


# abn 3.1.5
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.5)

> 2025-06-03


# abn 3.1.4
[diff](https://github.com/furrer-lab/abn/compare/3.1.3...3.1.4)

> 2025-05-29


# abn 3.1.3
[diff](https://github.com/furrer-lab/abn/compare/3.1.2...3.1.3)

> 2024-09-18

* JOSS paper final edits ([#152](https://github.com/furrer-lab/abn/issues/152))

# abn 3.1.2
[diff](https://github.com/furrer-lab/abn/compare/3.1.1...3.1.2)

> 2024-08-07

### Fix

* Simplified installation process ([#143](https://github.com/furrer-lab/abn/issues/143), [#142](https://github.com/furrer-lab/abn/issues/142), [#141](https://github.com/furrer-lab/abn/issues/141), [#125](https://github.com/furrer-lab/abn/issues/125), [#124](https://github.com/furrer-lab/abn/issues/124), [#123](https://github.com/furrer-lab/abn/issues/123), [#117](https://github.com/furrer-lab/abn/issues/117), [#114](https://github.com/furrer-lab/abn/issues/114))
* Minor fixes in JOSS paper ([#147](https://github.com/furrer-lab/abn/issues/147), [#129](https://github.com/furrer-lab/abn/issues/129))
* Some fixes and improvements in automation pipelines ([#140](https://github.com/furrer-lab/abn/issues/140), [#138](https://github.com/furrer-lab/abn/issues/138), [#132](https://github.com/furrer-lab/abn/issues/132), [#127](https://github.com/furrer-lab/abn/issues/127), [#111](https://github.com/furrer-lab/abn/issues/111))

# abn 3.1.1
[diff](https://github.com/furrer-lab/abn/compare/3.1.0...3.1.1)

> 2024-05-22

* Fixed URL in README
* Added configuration for lintr
* Removed non-standard files from being included in the package archive.

# abn 3.1.0
[diff](https://github.com/furrer-lab/abn/compare/3.0.7...3.1.0)

> 2024-05-20

* automated version release process
* JOSS paper content
* some checks should not run on release- branches
* various typos

# abn v3.0.7

2024-03-27 Matteo Delucchi

  * Reduced package size by compressions.
  * Speed up the checks on CRAN by reducing the number of examples to be run.
  * Fixed the checks on CRAN by resetting user options within functions and examples.
  * Fixed failing tests on CRAN due to the use of `INLA` package.
  * Reduced the number of empty tests by gracefully handling the specific conditions.
  * Improved documentation, fixed typos and added return values for each function.
  * Added Jonas Liechti as a contributor.

# abn 3.0.6

2024-03-21 Matteo Delucchi

  * New website
  * added JOSS paper
  * Changed `ChangeLog` into `NEWS.md`
  * Revised CITATION file
  * Updated README
  * added logo
  * Revised class names and methods
  * Added functional test for `getmarginals()`
  * Check for validity and sensibility of the `compute.fixed` argument.
  * Check for validity and sensibility of the `check` argument in `buildscoreCache()`.
  * Check codecoverage
  * harmonised plotting across different classes
  * Created new vignettes
  * relaxed tests to allow for complete and partial scientific notation in string outputs from BUGS files.
  * avoid examples and tests on CRAN that require INLA
  * Restructured the private developement repository to have only one public repo.

# abn 3.0.5

2024-01-18 Matteo Delucchi

  * Parallelization of `buildScoreCache(method='bayes')` and `fitAbn(method='bayes')` with `foreach()`.
  * Controlling the number of cores used by INLA. Before it was always using all available cores.
  * Specifying the type of parallelisation in `build.control()` and `fit.control()`.
  * Fixed bug in test of `simulateAbn()` that was failing on CRAN.
  * Increase package availability by skipping tests on CRAN that require INLA. This also speeds up the CRAN checks.

# abn 3.0.4

2023-11-03 Matteo Delucchi

  * Fixed format issues in source code that caused significant warnings with GCC and clang.
  * Temporary workaround to install the correct versions of INLA dependencies. This will be hopefully fixed in the next release of INLA.

# abn 3.0.3

2023-11-03 Matteo Delucchi

  * Temporarily fixed tests on CRAN Fedora flavour.

# abn 3.0.2

2023-10-23 Matteo Delucchi

  * Temporarily fixed tests on CRAN for multithreading.
  * Replaced README.md file.
  * Fixed bug in plotAbn() causing non-sensical edge labels.
  * New CoC and Licence files.
  * Updated examples in documentation.

# abn 3.0.1

2023-09-29 Matteo Delucchi

  * Patched memory access error in source code.
  * Fixed reverse dependency to `mcmcabn()`.

# abn 3.0.0

2023-09-01 Matteo Delucchi

  * Implemented random-intercept model for `method="mle"` with argument `group.var`.
  * simulateAbn(): Fixed to run with all arguments. Enhancement: Simulates a hierarchical model if `group.var` is provided.
  * Fixed parallel computing for `method="mle"` that did not always run smoothly.
  * Beware, the structure of some objects has changed!
  * Updated documentation to Roxygen.
  * Rigorous checking of function arguments.
  * Upgraded unit test for testthat e3.
  * Implemented extensive unit tests and CI/CD pipelines.
  * Many (!) bugs fixed that mainly flew under the radar.
  * Removed deprecated functions.
  * Refactoring here and there. Especially to harmonise `buildScoreCache()` and `fitAbn()` for both methods "mle" and "bayes".

# abn 2.7-0

2022-04-02 Reinhard Furrer

  * plotAbn(): new arguments, now fully based on class `graphAM`, proper handling of lty, lwd etc.
  * options of `fitAbn()` and `buildScoreCache()` are now handled with a list passed to argument `control=`. Some of the default values have been slightly adapted.
  * Default list available via `fit.control()` and `build.control()`, including more options!
  * Significantly improved help of several functions and examples.
  * Improved display of print methods.
  * Possiblity to use single threaded with `for()` loop and with `foreach()`.
  * Due to proper method handling, the arguments `create.graph` in `searchHillClimber()` and `fitAbn(method="bayes")` has been removed.
  * Argument `verbose` properly implemented (including passing to JAGS functions).
  * New function `compareEG()` to compare two essential graphs.
  * Additional argument checking an cleaning of many error messages.


  ### Internal changes:

  * Option `error.verbose` properly propagated to C functions.
	* New option `trace` passed to "lbfgsb".
	* Less ARMA messages, using now `#define ARMA_WARN_LEVEL 1`.
	* Some GSL error now caught at R level.
	* Severe R and moderate C code cleaning.

# abn 2.6.0

2021-06-08 Gilles Kratzer, Reinhard Furrer & Kalina Cherneva

  * Bug fixes:
     . initial values for several chains in `simulateAbn()`
     . new default value for `method` in `fitAbn()`
     . new default value `method='bic'` in `mostProbable()`
     . minor bug fixes in `plotAbn()`
     . minor bug fix in `help(plotAbn)`
     . number of parent specification for `buildScoreCache( , method='mle')`
  * Renaming many function arguments for consistency
  * Reordering, purging and adding function argments
  * Parallelization of `buildScoreCache(, method='mle')`
  * Cleaning of RCpp with improved default parameters (speed up 5%)
  * Adapting authorship of package
  * Minor cleaning of `Examples` in the help files

# abn 2.5.0

2021-04-23 Gilles Kratzer & Reinhard Furrer

  * Renaming many functions to CamelCase for consistency
  * Renaming many arguments for consistency
  * Reordering, purging and adding function argments
  * Major changes in `plotAbn()`
  * Cleaning of several help files

# abn 2.4.0

2020-10-31 Gilles Kratzer & Reinhard Furrer

  * Update on citations
  * Update on ORCID number

# abn 2.3.0

2020-10-22 Gilles Kratzer & Reinhard Furrer

  * Requirement for use of R>=4.0.0

# abn 2.2.1

2020-06-13 Gilles Kratzer & Reinhard Furrer

  * Patch for S4 exports not stated in the package
  * Testing automatic (testhat.R in tests repository)
  * Updated version of the vignette

# abn 2.2

2019-11-22 Gilles Kratzer & Reinhard Furrer

  * New vignette describing new functionalities!
  * Old vignette is store in the﻿ directories inst/bootstrapping_example and inst/old_vignette
  * Dataset: FCV
  * Methods related to the class abnFit: family, logLik, AIC, BIC, nobs, coef
  * Bugs corrected in plotabn()
  * Compilation warnings corrected
  * Typos in man help files corrected
  * Compiler flags that were not compatible with the Solaris OS

# abn 2.1

2019-07-04 Gilles Kratzer & Reinhard Furrer

  * Compiler flags that were not compatible with the Solaris system
  * Note related to INLA

# abn 2.0

2019-07-01 Gilles Kratzer & Reinhard Furrer

  * Implementation of S3 methods (print, summary, plot)
  * New classes:
  		* abnDag, abnCache, abnHeuristic, abnMostprobable, abnHillClimber, abnFit
  * renaming of functions
  		* link.strength() -> linkStrength()
  		* search.heuristic() -> searchHeuristic()
  		* backward compatibility is not ensured for those functions
  * New tests implemented for internal validation and some backward compatibility
  * Dependence with Cairo removed
  * modification of search.hillclimber(): no interactive plotting possible
  * buildscorecache() becomes a parametrisable wrapper for the bayesian and the mle methods
  * fitabn() becomes a parametrisable wrapper for the bayesian and the mle methods
  * Dataset: adg
  * Multiple Bugs corrected

# abn 1.2

2018-07-23 Gilles Kratzer

  * (many) new functions
  * buildscorecache.mle()
  * fitabn.mle()
  * plotabn()
  * simulateabn()
  * expit(), logit(), odds() and or()
  * mb()
  * link.strength()
  * discretization()
  * essentialGraph()
  * compareDag()
  * infoDag()
  * simulateDag()
  * search.heuristic()
  * entropyData()
  * miData()

  * modification of existing functions:
  * mostprobable()

# abn 1.0.2

2016-11-09 Gilles Kratzer

  * Default ban and retain matrix
  * Automatic testing procedure of the package

# abn 1.0

2016-01-11 Pittavino Marta

  * New vignette uploaded.

# abn 0.86

2015-12-22 Pittavino Marta

  * INLA routine updated, ERROR and NOTE addressed.

# abn 0.85

2014-12-22 Pittavino Marta

  * Published Version v0.85 on CRAN

# abn 0.85

2014-12-09 Pittavino Marta

  * cleanup file to delete additional Makevars added

# abn 0.85

2014-12-08 Pittavino Marta

  * FOSS license --as-cran NOTE addressed, eliminated BuildVignettes NO in the DESCRIPTION file

# abn 0.84

2014-12-06 Pittavino Marta

  * v0.84 submitted to CRAN

# abn 0.84

2014-12-05 Pittavino Marta

  * dependencies --as-cran NOTE addressed, with .Rprofile
  * graphical issue with Sweave solved

# abn 0.84

2014-12-04 Pittavino Marta

  * finished to edit vignette

# abn 0.84

2014-11-23/24 Pittavino Marta

  * updating the vignette

# abn 0.84

2014-11-03/04 Pittavino Marta

  * issue with Rd2pdf and LaTeX addressed

# abn 0.84

2014-10-27 Pittavino Marta

  * installation of local version of abn to text and work with vignette

# abn 0.84

2014-10-23 Pittavino Marta

  * confirmation of the change of maintainer

# abn 0.84

2014-10-13 Pittavino Marta

  * updated Rclean for file .log and .status
  * check length characters in manual

# abn 0.84

2014-10-06 Pittavino Marta
  * compressed Data with resaveRdaFiles(,compress=c("xz"))
  * solved installing issue (Rclean update)

# abn 0.84

2014-10-03 Pittavino Marta

  * alias Rclean extended
  * eliminate old and very old version
  * comment out require(Cairo) line-code

# abn 0.84

2014-10-02 Reinhard Furrer

  * Cairo NOTE fixed, moved from Suggests to Depends
  * DESCRIPTION and NAMESPACE update
  * check if last amendments works
  * .gitignore create

# abn 0.84

2014-09-29 Pittavino Marta

  * update version
  * rename R function as the one in manual

# abn 0.83

2014-09-26 Pittavino Marta

  * update DESCRIPTION

# abn 0.83

2014-09-19 Pittavino Marta

  * update DESCRIPTION
  * update tographviz documentation
  * create ChangeLog file
