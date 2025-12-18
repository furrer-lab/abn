
# abn: Additive Bayesian Networks <a href="https://r-bayesian-networks.org/"><img src="man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->
<!-- WARNING: -->
<!-- The ?branch=release-x.y.y is updated automatically by the initiate_version_release workflow -->
[![status](https://joss.theoj.org/papers/1bbc43a2be86f5d3f831cedb5cf81812/status.svg)](https://joss.theoj.org/papers/10.21105/joss.06822)
[![On Label CRAN Checks](https://github.com/furrer-lab/abn/actions/workflows/onlabel_CRAN_checks.yml/badge.svg?branch=release-3.1.13)](https://github.com/furrer-lab/abn/actions/workflows/onlabel_CRAN_checks.yml)
[![Codecov](https://img.shields.io/codecov/c/github/furrer-lab/abn)](https://app.codecov.io/gh/furrer-lab/abn)
[![GitHub R package version](https://img.shields.io/github/r-package/v/furrer-lab/abn)](https://github.com/furrer-lab/abn/tags)
![cran](https://www.r-pkg.org/badges/version-ago/abn) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/abn) 
![LICENCE](https://img.shields.io/cran/l/abn)
<!-- badges: end -->

The R package `abn` is a tool for Bayesian network analysis, a form of probabilistic graphical model.
It derives a directed acyclic graph (DAG) from empirical data that describes the dependency structure between random variables.
The package provides routines for structure learning and parameter estimation of additive Bayesian network models.  

# Installation
[![Ubuntu Install](https://github.com/furrer-lab/abn/actions/workflows/Ubuntu_setup.yml/badge.svg?branch=main)](https://github.com/furrer-lab/abn/actions/workflows/Ubuntu_setup.yml)
[![Fedora Install](https://github.com/furrer-lab/abn/actions/workflows/Fedora_setup.yml/badge.svg?branch=main)](https://github.com/furrer-lab/abn/actions/workflows/Fedora_setup.yml)
[![MacOS Install](https://github.com/furrer-lab/abn/actions/workflows/Macos_setup.yml/badge.svg?branch=main)](https://github.com/furrer-lab/abn/actions/workflows/Macos_setup.yml)
[![Windows Install](https://github.com/furrer-lab/abn/actions/workflows/Windows_setup.yml/badge.svg?branch=main)](https://github.com/furrer-lab/abn/actions/workflows/Windows_setup.yml)

`abn` and its installation process relies on various software that might, or might not, be present in your system.

## Prior to installing

In order for `abn` to work correctly on your system some dependencies need to be installed.
If you are on a Linux based system (most of) these dependencies are installed automatically for you when following the [pak](https://pak.r-lib.org/)-based installation procedure described in the [Installing from GitHub](#installing-from-github-recommended) section.

For MacOS and Windows based system some more preparatory steps are required.

The following paragraphs provide detailed instructions for the most common operating systems on the steps that need to be carried out prior to installing `abn`.

<details>
<summary><b><i>Ubuntu</i></b></summary>

You presumably have R installed already, if not, open a terminal and type:

```bash
apt-get install r-base
```

_**Note:** You might need to prepend `sudo ` to this command._

All you need for the installation is to have the R-package [pak](https://pak.r-lib.org/) installed.
`pak` is installed like any other R-package, however, it relies on `curl` being present on your system, so we make sure it is there:

  ```bash
  apt-get install libcurl4-openssl-dev
  ```
  Now, to install `pak` we start an R session and write:

  ```R
  install.packages('pak', repos=c(CRAN="https://cran.r-project.org"))
  ```

  With that you should be ready to [install `abn` from GitHub](#installing-from-github-recommended).

</details>

<details>
<summary><b><i>Fedora</i></b></summary>

  You presumably have R installed already, if not, open a terminal and type:
 
  ```bash
  dnf install R
  ```

  _**Note:** You might need to prepend `sudo ` to this command._

  For the installation you need to have the R-package [pak](https://pak.r-lib.org/) installed.
  `pak` is installed like any other R-package, however, it relies on `curl` being installed on your system, so we make sure it is there:

  ```bash
  dnf install libcurl-devel 
  ```
  Now, to install `pak` we start an R session and write:

  ```R
  install.packages('pak', repos=c(CRAN="https://cran.r-project.org"))
  ```

  There is one more thing we need to do before we can install `abn`:

  **Install JAGS from source**

  [JAGS](https://mcmc-jags.sourceforge.io/), _Just Another Gibbs Sampler_, is a program for analyzing Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation. [rjags](https://cran.r-project.org/package=rjags) is R's interface to the `JAGS` library.
  `JAGS` is required in some simulations `abn` can perform.

  The steps needed to install `JAGS 4.3.2` are:
  
  ```bash
  wget -O /tmp/jags.tar.gz https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.2.tar.gz/download
  cd /tmp
  tar -xf jags.tar.gz
  cd /tmp/JAGS-4.3.2 
  ./configure --libdir=/usr/local/lib64
  make
  sudo make install
  ```
  
  _**Note:**_
  _If you are on a 64bit system (you likely are) mind the `--libdir=/usr/local/lib64` argument when launching `./configure`.)_
  _Omitting this argument will lead to `rjags` "not seeing" `jags`._
  
  On Fedora `rjags` might need some special configuration for it to link properly to the `JAGS` library.
  Also, it might be needed to add the path to the `JAGS` library to the linker path (see [rjags INSTALL file](https://github.com/cran/rjags/blob/master/INSTALL) for further details).
  
  In order to add the `JAGS` library to the linker path, run the following commands:
  
  ```bash
  sudo echo "/usr/local/lib64" > /etc/ld.so.conf.d/jags.conf
  sudo /sbin/ldconfig
  ``` 

  _**Note:**_
  _These commands might not be needed, you might first try to install the R-package `rjags` and only run them if you encounter a `configure: error: Runtime link error`._

  With that you should be ready to [install `abn` from GitHub](#installing-from-github-recommended).
  
</details>

<details>
<summary><b><i>MacOS</i></b></summary>

  Most likely you have R installed already but if not run:
 
  ```bash
  brew install R
  ```

  For the installation you need to have the R-package [pak](https://pak.r-lib.org/) installed.
  `pak` is installed like any other R-package, we start an R session and write:

  ```R
  install.packages('pak', repos=c(CRAN="https://cran.r-project.org"))
  ```

  We will install the system dependencies with [Homebrew](https://brew.sh/).
  Head over to their site to see the installation process or simply open a terminal and run:

  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

  To correctly link to installed libraries and to build them, we need `pkg-config` and `automake`:
 
  ```bash
  brew install pkg-config
  brew install automake  # needed to run autoconf
  ```
  
  We will use `wget` to download `JAGS` later, as well as, the development headers `openssl`:
 
  ```bash
  brew install wget
  brew install openssl@1.1
  ```

  <ins>**Dependencies**</ins>
  
  On MacOS we need to install some system dependencies separately:
  
  
  - **GSL**
  
    [GSL](https://www.gnu.org/software/gsl/), the _GNU Scientific Library_, is a numerical library for C/C++.
    It is required to compile `abn`'s C/C++ code.
  
    With Homebrew you can install the `GSL` binaries directly:
    
    ```
    brew install gsl
    ```
  
  - **JAGS & rjags**
    
    [JAGS](https://mcmc-jags.sourceforge.io/), _Just Another Gibbs Sampler_, is a program for analyzing Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation. [rjags](https://cran.r-project.org/package=rjags) is R's interface to the `JAGS` library.
    `JAGS` is required in some simulations `abn` can perform.
  
    With Homebrew you can install the `JAGS` binaries directly:
    
    ```
    brew install jags
    ```
    
    And now to install `rjags`, open an R session and type:
    
    ```R
    install.packages("rjags", type="source", repos=c(CRAN="https://cran.r-project.org"))
    library("rjags")
    ```
  
  - **INLA**
  
    [INLA](https://www.r-inla.org/) is an R package that is not hosted on CRAN and thus needs to be installed separately.
    `abn` uses `INLA` to fit some models. 
    
    `INLA` relies on various other R packages and C/C++ libraries.
    It thus needs some additional installation steps:
    
    ```bash
    brew install udunits 
    brew install gdal  # installs also geos as dependency
    brew install proj
    ```
  
    Now, to install `INLA` itself, simply start an R session and run:
    
    ```R
    install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
    ```
    
    If you run into trouble, please see also [INLA's installation instructions](https://www.r-inla.org/download-install) for further details.

</details>

<details>
<summary><b><i>Windows</i></b></summary>

  For the installation you need to have the R-package [pak](https://pak.r-lib.org/) installed.
  `pak` is installed like any other R-package, we start an R session and write:

  ```R
  install.packages('pak', repos=c(CRAN="https://cran.r-project.org"))
  ```


  <ins>**Dependencies**</ins>

  On Windows we need to install some system dependencies separately:


  - **GSL**
  
    [GSL](https://www.gnu.org/software/gsl/), the _GNU Scientific Library_, is a numerical library for C/C++.
    It is required to compile `abn`'s C/C++ code.
  
    In Windows `GSL` is available a.o. through [cygwin](https://cygwin.com/index.html), which has a straight forward installation process.
    Either head over to the website, download and install the `setup-x86_64.exe` file or use PowerShell:
    
    ```powershell
    Import-Module bitstransfer
    New-Item -ItemType Directory -Force -Path "C:\Program Files\cygwin"
    start-bitstransfer -source https://cygwin.com/setup-x86_64.exe "C:\Program Files\cygwin\setup-x86_64.exe"
    Start-Process -Wait -FilePath "C:\Program Files\cygwin\setup-x86_64.exe" -ArgumentList "/S" -PassThru
    ```
  
  - **JAGS & rjags**
    
    [JAGS](https://mcmc-jags.sourceforge.io/), _Just Another Gibbs Sampler_, is a program for analyzing Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation. [rjags](https://cran.r-project.org/package=rjags) is R's interface to the `JAGS` library.
    `JAGS` is required in some simulations `abn` can perform.
  
    
    You can either head over to the [JAGS download page](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/), download and execute the installable, or use PowerShell.
    The following instructions will download and install `JAGS 4.3.1` in PowerShell:
    
    ```powershell
    Import-Module bitstransfer
    New-Item -ItemType Directory -Force -Path "C:\Program Files\JAGS\JAGS-4.3.1"
    start-bitstransfer -source https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/JAGS-4.3.1.exe/download  "C:\Program Files\JAGS\JAGS-4.3.1\JAGS-4.3.1.exe"
    Start-Process -Wait -FilePath "C:\Program Files\JAGS\JAGS-4.3.1\JAGS-4.3.1.exe" -ArgumentList "/S" -PassThru
    ```
    
    In order to make sure `rjags` finds `JAGS` we set the environment variable `JAGS_HOME` before installing `rjags`.
    To do so, open your R session and type:
    
    ```R
    Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
    install.packages("rjags", repos=c(CRAN="https://cran.r-project.org"))
    library("rjags")
    ```
  
  - **INLA**
  
    [INLA](https://www.r-inla.org/) is an R package that is not hosted on CRAN and thus needs to be installed separately.
    `abn` uses `INLA` to fit some models. 
    
    The installation is straight forward, simply start an R session and run:
    
    ```R
    install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
    ```
    
    If you run into trouble, please see also [INLA's installation instructions](https://www.r-inla.org/download-install) for further details.
  
</details>

_Click on your operating system to see the specific installation instructions_

## R version support

Officially supported is R version >= 4.4

## Installing from GitHub (recommended)

From GitHub you can install any version and/or state of the `abn` repository you want.
We recommend to not directly install `main`, but to a specific version.
Head over to our [version list](https://github.com/furrer-lab/abn/releases) to see which one is the latest version.
Here we assume the version is `3.1.2`.

We use [pak](https://pak.r-lib.org/) for the installation process.
If you followed the [Prior to installing](#prior-to-installing) section `pak` should already be installed.
<details><summary>If not, install it first.</summary> Open an R session and type:

```R
install.packages('pak', repos=c(CRAN="https://cran.r-project.org"))
```
</details>

To install `abn` run in your R session:

```R
pak::repo_add(INLA = "https://inla.r-inla-download.org/R/stable/")
pak::pkg_install("furrer-lab/abn@3.1.2", dependencies=TRUE)
```
_**Note:** The first command can be skipped on MacOS or Windows._

## Installing from CRAN

> [!NOTE]
> When installing from CRAN you might not get the latest version of `abn`.
> If you want the latest version follow the instructions from [Installing from GitHub](#installing-from-github-recommended).

In order to install the `abn` version on CRAN, open an R session and type:

```R
pak::repo_add(INLA = "https://inla.r-inla-download.org/R/stable/")
pak::pkg_install("abn", dependencies=TRUE)
```
_**Note:** The first command can be skipped on MacOS or Windows._

`abn` has several dependencies that are not available on CRAN.
This is why we rely on [pak](https://pak.r-lib.org/) for the installation and the [Prior to installing](#prior-to-installing) section should be followed through before installing `abn` from CRAN. [^1]

[^1]: The `abn` package includes certain features, such as multiprocessing and integration with the `INLA` package, which are limited or available only on specific CRAN flavors. 
While it is possible to relax the testing process by, e.g., excluding tests of these functionalities, we believe that rigorous testing is important for reliable software development, especially for a package like `abn` that includes complex functionalities. 
We have implemented a rigorous testing framework similar to CRAN's to validate these functionalities in our development process. 
Our aim is to maximize the reliability of the `abn` package under various conditions, and we are dedicated to providing a robust and reliable package. 
We appreciate your understanding as we work towards making `abn` available on CRAN soon.

## Installing from source

It is also possible to clone this repository and install `abn` from source.

> [!NOTE]
> Also in this case you need to first prepare your system by following the [Prior to installing](#prior-to-installing) section.

Installing from source is done with the following steps:

1. Clone the repository and go to the root directory of the repo:

   ```bash
   git clone https://github.com/furrer-lab/abn
   cd abn
   ```

2. Deactivate `abn`'s development environment (a [renv](https://rstudio.github.io/renv/articles/renv.html) virtual environment):

   ```R
   renv::deactivate()
   ```

3. Build and install the local content with dependencies:

   ```R
   pak::repo_add(INLA = "https://inla.r-inla-download.org/R/stable/")
   pak::local_install(dependencies=TRUE)
   ```
  _**Note:** The first command can be skipped on MacOS or Windows._

# Quickstart

Explore the basics of data analysis using additive Bayesian networks with the [`abn` package](https://CRAN.R-project.org/package=abn) through our [simple example](vignettes/quick_start_example.Rmd).
The datasets required for these examples are included within the [`abn` package](https://CRAN.R-project.org/package=abn).

For a deeper understanding, refer to the manual pages on the [`abn` homepage](https://r-bayesian-networks.org/), which include numerous examples. 
Key pages to visit are `fitAbn()`, `buildScoreCache()`, `mostProbable()`, and `searchHillClimber()`.
Also, see the [examples](#examples) below for a quick overview of the package's capabilities.

# Features

The R package `abn` provides routines for determining optimal additive Bayesian network models for a given data set. 
The core functionality is concerned with model selection - determining the most likely model of data from interdependent variables. 
The model selection process can incorporate expert knowledge by specifying structural constraints, such as which arcs are banned or retained.

The general workflow with `abn` follows a three-step process:

1. **Determine the model search space**: The function `buildScoreCache()` builds a cache of pre-computed scores for each possible DAG.
For this, it's required to specify the data types of the variables in the data set and the structural constraints of the model (e.g. which arcs are banned or retained and the maximum number of parents per node).

2. **Structure learning**: `abn` offers different structure learning algorithms:
    - The exact structure learning algorithm from [Koivisto and Sood (2004)](https://www.jmlr.org/papers/volume5/koivisto04a/koivisto04a.pdf) is implemented in `C` and can be called with the function `mostProbable()`, which finds the most probable DAG for a given data set.
    The function `searchHeuristic()` provides a set of heuristic search algorithms. These include the hill-climber, tabu search, and simulated annealing algorithms implemented in `R`.
    `searchHillClimber()` searches for high-scoring DAGs using a random re-start greedy hill-climber heuristic search and is implemented in `C`. It slightly deviates from the method initially presented by [Heckerman et al. 1995](https://doi.org/10.1023/A:1022623210503) (for details consult the respective help page `?abn::searchHillClimber()`).

3. **Parameter estimation**: The function `fitAbn()` estimates the model's parameters based on the DAG from the previous step.

`abn` allows for two different model formulations, specified with the argument `method`:

- `method = "mle"` fits a model under the frequentist paradigm using information-theoretic criteria to select the best model.

- `method = "bayes"` estimates the posterior distribution of the model parameters based on two Laplace approximation methods, that is, a method for Bayesian inference and an alternative to Markov Chain Monte Carlo (MCMC): A standard Laplace approximation is implemented in the `abn` source code but switches in specific cases (see help page `?fitAbn`) to the Integrated Nested Laplace Approximation from the [INLA package](https://www.r-inla.org) requiring the installation thereof.

To generate new observations from a fitted ABN model, the function `simulateAbn()` simulates data based on the DAG and the estimated parameters from the previous step. `simulateAbn()` is available for both `method = "mle"` and `method = "bayes"` and requires the installation of the [JAGS package](https://mcmc-jags.sourceforge.io). 

## Supported Data types

The `abn` package supports the following distributions for the variables in the network:

- Gaussian distribution for continuous variables.

- Binomial distribution for binary variables.

- Poisson distribution for variables with count data.

- Multinomial distribution for categorical variables (only available with `method = "mle"`).

Unlike other packages, `abn` does not restrict the combination of parent-child distributions.

## Multilevel Models for Grouped Data Structures

The analysis of "hierarchical" or "grouped" data, in which observations are nested within higher-level units, requires statistical models with parameters that vary across groups (e.g. mixed-effect models).

`abn` allows to control for one-layer clustering, where observations are grouped into a single layer of clusters that are themself assumed to be independent, but observations within the clusters may be correlated (e.g. students nested within schools, measurements over time for each patient, etc).
The argument `group.var` specifies the discrete variable that defines the group structure. The model is then fitted separately for each group, and the results are combined. 

For example, studying student test scores across different schools, a varying intercept model would allow for the possibility that average test scores (the intercept) might be higher in one school than another due to factors specific to each school. This can be modeled in `abn` by setting the argument `group.var` to the variable containing the school names. The model is then fitted as a varying intercept model, where the intercept is allowed to vary across schools, but the slope is assumed to be the same for all schools.

Under the frequentist paradigm (`method = "mle"`), `abn` relies on the `lme4` package to fit generalized linear mixed models (GLMMs) for Binomial, Poisson, and Gaussian distributed variables. For multinomial distributed variables, `abn` fits a multinomial baseline category logit model with random effects using the `mclogit` package. Currently, only one-layer clustering is supported (e.g., for `method = "mle"`, this corresponds to a random intercept model).

With a Bayesian approach (`method = "bayes"`), `abn` relies on its own implementation of the Laplace approximation and the package `INLA` to fit a single-level hierarchical model for Binomial, Poisson, and Gaussian distributed variables. Multinomial distributed variables in general (see Section [Supported Data Types](#supported-data-types)) are not yet implemented with `method = "bayes"`.

# Basic Background

Bayesian network modeling is a data analysis technique ideally suited to messy, highly correlated and complex datasets. 
This methodology is rather distinct from other forms of statistical modeling in that its focus is on structure discovery—determining an optimal graphical model that describes the interrelationships in the underlying processes that generated the data. 
It is a **multivariate** technique and can be used for one or many dependent variables. 
This is a data-driven approach, as opposed to relying only on subjective expert opinion to determine how variables of interest are interrelated (for example, structural equation modeling). 

[Below](#examples) and on the [package's website](https://r-bayesian-networks.org/), we provide some [cookbook](#examples)-type examples of how to perform Bayesian network **structure discovery** analyses with observational data. 
The particular type of Bayesian network models considered here are **additive Bayesian networks**. 
These are rather different, mathematically speaking, from the standard form of Bayesian network models (for binary or categorical data) presented in the academic literature, which typically use an analytically elegant but arguably interpretation-wise opaque contingency table parametrization. 
An additive Bayesian network model is simply a **multidimensional regression model**, e.g., directly analogous to generalized linear modeling but with all variables potentially dependent. 

An example can be found in the [American Journal of Epidemiology](https://doi.org/10.1093/aje/kws183), where this approach was used to investigate risk factors for child diarrhea. 
A special issue of **Preventive Veterinary Medicine** on graphical modeling features several articles that use [abn](https://CRAN.R-project.org/package=abn) to fit epidemiological data (e.g., [Ludwig et al., 2013](https://doi.org/10.1016/j.prevetmed.2013.02.005)). 
Introductions to this methodology can be found in [Emerging Themes in Epidemiology](https://link.springer.com/journal/12982) and in [Computers in Biology and Medicine](https://doi.org/10.1016/j.compbiomed.2022.105740) where it is compared to other approaches.

## What is an additive Bayesian network?

Additive Bayesian network (ABN) models are statistical models that use the principles of Bayesian statistics and graph theory. 
They provide a framework for representing data with multiple variables, known as multivariate data.

ABN models are a graphical representation of (Bayesian) multivariate regression. 
This form of statistical analysis enables the prediction of multiple outcomes from a given set of predictors while simultaneously accounting for the relationships between these outcomes.

In other words, additive Bayesian network models extend the concept of generalized linear models (GLMs), which are typically used to predict a single outcome, to scenarios with multiple dependent variables. 
This makes them a powerful tool for understanding complex, multivariate datasets.

## The term Bayesian network is interpreted differently across various fields.

Bayesian network models often involve binary nodes, arguably the most frequently used type of Bayesian network. 
These models typically use a contingency table instead of an additive parameter formulation. 
This approach allows for mathematical elegance and enables key metrics like model goodness of fit and marginal posterior parameters to be estimated analytically (i.e., from a formula) rather than numerically (an approximation). 
However, this parametrization may not be parsimonious, and the interpretation of the model parameters is less straightforward than the usual Generalized Linear Model (GLM) type models, which are prevalent across all scientific disciplines.

While this is a crucial practical distinction, it’s a relatively low-level technical one, as the primary aspect of BN modeling is that it’s a form of graphical modeling – a model of the data’s joint probability distribution. 
This joint – multidimensional – aspect makes this methodology highly attractive for complex data analysis and sets it apart from more standard regression techniques, such as GLMs, GLMMs, etc., which are only one-dimensional as they assume all covariates are independent. 
While this assumption is entirely reasonable in a classical experimental design scenario, it’s unrealistic for many observational studies in fields like medicine, veterinary science, ecology, and biology.

# Examples

- [Example 1: Basic usage](#example-1-basic-usage)
- [Example 2: Restrict model search space](#example-2-restrict-model-search-space)
- [Example 3: Grouped Data Structures](#example-3-grouped-data-structures)
- [Example 4: Using INLA vs internal Laplace approximation](#example-4-using-inla-vs-internal-laplace-approximation)

## Example 1: Basic Usage

This is a basic example which shows the basic workflow:

``` r
library(abn)

# Built-in toy dataset with two Gaussian variables G1 and G2, two Binomial variables B1 and B2, and one multinomial variable C
str(g2b2c_data)

# Define the distributions of the variables
dists <- list(G1 = "gaussian",
              B1 = "binomial",
              B2 = "binomial",
              C = "multinomial",
              G2 = "gaussian")


# Build the score cache
cacheMLE <- buildScoreCache(data.df = g2b2c_data,
                         data.dists = dists,
                         method = "mle",
                         max.parents = 2)

# Find the most probable DAG
dagMP <- mostProbable(score.cache = cacheMLE)

# Print the most probable DAG
print(dagMP)

# Plot the most probable DAG
plot(dagMP)

# Fit the most probable DAG
myfit <- fitAbn(object = dagMP,
                method = "mle")

# Print the fitted DAG
print(myfit)
```

## Example 2: Restrict Model Search Space

Based on [example 1](#example-1-basic-usage), we may know that the arc G1->G2 is not possible and that the arc from C -> G2 must be present.
This "expert knowledge" can be included in the model by banning the arc from G1 to G2 and retaining the arc from C to G2.

The retain and ban matrices are specified as an adjacency matrix of 0 and 1 entries, where 1 indicates that the arc is banned or retained, respectively. 
Row and column names must match the variable names in the data set. 
The corresponding column is a parent of the variable in the row.
Each column represents the parents, and the row is the child. For example, the first row of the ban matrix indicates that G1 is banned as a parent of G2.

Further, we can restrict the maximum number of parents per node to 2.

```r

# Ban the edge G1 -> G2
banmat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
banmat[1, 5] <- 1

# retain always the edge C -> G2
retainmat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
retainmat[5, 4] <- 1

# Limit the maximum number of parents to 2
max.par <- 2

# Build the score cache
cacheMLE_small <- buildScoreCache(data.df = g2b2c_data,
                            data.dists = dists,
                            method = "mle",
                            dag.banned = banmat,
                            dag.retained = retainmat,
                            max.parents = max.par)
print(paste("Without restrictions from example 1: ", nrow(cacheMLE$node.defn)))
print(paste("With restrictions as in example 2: ", nrow(cacheMLE_small$node.defn)))

```

## Example 3: Grouped Data Structures

Depending on the data structure, we may want to control for one-layer clustering, where observations are grouped into a single layer of clusters that are themselves assumed to be independent, but observations within the clusters may be correlated (e.g., students nested within schools, measurements over time for each patient, etc.).

Currently, `abn` supports only one layer clustering. 

```r

# Built-in toy data set
str(g2pbcgrp)

# Define the distributions of the variables
dists <- list(G1 = "gaussian",
              P = "poisson",
              B = "binomial",
              C = "multinomial",
              G2 = "gaussian") # group is not among the list of variable distributions

# Ban arcs such that C has only B and P as parents
ban.mat <- matrix(0, nrow = 5, ncol = 5, dimnames = list(names(dists), names(dists)))
ban.mat[4, 1] <- 1
ban.mat[4, 4] <- 1
ban.mat[4, 5] <- 1

# Build the score cache
cache <- buildScoreCache(data.df = g2pbcgrp,
                         data.dists = dists,
                         group.var = "group",
                         dag.banned = ban.mat,
                         method = "mle",
                         max.parents = 2)

# Find the most probable DAG
dag <- mostProbable(score.cache = cache)

# Plot the most probable DAG
plot(dag)

# Fit the most probable DAG
fit <- fitAbn(object = dag,
              method = "mle")

# Plot the fitted DAG
plot(fit)

# Print the fitted DAG
print(fit)

```


## Example 4: Using INLA vs internal Laplace approximation

Under a Bayesian approach, `abn` automatically switches to the Integrated Nested Laplace Approximation from the [INLA package](https://www.r-inla.org) if the internal Laplace approximation fails to converge. 
However, we can also force the use of INLA by setting the argument `control=list(max.mode.error=100)`.

The following example shows that the results are very similar. It also shows how to constrain arcs as formula objects and how to specify different parent limits for each node separately.

``` r
library(abn)

# Subset of the build-in dataset, see  ?ex0.dag.data
mydat <- ex0.dag.data[,c("b1","b2","g1","g2","b3","g3")] ## take a subset of cols

# setup distribution list for each node
mydists <- list(b1="binomial", b2="binomial", g1="gaussian",
                g2="gaussian", b3="binomial", g3="gaussian")

# Structural constraints
## ban arc from b2 to b1
## always retain arc from g2 to g1
## parent limits - can be specified for each node separately
max.par <- list("b1"=2, "b2"=2, "g1"=2, "g2"=2, "b3"=2, "g3"=2)

# now build the cache of pre-computed scores according to the structural constraints
res.c <- buildScoreCache(data.df=mydat, data.dists=mydists,
                         dag.banned= ~b1|b2, 
                         dag.retained= ~g1|g2, 
                         max.parents=max.par)


# repeat but using R-INLA. The mlik's should be virtually identical.
if(requireNamespace("INLA", quietly = TRUE)){
  res.inla <- buildScoreCache(data.df=mydat, data.dists=mydists,
                              dag.banned= ~b1|b2, # ban arc from b2 to b1
                              dag.retained= ~g1|g2, # always retain arc from g2 to g1
                              max.parents=max.par,
                              control=list(max.mode.error=100)) # force using of INLA
  
  ## comparison - very similar
  difference <- res.c$mlik - res.inla$mlik
  summary(difference)
}
```

# Contributing

We greatly appreciate contributions from the community and are excited to welcome you to the development process of the `abn` package. 
Here are some guidelines to help you get started:

1. **Seeking Support:** 
If you need help with using the `abn` package, you can seek support by creating a new issue on our GitHub repository. 
Please describe your problem in detail and include a minimal reproducible example if possible.

2. **Reporting Issues or Problems:** 
If you encounter any issues or problems with the software, please report them by creating a new issue on our GitHub repository.
When reporting an issue, try to include as much detail as possible, including steps to reproduce the issue, your operating system and R version, and any error messages you received.

3. **Software Contributions:**
We encourage contributions directly via pull requests on our GitHub repository.
Before starting your work, please first create an issue describing the contribution you wish to make. 
This allows us to discuss and agree on the best way to integrate your contribution into the package.

By participating in this project, you agree to abide by our [code of conduct](#code-of-conduct). 
We are committed to making participation in this project a respectful and harassment-free experience for everyone.

# Citation

If you use `abn` in your research, please cite it as follows:

``` r
> citation("abn")
To cite the software implementation of the R package 'abn' use:

  Delucchi M, Furrer R, Kratzer G, Lewis F, Liechti J, Pittavino M, Cherneva K (2024). _abn: Modelling Multivariate Data with Additive Bayesian Networks_. R package version 3.1.3, <https://CRAN.R-project.org/package=abn>.

To cite the methodology of the R package 'abn' use:

  Kratzer G, Lewis F, Comin A, Pittavino M, Furrer R (2023). “Additive Bayesian Network Modeling with the R Package abn.” _Journal of Statistical Software_, *105*(8), 1-41. doi:10.18637/jss.v105.i08 <https://doi.org/10.18637/jss.v105.i08>.

To cite the application of mixed-effects ABN use:

  Delucchi M, Liechti J, Spinner G, Furrer R (2024). “abn: Additive Bayesian Networks.” _Journal of Open Source Software_, *9*(101), 6822. R package version 3.1.3, <https://doi.org/10.21105/joss.06822>.

To cite an example of a typical ABN analysis use:

  Kratzer, G., Lewis, F.I., Willi, B., Meli, M.L., Boretti, F.S., Hofmann-Lehmann, R., Torgerson, P., Furrer, R. and Hartnack, S. (2020). Bayesian Network Modeling Applied to Feline Calicivirus Infection Among Cats in Switzerland. Frontiers
  in Veterinary Science, 7, 73
```

# License

The `abn` package is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Code of Conduct

Please note that the `abn` project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

# Applications

The [abn website](https://r-bayesian-networks.org/) provides a comprehensive set of documented case studies, numerical accuracy/quality assurance exercises, and additional documentation.

## Technical articles

- Delucchi et al. (2024): [Additive Bayesian Networks. Journal of Open Source Software, 9(101), 6822](https://doi.org/10.21105/joss.06822)

- Kratzer et al. (2023): [Additive Bayesian Network Modeling with the R Package abn](https://doi.org/10.18637/jss.v105.i08)

- Kratzer et al. (2020) [Bayesian Networks modeling applied to Feline Calicivirus infection among cats in Switzerland](https://doi.org/10.3389/fvets.2020.00073)

- Kratzer et al. (2018): [Comparison between Suitable Priors for Additive Bayesian Networks](https://arxiv.org/pdf/1809.06636)

- Koivisto et al. (2004): [Exact Bayesian structure discovery in Bayesian networks](https://www.jmlr.org/papers/volume5/koivisto04a/koivisto04a.pdf)

- Friedman et al. (2003): [Being Bayesian about network structure. A Bayesian approach to structure discovery in Bayesian networks](https://doi.org/10.1023/A:1020249912095)

- Friedman et al. (1999): [Data analysis with Bayesian networks: A bootstrap approach](https://arxiv.org/abs/1301.6695)

- Heckerman et al. (1995): [Learning Bayesian Networks – The Combination of Knowledge And Statistical-Data](http://maxchickering.com/publications/ml95.pdf)

## Application articles

- Delucchi et al. (2022): [Bayesian network analysis reveals the interplay of intracranial aneurysm rupture risk factors](https://doi.org/10.1016/j.compbiomed.2022.105740)

- Guinat et al. (2020) [Biosecurity risk factors for highly pathogenic avian influenza (H5N8) virus infection in duck farms, France](https://onlinelibrary.wiley.com/doi/10.1111/tbed.13672)

- Hartnack et al. (2019) [Additive Bayesian networks for antimicrobial resistance and potential risk factors in non-typhoidal Salmonella isolates from layer hens in Uganda](https://doi.org/10.1186/s12917-019-1965-y)

- Ruchti et al. (2019): [Progression and risk factors of pododermatitis in part-time group housed rabbit does in Switzerland](https://doi.org/10.1016/j.prevetmed.2019.01.013)

- Comin et al. (2019) [Revealing the structure of the associations between housing system, facilities, management and welfare of commercial laying hens using Additive Bayesian Networks](https://doi.org/10.1016/j.prevetmed.2019.01.004)

- Ruchti et al. (2018): [Pododermatitis in group housed rabbit does in Switzerland – prevalence, severity and risk factors](https://doi.org/10.1016/j.prevetmed.2018.06.011)

- Pittavino et al. (2017): [Comparison between generalised linear modelling and additive Bayesian network; identification of factors associated with the incidence of antibodies against Leptospira interrogans sv Pomona in meat workers in New Zealand](https://doi.org/10.1016/j.actatropica.2017.04.034)

- Hartnack et al. (2017): [Attitudes of Austrian veterinarians towards euthanasia in small animal practice: impacts of age and gender on views on euthanasia](https://doi.org/10.1186/s12917-016-0649-0)

- Lewis et al. (2012): [Revealing the Complexity of Health Determinants in Resource-poor Settings](https://doi.org/10.1093/aje/kws183)

- Lewis et al. (2011): [Structure discovery in Bayesian networks: An analytical tool for analysing complex animal health data](https://doi.org/10.1016/j.prevetmed.2011.02.003)

## Workshops

### Causality:

- 4 December 2018, Beate Sick & Gilles Kratzer of the [1st Causality workshop](https://bsick.github.io/causality_workshop/) **talk**, Bayesian Networks meet Observational data. (UZH, Switzerland)

### ABN modeling

- 07 July 2021, **workshop** at the [UseR!](https://user2021.r-project.org/) Conference on [Additive Bayesian Networks Modeling](https://gilleskratzer.github.io/ABN-UseR-2021/). (Online)

- 29 March 2019, **workshop** at the SVEPM conference on [Multivariate analysis using Additive Bayesian Networks](https://gilleskratzer.github.io/SVEPM2019/). (Utrecht, Netherland)

## Presentations

- 4 October 2018, **talk** in Nutricia (Danone). Multivariable analysis: variable and model selection in system epidemiology. (Utrecht, Netherland)

- 30 May 2018. [Brown Bag Seminar](https://tensorchiefs.github.io/bbs/) in ZHAW. **Presentation**: Bayesian Networks Learning in a Nutshell. (Winterthur, Switzerland)
