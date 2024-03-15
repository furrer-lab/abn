#' @title Dataset related to average daily growth performance and abattoir findings in pigs commercial production.
#' @usage adg
#' @description
#' The case study dataset is about growth performance and abattoir findings in pigs commercial production in a selected set of 15 Canadian farms collected in March 1987.
#' @format
#' An adapted data frame of the original dataset which  consists of 341 observations of 8 variables and a grouping variable (farm).
#' \describe{
#'   \item{AR}{presence of atrophic rhinitis.}
#'   \item{pneumS}{presence of moderate to severe pneumonia.}
#'   \item{female}{sex of the pig (1=female, 0=castrated). }
#'   \item{livdam}{presence of liver damage (parasite-induced white spots).}
#'   \item{eggs}{presence of fecal/gastrointestinal nematode eggs at time of slaughter.}
#'   \item{wormCount}{count of nematodes in small intestine at time of slaughter.}
#'   \item{age}{days elapsed from birth to slaughter (days).}
#'   \item{adg}{average daily weight gain (grams).}
#'   \item{farm}{farm ID.}
#' }
#' @details
#' When using the data to fit an additive Bayesian network,
#' the variables \code{AR}, \code{pneumS}, \code{female}, \code{livdam},
#' \code{eggs} are considered binomial, \code{wormcount} Poisson,
#' \code{age} and \code{adg} Gaussian.
#' The variable \code{farm} can be used to adjust for grouping.
#'
#' @references
#' Kratzer, G., Lewis, F.I., Comin, A., Pittavino, M. and Furrer, R. (2019). "Additive Bayesian Network Modelling with the R Package abn". arXiv preprint arXiv:1911.09006.
#' Dohoo, Ian Robert, Wayne Martin, and Henrik Stryhn. Veterinary epidemiologic research. No. V413 DOHv. Charlottetown, Canada: AVC Incorporated, 2003.
#' @keywords datasets internal
"adg"

#' @title Synthetic validation data set for use with abn library examples
#' @description
#' 300 observations simulated from a DAG with 10 binary variables, 10 Gaussian variables and 10 poisson variables.
#' @usage ex0.dag.data
#'
#' @format A data frame, binary variables are factors.
#' The relevant formulas are given below (note these do not give parameter
#' estimates just the form of the relationships, e.g. logit()=1 means a logit
#' link function and comprises of only an intercept term).
#'
#' \describe{
#'   \item{b1}{binary, logit()=1 }
#'   \item{b2}{binary, logit()=1 }
#'   \item{b3}{binary, logit()=1 }
#'   \item{b4}{binary, logit()=1 }
#'   \item{b5}{binary, logit()=1 }
#'   \item{b6}{binary, logit()=1 }
#'   \item{b7}{binary, logit()=1 }
#'   \item{b8}{binary, logit()=1 }
#'   \item{b9}{binary, logit()=1 }
#'   \item{b10}{binary, logit()=1 }
#'   \item{g1}{gaussian, identity()=1 }
#'   \item{g2}{gaussian, identity()=1 }
#'   \item{g3}{gaussian, identity()=1 }
#'   \item{g4}{gaussian, identity()=1 }
#'   \item{g5}{gaussian, identity()=1 }
#'   \item{g6}{gaussian, identity()=1 }
#'   \item{g7}{gaussian, identity()=1 }
#'   \item{g8}{gaussian, identity()=1 }
#'   \item{g9}{gaussian, identity()=1 }
#'   \item{g10}{gaussian, identity()=1 }
#'   \item{p1}{poisson, log()=1 }
#'   \item{p2}{poisson, log()=1 }
#'   \item{p3}{poisson, log()=1 }
#'   \item{p4}{poisson, log()=1 }
#'   \item{p5}{poisson, log()=1 }
#'   \item{p6}{poisson, log()=1 }
#'   \item{p7}{poisson, log()=1 }
#'   \item{p8}{poisson, log()=1 }
#'   \item{p9}{poisson, log()=1 }
#'   \item{p10}{poisson, log()=1 }
#'}
#' @examples
#' \dontrun{
#' ## The dataset was (essentially) generated using the following code:
#' datasize <- 300
#' tmp <- c(rep("y", as.integer(datasize/2)), rep("n", as.integer(datasize/2)))
#' set.seed(1)
#'
#' ex0.dag.data <- data.frame(b1=sample(tmp, size=datasize, replace=TRUE),
#'                            b2=sample(tmp, size=datasize, replace=TRUE),
#'                            b3=sample(tmp, size=datasize, replace=TRUE),
#'                            b4=sample(tmp, size=datasize, replace=TRUE),
#'                            b5=sample(tmp, size=datasize, replace=TRUE),
#'                            b6=sample(tmp, size=datasize, replace=TRUE),
#'                            b7=sample(tmp, size=datasize, replace=TRUE),
#'                            b8=sample(tmp, size=datasize, replace=TRUE),
#'                            b9=sample(tmp, size=datasize, replace=TRUE),
#'                            b10=sample(tmp, size=datasize, replace=TRUE),
#'                            g1=rnorm(datasize, mean=0,sd=1),
#'                            g2=rnorm(datasize, mean=0,sd=1),
#'                            g3=rnorm(datasize, mean=0,sd=1),
#'                            g4=rnorm(datasize, mean=0,sd=1),
#'                            g5=rnorm(datasize, mean=0,sd=1),
#'                            g6=rnorm(datasize, mean=0,sd=1),
#'                            g7=rnorm(datasize, mean=0,sd=1),
#'                            g8=rnorm(datasize, mean=0,sd=1),
#'                            g9=rnorm(datasize, mean=0,sd=1),
#'                            g10=rnorm(datasize, mean=0,sd=1),
#'                            p1=rpois(datasize, lambda=10),
#'                            p2=rpois(datasize, lambda=10),
#'                            p3=rpois(datasize, lambda=10),
#'                            p4=rpois(datasize, lambda=10),
#'                            p5=rpois(datasize, lambda=10),
#'                            p6=rpois(datasize, lambda=10),
#'                            p7=rpois(datasize, lambda=10),
#'                            p8=rpois(datasize, lambda=10),
#'                            p9=rpois(datasize, lambda=10),
#'                            p10=rpois(datasize, lambda=10))
#' }
#' @keywords datasets internal
"ex0.dag.data"

#' @title Synthetic validation data set for use with abn library examples
#' @description
#' 10000 observations simulated from a DAG with 10 variables from Poisson, Bernoulli and Gaussian distributions.
#' @usage ex1.dag.data
#'
#' @format A data frame, binary variables are factors.
#' The relevant formulas are given below (note these do not give parameter
#' estimates just the form of the relationships, like in glm(),
#' e.g. logit()=1+p1 means a logit link function and comprises of an
#' intercept term and a term involving p1).
#' \describe{
#'   \item{b1}{binary, logit()=1 }
#'   \item{p1}{poisson, log()=1 }
#'   \item{g1}{gaussian, identity()=1 }
#'   \item{b2}{binary, logit()=1}
#'   \item{p2}{poisson, log()=1+b1+p1 }
#'   \item{b3}{binary, logit()=1+b1+g1+b2 }
#'   \item{g2}{gaussian, identify()=1+p1+g1+b2 }
#'   \item{b4}{binary, logit()=1+g1+p2}
#'   \item{b5}{binary, logit()=1+g1+g2 }
#'   \item{g3}{gaussian, identity()=1+g1+b2 }
#' }
#' @examples
#' ## The data is one realisation from the the underlying DAG:
#' ex1.true.dag <- matrix(data=c(
#'   0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,
#'   1,1,0,0,0,0,0,0,0,0,
#'   1,0,1,1,0,0,0,0,0,0,
#'   0,1,1,1,0,0,0,0,0,0,
#'   0,0,1,0,1,0,0,0,0,0,
#'   0,0,1,0,0,0,1,0,0,0,
#'   0,0,1,1,0,0,0,0,0,0), ncol=10, byrow=TRUE)
#'
#' colnames(ex1.true.dag) <- rownames(ex1.true.dag) <-
#'   c("b1","p1","g1","b2","p2","b3","g2","b4","b5","g3")
#' @keywords datasets internal
"ex1.dag.data"

#' @title Synthetic validation data set for use with abn library examples
#' @description
#' 10000 observations simulated from a DAG with 18 variables three sets each from Poisson, Bernoulli and Gaussian distributions.
#' @usage ex2.dag.data
#' @format A data frame, binary variables are factors.
#' The relevant formulas are given below (note these do not give parameter
#' estimates just the form of the relationships, e.g. logit()=1
#' means a logit link function and comprises of only an intercept term).
#' \describe{
#'   \item{b1}{binary,logit()=1+g1+b2+b3+p3+b4+g4+b5}
#'   \item{g1}{gaussian,identity()=1}
#'   \item{p1}{poisson,log()=1+g6}
#'   \item{b2}{binary,logit()=1+p3+b4+p6}
#'   \item{g2}{gaussian,identify()=1+b2}
#'   \item{p2}{poisson,log()=1+b2}
#'   \item{b3}{binary,logit()=1+g1+g2+p2+g3+p3+g4}
#'   \item{g3}{gaussian,identify()=1+g1+p3+b4}
#'   \item{p3}{poisson,log()=1}
#'   \item{b4}{binary,logit()=1+g1+p3+p5}
#'   \item{g4}{gaussian,identify()=1+b4;}
#'   \item{p4}{poisson,log()=1+g1+b2+g2+b5}
#'   \item{b5}{binary,logit()=1+b2+g2+b3+p3+g4}
#'   \item{g5}{gaussian,identify()=1}
#'   \item{p5}{poisson,log()=1+g1+g5+b6+g6}
#'   \item{b6}{binary,logit()=1}
#'   \item{g6}{gaussian,identify()=1}
#'   \item{p6}{poisson,log()=1+g5}
#' }
#' @examples
#' ## The true underlying stochastic model has DAG - this data is a single realisation.
#' ex2.true.dag <- matrix(data = c(
#'   0,1,0,1,0,0,1,0,1,1,1,0,1,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
#'   0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,
#'   0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,1,0,0,1,1,0,1,1,0,1,0,0,0,0,0,0,0,
#'   0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,
#'   0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
#'   0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,
#'   0,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'   0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0
#' ), ncol = 18, byrow = TRUE)
#'
#' colnames(ex2.true.dag) <- rownames(ex2.true.dag) <- c("b1","g1","p1","b2",
#'                                                       "g2","p2","b3","g3",
#'                                                       "p3","b4","g4","p4",
#'                                                       "b5","g5","p5","b6",
#'                                                       "g6","p6")
#' @keywords datasets internal
"ex2.dag.data"

#' @title Validation data set for use with abn library examples
#' @description
#' 1000 observations across with 13 binary variables and one grouping variable. Real (anonymised) data of unknown structure.
#' @usage ex3.dag.data
#' @format A data frame with 14 columns, where
#' \code{b1,b2,\dots,b13} are binary variables encoded as factors and
#' \code{group} is a factor with 100 factors  defining the sampling
#' groups (10 observations each).
#'
#' @keywords datasets internal
"ex3.dag.data"

#' @title Valdiation data set for use with abn library examples
#' @description
#' 2000 observations across with 10 binary variables and one grouping variable. Real (anonymised) data of unknown structure.
#' @usage ex4.dag.data
#' @format A data frame with eleven columns:
#' \code{group} factor with 85 levels defining sampling groups;
#' \code{b1,\dots,b10} binary variables encoded as factors.
#'
#' @keywords datasets internal
"ex4.dag.data"

#' @title Valdiation data set for use with abn library examples
#' @description
#' 434 observations across with 18 variables, 6 binary and 12 continuous, and one grouping variable. Real (anonymised) data of unknown structure.
#' @usage ex5.dag.data
#' @format A data frame with 19 columns: \code{b1,\dots,b6} binary
#' variables, encoded as factors;
#' \code{g1,\dots,g12} continuous variables. Finally, the column
#' \code{group} defines sampling groups (encoded as a factor as well).
#'
#' @keywords datasets internal
"ex5.dag.data"

#' @title Valdiation data set for use with abn library examples
#' @description
#' 800 observations across with 8 variables, 1 count, 2 binary and 4 continuous, and 1 grouping variable. Real (anonymised) data of unknown structure.
#' @usage ex6.dag.data
#' @format A data frame with eight columns. Binary variables are factors
#' \describe{
#'   \item{p1}{count}
#'   \item{g1}{continuous}
#'   \item{g2}{continuous}
#'   \item{b1}{binary}
#'   \item{b2}{binary}
#'   \item{g3}{continuous}
#'   \item{g4}{continuous}
#'   \item{group}{factor,defines sampling groups}
#' }
#'
#' @keywords datasets internal
"ex6.dag.data"

#' @title Valdiation data set for use with abn library examples
#' @description
#' 10648 observations across with 3 variables, 2 binary and 1 grouping variable.
#' Real (anonymised) data of unknown structure.
#' @usage ex7.dag.data
#' @format A data frame, binary variables are factors
#' \describe{
#'   \item{b1}{binary}
#'   \item{b2}{binary}
#'   \item{group}{factor, defines sampling groups}
#' }
#'
#' @keywords datasets internal
"ex7.dag.data"

#' @title Dataset related to Feline calicivirus infection among cats in Switzerland.
#' @description
#' The dataset is about the Feline calicivirus (FCV) infection among cats in Switzerland.
#' FCV is a virus that occurs worldwide in domestic cats but also in exotic felids.
#' FCV is a highly contagious virus that is the major cause of upper respiratory
#' disease or cat flue that affects felids. This is a complex disease caused by
#' different viral and bacterial pathogens, i.e., FCV, FHV-1, \emph{Mycoplasma felis},
#' \emph{Chlamydia felis} and \emph{Bordetella bronchiseptica}.
#' It can be aggravated by retrovirus infections such as FeLV and FIV.
#' This composite dynamic makes it very interesting for a BN modeling approach.
#' The data were collected between September 2012 and April 2013.
#' @usage FCV
#' @format An adapted data frame of the original dataset, which consists of 300 observations of 15 variables.
#' \describe{
#'  \item{FCV}{Feline Calici Virus status (0/1).}
#'  \item{FHV_1}{Feline Herpes Virus 1 status (0/1). }
#'  \item{C_felis}{C-felis and Chlamydia felis status (0/1).}
#'  \item{M_felis}{Mycoplasma felis status (0/1).}
#'  \item{B_bronchiseptica}{B-bronchiseptica & Bordetella bronchispetica status (0/1).}
#'  \item{FeLV}{feline leukosis virus status (0/1).}
#'  \item{FIV}{feline immunodeficiency virus status (0/1).}
#'  \item{Gingivostomatitis}{gingivostomatitis complex status (0/1).}
#'  \item{URTD}{URTD complex (upper respiratory complex) (0/1).}
#'  \item{Vaccinated}{vaccination status (0/1).}
#'  \item{Pedigree}{pedigree (0/1).}
#'  \item{Outdoor}{outdoor access (0/1).}
#'  \item{Sex}{sex and castrated status (M, MN, F, FS).}
#'  \item{GroupSize}{number of cats in the group (counts).}
#'  \item{Age}{age in year (continuous)\.}
#' }
#' @references Berger, A., Willi, B., Meli, M. L., Boretti, F. S., Hartnack, S., Dreyfus, A., ... and Hofmann-Lehmann, R. (2015). Feline calicivirus and other respiratory pathogens in cats with Feline calicivirus-related symptoms and in clinically healthy cats in Switzerland. BMC Veterinary Research, 11(1), 282.
#'
#' @keywords datasets internal
"FCV"

#' @title Dataset related to diseases present in `finishing pigs', animals about to enter the human food chain at an abattoir.
#' @description
#' The data we consider here comprise of a randomly chosen batch of 50 pigs from
#' each of 500 randomly chosen pig producers in the UK.
#' The dataset consists of 25000 observations, 10 binary variables, and a grouping variable.
#' These are `finishing pigs', animals about to enter the human food chain at an abattoir.
#' Further description of the data set is present on the vignette.
#' @usage pigs.vienna
#' @format A data frame with a mixture of 10 discrete variables, each of which is set as a factor, and a grouping variable.
#' \describe{
#'  \item{PC}{Binary.}
#'  \item{PT}{Binary. }
#'  \item{MS}{Binary.}
#'  \item{HS}{Binary.}
#'  \item{TAIL}{Binary.}
#'  \item{Abscess}{Binary.}
#'  \item{Pyaemia}{Binary.}
#'  \item{EPcat}{Binary.}
#'  \item{PDcat}{Binary.}
#'  \item{plbinary}{Binary.}
#'  \item{batch}{Group variable, corresponding to the 500 different pig producers}
#' }
#' @details
#' This dataset was used in an older version of the vignette.
#' See also the files provided in the package directories
#' \code{inst/bootstrapping_example} and \code{inst/old_vignette} give a
#' detailed analysis of the dataset and provide more details for a
#' bootstrapping example thereof.
#' @references Hartnack, S., et al. (2016) "Attitudes of Austrian veterinarians towards euthanasia in small animal practice: impacts of age and gender on views on euthanasia." BMC Veterinary Research 12.1: 26.
#'
#' @keywords datasets internal
"pigs.vienna"

#' @title simulated dataset from a DAG comprising of 33 variables
#' @description
#' 250 observations simulated from a DAG with 17 binary variables and 16 continuous.
#' A DAG of this data features in the vignette.
#' Note that the conditional dependence relations given are those in the population
#' and may differ in the realization of 250 observations.
#' @usage var33
#' @format A data frame with a mixture of discrete variables each of which is
#' set as a factor and continuous variables.
#' Joint distribution structure used to generate the data.
#' \describe{
#'  \item{v1}{Binary, independent. }
#'  \item{v2}{Gaussian, conditionally dependent upon v1. }
#'  \item{v3}{Binary, independent. }
#'  \item{v4}{Binary, conditionally dependent upon v3. }
#'  \item{v5}{Gaussian, conditionally dependent upon v6. }
#'  \item{v6}{Binary, conditionally dependent upon v4 and v7. }
#'  \item{v7}{Gaussian,  conditionally dependent upon v8. }
#'  \item{v8}{Gaussian,  conditionally dependent upon v9. }
#'  \item{v9}{Binary, conditionally dependent upon v10. }
#'  \item{v10}{Binary, independent. }
#'  \item{v11}{Binary,  conditionally dependent upon v10, v12 and v19. }
#'  \item{v12}{Binary, independent.}
#'  \item{v13}{Gaussian, independent.}
#'  \item{v14}{Gaussian,  conditionally dependent upon v13. }
#'  \item{v15}{Binary, conditionally dependent upon v14 and v21. }
#'  \item{v16}{Gaussian, independent. }
#'  \item{v17}{Gaussian, conditionally dependent upon v16 and v20. }
#'  \item{v18}{Binary,  conditionally dependent upon v20. }
#'  \item{v19}{Binary,  conditionally dependent upon v20. }
#'  \item{v20}{Binary, independent. }
#'  \item{v21}{Binary, conditionally dependent upon v20. }
#'  \item{v22}{Gaussian, conditionally dependent upon v21. }
#'  \item{v23}{Gaussian, conditionally dependent upon v21. }
#'  \item{v24}{Gaussian, conditionally dependent upon v23. }
#'  \item{v25}{Gaussian, conditionally dependent upon v23 and v26. }
#'  \item{v26}{Binary, conditionally dependent upon v20. }
#'  \item{v27}{Binary, independent. }
#'  \item{v28}{Binary, conditionally dependent upon v27, v29 and v31. }
#'  \item{v29}{Gaussian, independent. }
#'  \item{v30}{Gaussian,  conditionally dependent upon v29. }
#'  \item{v31}{Gaussian, independent. }
#'  \item{v32}{Binary,  conditionally dependent upon v21, v29 and v31. }
#'  \item{v33}{Gaussian,  conditionally dependent upon v31. }
#' }
#' @examples
#' ## Constructing the DAG of the dataset:
#' dag33 <- matrix(0, 33, 33)
#' dag33[2,1] <- 1
#' dag33[4,3] <- 1
#' dag33[6,4] <- 1; dag33[6,7] <- 1
#' dag33[5,6] <- 1
#' dag33[7,8] <- 1
#' dag33[8,9] <- 1
#' dag33[9,10] <- 1
#' dag33[11,10] <- 1; dag33[11,12] <- 1; dag33[11,19] <- 1;
#' dag33[14,13] <- 1
#' dag33[17,16] <- 1; dag33[17,20] <- 1
#' dag33[15,14] <- 1; dag33[15,21] <- 1
#' dag33[18,20] <- 1
#' dag33[19,20] <- 1
#' dag33[21,20] <- 1
#' dag33[22,21] <- 1
#' dag33[23,21] <- 1
#' dag33[24,23] <- 1
#' dag33[25,23] <- 1; dag33[25,26] <- 1
#' dag33[26,20] <- 1
#' dag33[33,31] <- 1
#' dag33[33,31] <- 1
#' dag33[32,21] <- 1; dag33[32,31] <- 1; dag33[32,29] <- 1
#' dag33[30,29] <- 1
#' dag33[28,27] <- 1; dag33[28,29] <- 1; dag33[28,31] <- 1
#'
#' @keywords datasets internal
"var33"

#' @title Toy Data Set for Examples in README
#' @description
#' 1000 observations with 5 variables: 2 continuous, 2 binary and 1 categorical.
#' @usage g2b2c_data
#' @format A data frame with five columns. Binary and categorical variables are factors.
#' \describe{
#'   \item{G1}{gaussian}
#'   \item{B1}{binomial}
#'   \item{B2}{binomial}
#'   \item{C}{categorical}
#'   \item{G2}{gaussian}
#' }
#'
#' @keywords datasets internal
"g2b2c_data"

#' @title Toy Data Set for Examples in README
#' @description
#' 10000 observations with 6 variables: 2 continuous, 1 binary, 1 count, 1 categorical and 1 grouping factor.
#' @usage g2pbcgrp
#' @format A data frame with six columns. Binary and categorical variables are factors.
#' \describe{
#'   \item{G1}{gaussian}
#'   \item{P}{poisson}
#'   \item{B}{binomial}
#'   \item{C}{categorical}
#'   \item{G2}{gaussian}
#'   \item{group}{categorical}
#' }
#'
#' @keywords datasets internal
"g2pbcgrp"
