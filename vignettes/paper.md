---
title: "Additive Bayesian Networks"
tags:
- data science
- R
- mixed-effects models
- Bayesian networks
- graphical models
authors:
- name: Matteo Delucchi
  orcid: 0000-0002-9327-1496
  affiliation: "1, 2"
- name: Jonas I. Liechti
  orcid: 0000-0003-3447-3060
  affiliation: "3"
- name: Georg R. Spinner
  orcid: 0000-0001-9640-8155
  affiliation: "2"
- name: Reinhard Furrer
  orcid: 0000-0002-6319-2332
  corresponding: true
  affiliation: "1"
affiliations:
 - name: Department of Mathematical Modeling and Machine Learning, University of Zurich, Zürich, Switzerland
   index: 1
 - name: Centre for Computational Health, Institute of Computational Life Sciences, Zurich University of Applied Sciences (ZHAW), Wädenswil, Switzerland
   index: 2
 - name: www.T4D.ch, T4D GmbH, Zurich, Switzerland
   index: 3
date: 20. Mai 2024
bibliography: paper.bib
---

# Summary
The R package `abn` is a comprehensive tool for Bayesian Network (BN) analysis, a form of probabilistic graphical model. 
BNs are a type of statistical model that leverages the principles of Bayesian statistics and graph theory to provide a framework for representing complex multivariate data. 
They can derive a directed acyclic graph from empirical data to describe the dependency structure between random variables. 

Additive Bayesian Network (ABN) models extend the concept of generalized linear models, typically used for predicting a single outcome, to scenarios with multiple dependent variables (e.g., @kratzer_additive_2023).
This makes them a powerful tool for understanding complex, multivariate datasets.
This package provides routines for structure learning and parameter estimation of ABN models.

# Statment of need
The increasing complexity of data in various fields, ranging from healthcare research to environmental science and ecology, has resulted in a need for a tool like `abn`.
Researchers often face multivariate, tabular data where the relationships between variables are not straightforward. 
BN analysis becomes essential when traditional statistical methods fail to analyze multivariate data with intricate relationships, as it models these relationships graphically for more straightforward data interpretation.

Commonly used implementations of BN models, such as `bnlearn` [@bnlearn2010], `bnstruct` [@franzin_bnstruct_2017], `deal` [@boettcher_deal_2003], `gRain` [@hojsgaard_graphical_2012], `pcalg` [@kalisch_causal_2012] and `pchc` [@tsagris_new_2021], limit variable types, often allowing discrete variables to have only discrete parent variables, where a parent starts a directed edge in the graph.
This limitation can pose challenges when dealing with continuous or mixed-type data (i.e., data that includes both continuous and discrete variables) or when attempting to model complex relationships that do not fit these restricted categories.
For a comprehensive overview of structure learning algorithms, including those applicable to mixed-type data, we refer the reader to the works of @kitson_survey_2023 and @zanga_survey_2022.
In the context of patient data, the study from @delucchi_bayesian_2022 has discussed further details and strategies for handling this scenario, particularly in relation to the `abn` package and the widely used `bnlearn` package [@bnlearn2010].

The `abn` package overcomes this limitation through its additive model formulation, which generalizes the usual (Bayesian) multivariable regression to accommodate multiple dependent variables.
Additionally, the `abn` package offers a comprehensive suite of features for model selection, structure learning, and parameter estimation.
It includes exact and greedy search algorithms for structure learning and allows for integrating prior expert knowledge into the model selection process by specifying structural constraints.
For model selection, a Bayesian and an information-theoretic model scoring approach are available, allowing users to choose between a Bayesian and frequentist paradigm.
To our knowledge, this feature is not available in other software.
Furthermore, it supports mixed-effect models to control one-layer clustering, making it suitable, e.g., for handling data from different sources.

Previous versions of the `abn` package have been successfully used in various fields, including epidemiology [@pittavino_comparison_2017, @kratzer_information-theoretic_2018] and health [@hartnack_additive_2019, @kratzer_bayesian_2020, @delucchi_bayesian_2022].
Despite its promise, the `abn` package encountered historical obstacles.
Sporadic maintenance and an incomplete codebase hindered its full potential. 
Recognizing the need for enhancement, we undertook a substantial upgrade and meticulously addressed legacy issues, revamped the codebase, and introduced significant improvements. 
The latest version 3 of `abn` is now a robust and reliable tool for BN analysis.
Applying the latest standards for open-source software, we guarantee active maintenance of `abn`. 
Future updates are planned to enhance its functionality and user experience further. 
We highly value feedback from the user community, which will guide our ongoing developments.

In summary, `abn` sets itself apart by emphasizing ABNs and its exhaustive features for model selection and structure learning. 
Its unique contribution is the implementation of mixed-effect BN models, thereby extending its applicability to a broader range of complex, multivariate datasets of mixed, continuous, and discrete data.

# Implementation
As outlined in @kratzer_additive_2023, the package's comprehensive framework integrates the mixed-effects model for clustered data, considering data heterogeneity and grouping effects.
However, this was confined to a Bayesian context and was only a preliminary implementation.
With the release of `abn` major version 3, this was completed with an implementation under the information-theoretic (`method = "mle"`) setting.

Analyzing hierarchical or grouped data, i.e., observations nested within higher-level units, requires statistical models with group-varying parameters (e.g., mixed-effect models). 
The `abn` package facilitates single-layer clustering, where observations are grouped. 
These clusters are assumed to be independent, but intra-cluster observations may exhibit correlation (e.g., students within schools, patient-specific measurements over time, etc.). 
The ABN model is fitted independently as a varying intercept model, where the intercept can vary while the slope is assumed constant across all group levels.

Under the frequentist paradigm (`method = "mle"`), `abn` employs the `lme4` package [@lme42015] to fit generalized linear mixed models for each of the Binomial, Poisson, and Gaussian distributed variables. 
For multinomial distributed variables, `abn` fits a multinomial baseline category logit model with random effects using the `mclogit` package [@mclogit2022]. 
Currently, only single-layer clustering is supported (e.g., for `method = "mle"`, this corresponds to a random intercept model).

With a Bayesian approach (`method = "bayes"`), `abn` utilizes its own implementation of the Laplace approximation as well as the `INLA` package [@inla2013] to fit a single-level hierarchical model for Binomial, Poisson, and Gaussian distributed variables. 

Furthermore, the code base has been enhanced to be more efficient, reliable, and user-friendly through code optimization, regular reviews, and continuous integration practices. 
We have adhered to the latest open-source software standards, including active maintenance of `abn`. 
Future updates to augment its functionality are planned via a flexible roadmap.
User feedback is valued through open communication channels, which will steer our ongoing developments. 
Consequently, the latest version of `abn` is now a robust and reliable tool for BN analysis.

# Validation and Testing
A comprehensive set of documented case studies has been published to validate the `abn` package (see the `abn` [website](https://r-bayesian-networks.org/)).
The numerical accuracy and quality assurance exercises were demonstrated in @kratzer_additive_2023.
A rigorous testing framework is implemented using the `testthat` package [@testthat2011], which is executed as part of an extensive continuous integration pipeline designed explicitly for non-standard R packages that rely on `Rcpp` [@rcpp2023] and `JAGS` [@plummer_jags_2003].
Additional documentation and resources are available on the `abn` [website](https://r-bayesian-networks.org/) for further reference and guidance.

# Availability

The latest version of the `abn` package is available on [GitHub](https://github.com/furrer-lab/abn) and can be installed using the `devtools` package:

```r
devtools::install_github("furrer-lab/abn")
```

# Acknowledgments

The development of the `abn` package would not have been possible without the significant contributions of the former developers whose efforts have been instrumental in shaping this project. 
We acknowledge the contributions of Fraser Iain Lewis, Marta Pittavino, Gilles Kratzer, and Kalina Cherneva, in particular.
We extend our gratitude to the faculty staff at the [Department of Mathematical Modeling and Machine Learning](https://dm3l.uzh.ch/home), University of Zurich (UZH), and the [Department of Mathematics](https://www.math.uzh.ch/home), UZH, who maintain the research and teaching infrastructure.
Our appreciation also goes to the UZH and the ZHAW for their financial support. 
We want to highlight the funding from both the Zurich University of Applied Sciences (ZHAW) and the Digitalization Initiative of the Zurich Higher Education Institutions (DIZH), which were instrumental in realizing this project, particularly within the context of the "Modeling of multicentric and dynamic stroke health data" and "Stroke DynamiX" projects, respectively.
This work was conducted as part of M.D.'s PhD project, co-supervised by Prof. Dr. Sven Hirsch (ZHAW) and Prof. Dr. Reinhard Furrer (UZH).

# References

