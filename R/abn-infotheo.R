#' Discretization of a Possibly Continuous Data Frame of Random Variables based on their distribution
#'
#' @usage discretization(data.df = NULL,
#'                       data.dists = NULL,
#'                       discretization.method = "sturges",
#'                       nb.states = FALSE)
#'
#'
#' @param data.df a data frame containing the data to discretize, binary and multinomial variables must be declared as factors, others as a numeric vector. The data frame must be named.
#' @param data.dists a named list giving the distribution for each node in the network.
#' @param discretization.method a character vector giving the discretization method to use; see details. If a number is provided, the variable will be discretized by equal binning.
#' @param nb.states logical variable to select the output. If set to \code{TRUE} a list with the discretized data frame and the number of state of each variable is returned. If set to FALSE only the discretized data frame is returned.
#'
#' @description
#' This function discretizes a data frame of possibly continuous random variables
#' through rules for discretization. The discretization algorithms are unsupervised
#' and univariate. See details for the complete list of discretization rules (the number of state of each
#' random variable could also be provided).
#'
#' @details
#' \code{fd} Freedman Diaconis rule. \code{IQR()} stands for interquartile range.
#' The number of bins is given by \deqn{\frac{range(x) * n^{1/3}}{2 * IQR(x)}}{range(x) * n^{1/3} / 2 * IQR(x)}
#' The Freedman Diaconis rule is known to be less sensitive than the Scott's rule to outlier.
#'
#' \code{doane} Doane's rule.
#' The number of bins is given by \deqn{1 + \log_{2}{n} + \log_{2}{1+\frac{|g|}{\sigma_{g}}}}{1 + \log_{2}{n} + \log_{2}{1+\frac{|g|}{\sigma_{g}}}}
#' This is a modification of Sturges' formula, which attempts to improve its performance with non-normal data.
#'
#' \code{sqrt} The number of bins is given by: \deqn{\sqrt(n)}{\sqrt(n)}
#'
#' \code{cencov} Cencov's rule.
#' The number of bins is given by: \deqn{n^{1/3}}{n^{1/3}}
#'
#' \code{rice} Rice' rule.
#' The number of bins is given by: \deqn{2 n^{1/3}}{2 n^{1/3}}
#'
#' \code{terrell-scott} Terrell-Scott's rule.
#' The number of bins is given by: \deqn{(2 n)^{1/3}}{(2 n)^{1/3}}
#' It is known that Cencov, Rice, and Terrell-Scott rules over-estimates k, compared to other rules due to its simplicity.
#'
#' \code{sturges} Sturges's rule.
#' The number of bins is given by: \deqn{1 + \log_2(n)}{1 + \log_2(n)}
#'
#' \code{scott} Scott's rule.
#' The number of bins is given by: \deqn{range(x) / \sigma(x) n^{-1/3}}{range(x) / \sigma(x) n^{-1/3}}
#'
#' @return The discretized data frame or a list containing the table of
#' counts for each bin the discretized data frame.
#'
#' @references
#' Garcia, S., et al.  (2013). A survey of discretization techniques: Taxonomy and empirical analysis in supervised learning. \emph{IEEE Transactions on Knowledge and Data Engineering}, 25.4, 734-750.
#'
#' Cebeci, Z. and Yildiz, F. (2017). Unsupervised Discretization of Continuous Variables in a Chicken Egg Quality Traits Dataset. \emph{Turkish Journal of Agriculture-Food Science and Technology}, 5.4, 315-320.
#'
#' @export
#'
#' @examples
#' ## Generate random variable
#' rv <- rnorm(n = 100, mean = 5, sd = 2)
#' dist <- list("gaussian")
#' names(dist) <- c("rv")
#'
#' ## Compute the entropy through discretization
#' entropyData(freqs.table = discretization(data.df = rv, data.dists = dist,
#' discretization.method = "sturges", nb.states = FALSE))
#' @keywords utilities
#' @returns table of counts for each bin of the discretized data frame.
#' @importFrom stats sd IQR
discretization <- function(data.df = NULL, data.dists = NULL, discretization.method = "sturges", nb.states = FALSE) {

    # tests:

    if (is.character(discretization.method))
        discretization.method <- tolower(discretization.method)

    if (is.character(discretization.method))
        discretization.method <- c("fd", "doane", "cencov", "sturges", "rice", "scott", "sqrt", "terrell-scott")[pmatch(discretization.method, c("fd", "doane", "cencov", "sturges", "rice", "scott",
            "sqrt", "terrell-scott"))]

    if (!(discretization.method %in% c("fd", "doane", "cencov", "sturges", "rice", "scott", "sqrt", "terrell-scott") | is.numeric(discretization.method))) {
        stop("Wrong definition of the discretization.method's parameter")
    }

    df.tab <- list()
    nb.states.l <- list()
    data.df <- as.data.frame(data.df)

    ## =============================================================== Discrete cutoff provided ===============================================================

    if (is.numeric(discretization.method)) {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = discretization.method + 1)
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Freedman-Diaconis rule ===============================================================

    if (discretization.method == "fd") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = (abs(range(data.df[, i])[2] - range(data.df[, i])[1])/diff(range(data.df[, i]))/(2 * IQR(data.df[,
                  i])/length(data.df[, i])^(1/3))) + 1)
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Doane s formula ===============================================================

    if (discretization.method == "doane") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = 1 + log2(length(data.df[, i])) + log2(1 + (abs(skewness(data.df[, i]))/sqrt((6 * (length(data.df[,
                  i]) - 2))/((length(data.df[, i]) + 1) * (length(data.df[, i]) + 3))))))
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Square root choice ===============================================================

    if (discretization.method == "sqrt") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = ceiling(abs(sqrt(length(data.df[, i])) + 1)))
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Sturges ===============================================================

    if (discretization.method == "sturges") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = log(x = length(data.df[, i]), base = 2) + 2)
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Rice Rule ===============================================================

    if (discretization.method == "rice") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = 2 * length(data.df[, i])^(1/3) + 1)
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }

    ## =============================================================== Scott rule ===============================================================

    if (discretization.method == "scott") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = (abs(range(data.df[, i])[2] - range(data.df[, i])[1])/diff(range(data.df[, i]))/(3.5 *
                  sd(data.df[, i])/length(data.df[, i])^(1/3))) + 1)
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }
    ## =============================================================== Terrell-Scott ===============================================================

    if (discretization.method == "terrell-scott") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = ceiling(abs(2 * nobs)^(1/3)))
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }



    ## =============================================================== cencov ===============================================================

    if (discretization.method == "cencov") {

        for (i in 1:length(data.dists)) {

            if (data.dists[i] == "binomial" | data.dists[i] == "multinomial") {
                cut.index1 <- seq(from = range(as.numeric(data.df[, i]))[1], to = range(as.numeric(data.df[, i]))[2], length.out = length(levels(data.df[, i])) + 1)
            } else {
                cut.index1 <- seq(from = range(data.df[, i])[1], to = range(data.df[, i])[2], length.out = ceiling(abs(nobs^(1/3))))
            }

            df.tab[[i]] <- as.factor(cut(as.numeric(data.df[, i]), breaks = cut.index1, include.lowest = TRUE))
            nb.states.l[[i]] <- length(cut.index1)
        }
    }





    if (nb.states) {
        return(list(df = table(df.tab), nb.states = nb.states.l))
    } else {
        return(table(df.tab))
    }
}

#' Empirical Estimation of the Entropy from a Table of Counts
#'
#' This function empirically estimates the Mutual Information from a table of counts using the observed frequencies.
#'
#' @usage miData(freqs.table, method = c("mi.raw", "mi.raw.pc"))
#'
#' @param freqs.table a table of counts.
#' @param method a character determining if the Mutual Information should be normalized.
#'
#' @details
#' The mutual information estimation is computed from the observed frequencies through a plugin estimator based on entropy.
#'
#' The plugin estimator is \deqn{I(X, Y) = H (X) + H(Y) - H(X, Y)}, where \deqn{H()} is the entropy computed with \code{\link{entropyData}}.
#'
#' @return Mutual information estimate.
#' @seealso \code{\link{discretization}}
#' @references Cover, Thomas M, and Joy A Thomas. (2012). "Elements of Information Theory". John Wiley & Sons.
#' @export
#' @returns integer
#' @examples
#' ## Generate random variable
#' Y <- rnorm(n = 100, mean = 0, sd = 2)
#' X <- rnorm(n = 100, mean = 5, sd = 2)
#'
#' dist <- list(Y="gaussian", X="gaussian")
#'
#' miData(discretization(data.df = cbind(X,Y), data.dists = dist,
#'                       discretization.method = "fd", nb.states = FALSE),
#'                       method = "mi.raw")
#' @keywords utilities
miData <- function(freqs.table, method = c("mi.raw", "mi.raw.pc")) {

    # normalization
    freqs.joint <- prop.table(freqs.table)

    # xyz ...
    for (i in 1:length(dim(freqs.joint))) {
        if (i == 1) {
            freqs.xy <- margin.table(x = freqs.table, margin = 1)
        } else {
            freqs.xy <- outer(X = freqs.xy, Y = margin.table(x = freqs.table, margin = i), FUN = "*")
        }
    }

    freqs.xy <- freqs.xy/sum(freqs.xy)

    # computing log part
    log.part <- ifelse(freqs.joint > 0, log(freqs.joint/freqs.xy), 0)

    # computing MI
    mi <- sum(freqs.joint * log.part)

    if (method == "mi.raw") {
        return(mi)
    }

    if (method == "mi.raw.pc") {
        tmp <- ifelse((1/margin.table(x = freqs.table, margin = 1)) < Inf, log(1/margin.table(x = freqs.table, margin = 1), base = 2), 0)
        return(-(mi * 100)/sum(margin.table(x = freqs.table, margin = 1) * tmp))
        #
    }

    if (FALSE) {

        # normalization
        freqs.joint <- as.matrix(freqs.table/sum(freqs.table))

        # xy
        freqs.x <- rowSums(freqs.table)
        freqs.y <- colSums(freqs.table)
        freqs.xy <- freqs.x %o% freqs.y
        freqs.xy <- freqs.xy/sum(freqs.xy)

        # computing log part
        log.part <- ifelse(freqs.joint > 0, log(freqs.joint/freqs.xy), 0)

        # computing MI
        mi <- sum(freqs.joint * log.part)

        if (method == "mi.raw") {
            return(mi)
        }

        if (method == "mi.raw.pc") {
            tmp <- ifelse((1/freqs.x) < Inf, log(1/freqs.x, base = 2), 0)
            return(-(mi * 100)/sum(freqs.x * tmp))
            #
        }
    }

}

#' Computes an Empirical Estimation of the Entropy from a Table of Counts
#'
#' This function empirically estimates the Shannon entropy from a table of counts using the observed frequencies.
#'
#' @param freqs.table a table of counts.
#'
#' @details
#' The general concept of entropy is defined for probability distributions.
#' The \code{entropyData()} function estimates empirical entropy from data.
#' The probability is estimated from data using frequency tables.
#' Then the estimates are plug-in in the definition of the entropy to return
#' the so-called empirical entropy. A common known problem of empirical entropy
#' is that the estimations are biased due to the sampling noise.
#' It is also known that the bias will decrease as the sample size increases.
#'
#' @return Shannon's entropy estimate on natural logarithm scale.
#' @seealso \code{\link{discretization}}
#' @references Cover, Thomas M, and Joy A Thomas. (2012). "Elements of Information Theory". John Wiley & Sons.
#' @export
#' @returns integer
#' @examples
#' ## Generate random variable
#' rv <- rnorm(n = 100, mean = 5, sd = 2)
#' dist <- list("gaussian")
#' names(dist) <- c("rv")
#'
#' ## Compute the entropy through discretization
#' entropyData(freqs.table = discretization(data.df = rv, data.dists = dist,
#' discretization.method = "sturges", nb.states = FALSE))
#' @keywords utilities
entropyData <- function(freqs.table) {

    # normalization
    freqs.joint <- prop.table(freqs.table)

    # computing log part
    log.part <- ifelse(freqs.joint > 0, log(freqs.joint, base = 2), 0)

    # computing entropy
    entropy.out <- sum(freqs.joint * log.part)

    return(-entropy.out)

}
