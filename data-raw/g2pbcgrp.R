## code to prepare `g2pbcgrp` dataset goes here
# Generate a grouping factor to simulate a hierarchical structure in the data
n <- 10000
n_groups <- 4
grouping_factor <- sample(1:n_groups, n, replace = TRUE)

# Generate Gaussian distributed variable
group_means <- rnorm(n_groups, mean = 0, sd = 1) # one for each group
gaussian_var <- rnorm(n, mean = group_means[grouping_factor], sd = 1)

# Generate Poisson distributed variable influenced by the Gaussian variable
lambda <- exp(0.5 * gaussian_var)
poisson_var <- rpois(n, lambda)

# Generate Binomial distributed variable influenced by the Gaussian variable
prob <- pnorm(gaussian_var)
binomial_var <- rbinom(n, size = 1, prob = prob)

# Generate Categorical distributed variable influenced by the Poisson and Binomial variables
prob <- cbind(poisson_var, binomial_var, abs(gaussian_var)) / rowSums(cbind(poisson_var, binomial_var, abs(gaussian_var)))
categorical_var <- apply(prob, 1, function(x) sample(1:3, size = 1, prob = x))

# Generate another Gaussian variable that is directly dependent on the Categorical variable
gaussian_var2 <- rnorm(n, mean = categorical_var, sd = 1)


# Combine all variables into a data frame
g2pbcgrp <- data.frame("G1" = gaussian_var,
                   "P" = poisson_var,
                   "B" = factor(binomial_var),
                   "C" = factor(categorical_var),
                   "G2" = gaussian_var2,
                   "group" = factor(grouping_factor))

# Print the first few rows of the data frame
head(g2pbcgrp)

# Save the data frame to a file
usethis::use_data(g2pbcgrp, overwrite = TRUE)
