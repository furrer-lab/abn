## code to prepare `g2b2c_data` dataset goes here
# Set the seed for reproducibility
set.seed(123)

# Number of observations
n <- 1000

# Generate Gaussian distributed variable
gaussian_var <- rnorm(n, mean = 0, sd = 1)

# Generate Binomial distributed variable influenced by the Gaussian variable
prob <- pnorm(gaussian_var)
binomial_var <- rbinom(n, size = 1, prob = prob)

# Generate second Binomial distributed variable influenced by the Gaussian variable
prob <- pnorm(gaussian_var)
binomial_var2 <- rbinom(n, size = 1, prob = prob)

# Generate Categorical distributed variable influenced by the Poisson and Binomial variables
prob <- cbind(binomial_var2, binomial_var, abs(gaussian_var)) / rowSums(cbind(binomial_var2, binomial_var, abs(gaussian_var)))
categorical_var <- factor(apply(prob, 1, function(x) sample(letters[c(1:3)], size = 1, prob = x)))

# Generate another Gaussian variable that is directly dependent on the Categorical variable
gaussian_var2 <- rnorm(n, mean = as.integer(categorical_var), sd = 1)

# Combine all variables into a data frame
g2b2c_data <- data.frame("G1" = gaussian_var,
                         "B1" = factor(binomial_var),
                         "B2" = factor(binomial_var2),
                         "C" = factor(categorical_var),
                         "G2" = gaussian_var2)

# Save the data frame as an RDS file
usethis::use_data(g2b2c_data, overwrite = TRUE)
