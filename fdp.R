library(dplyr)
library(glmnet)
library(magrittr)

generateY <- function(betas, X) {
    n <- nrow(X)
    y <- X %*% betas + rnorm(n=n, mean=0, sd=1)
    return(y)
}

getDesignMatrix <- function(betas, n) {
    mat <- {}
    p <- length(betas)
    for(i in 1:n) {
        mat <- rbind(mat, generateX(betas))
    }
    return(mat)
}

generateX <- function(betas) {
    p <- length(betas)
    x <- rnorm(n=p, mean=0, sd=1)
    return(x)
}

getFDP <- function(beta_true, beta_hat) {
    numerator <- ifelse(beta_hat != 0 & beta_true == 0, 1, 0) %>%
        sum
    denominator <- beta_hat[beta_hat != 0] %>% 
        length
    return(numerator/denominator)
}

getTPP <- function(beta_true, beta_hat) {
    numerator <- ifelse(beta_hat != 0 & beta_true != 0, 1, 0) %>%
        sum
    denominator <- beta_true[beta_true != 0] %>% 
        length
    return(numerator/denominator)
}

n <- 1010
betas <- c(rep(4, 200), rep(0, 800))

X <- getDesignMatrix(betas, n)
y <- generateY(betas, X)

lambdas <- 1

nlambda <- 50000
max_lambda <- 5*sum(betas)
lambdas <- max_lambda*exp(-seq(0.001, 10, length.out=nlambda))
lasso_mod <- glmnet::glmnet(X, y, alpha=1, lambda=lambdas)
beta_lasso <- lasso_mod$beta

ols_mod <- lm(y ~ X)
beta_ols <- ols_mod$coef

ols_order <- beta_ols %>%
    order

fdp <- {}
tpp <- {}
for(i in 1:length(lambdas)) {
    beta_hat <- beta_lasso[, i]
    fdp[i] <- getFDP(betas, beta_hat)
    tpp[i] <- getTPP(betas, beta_hat)
}

# Reproduce Figure 1 of Su et al.
pdf("Figures/Figure1.pdf", width=10, height=6)
plot(fdp ~ tpp, pch=19, lwd=1)
dev.off()
