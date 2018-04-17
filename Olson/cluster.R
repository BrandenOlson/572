library(data.table)
library(dplyr)
library(magrittr)
library(mclust) # For fitting BIC 
library(mvtnorm)

readData <- function(filename, change_names=TRUE) {
    dat <- fread(filename)
    if(change_names) {
        dat_names <- c("x", "y")
        if(ncol(dat) == 3) {
            dat_names <- c(dat_names, "z")
        } 
        names(dat) <- dat_names
    }
    return(dat)
}

getSlot <- function(mixmod_object, slot_name) {
    slot_value <- mixmod_object %>%
        slot("bestResult") %>%
        slot("parameters") %>%
        slot(slot_name)
}

getModelParameters <- function(mod) {
    params <- mod %$% 
        parameters
    means <- params$mean
    variances <- params$variance$sigma
    props <- params$pro
    return(list(Mean=means,
                Variance=variances,
                Prop=props))
}

getModelLikelihood <- function(dat, params) {
    ell <- 0
    prop <- params$Prop
    means <- params$Mean
    variances <- params$Variance
    K <- params$Prop %>% length
    n <- nrow(dat)
    ell <- 0
    for(i in 1:n) {
        log_sum <- 0
        for(k in 1:K) {
            log_sum <- log_sum + prop[k]*dmvnorm(dat[i, ], 
                                                 mean=means[, k],
                                                 sigma=variances[, , k])
        }
        ell <- ell + log(log_sum)
    }
    return(ell)
}

getConditionalProbability <- function(x_i, params, K) {
    prop <- params$Prop
    means <- params$Mean
    variances <- params$Variance
    probs <- {}
    for(k in 1:K) {
        probs[k] <- prop[k]*dmvnorm(x_i,
                                     mean=means[, k],
                                     sigma=variances[, , k]
                                    ) 
    }
    t_ik <- probs/sum(probs)
    z_i <- which.max(t_ik)
    return(z_i)
}

getTs <- function(dat, params, K) {
    ts <- {}
    for(i in 1:nrow(dat)) {
        ts[i] <- getConditionalProbability(dat[i, ],
                                         params,
                                         K)
    }
    return(ts)
}

BIC <- function(log_lik, nu, n) {
    bic <- log_lik - nu*log(n)/2
    return(bic)
}

d1 <- readData("Data/4.1.csv")
d2 <- readData("Data/4.2.csv")
d3 <- readData("Data/4.3.csv")
d41 <- readData("Data/4.4.1.csv")
d42 <- readData("Data/4.4.2.csv")

d51_raw <- readData("Data/GvHD+.csv", change_names=FALSE)
names(d51_raw) <- c("CD4", "CD8beta", "CD3", "CD8")
d51 <- d51_raw[d51_raw$CD3 > 280, ]



d52_raw <- readData("Data/GvHD-.csv", change_names=FALSE)
names(d52_raw) <- c("CD4", "CD8beta", "CD3", "CD8")
d52 <- d52_raw[d52_raw$CD3 > 280, ]

K_min <- 1
K_max <- 10

mc_BIC <- mclustBIC(d1)
mc <- Mclust(d1, x=mc_BIC)

mc_params <- mc %>% getModelParameters


