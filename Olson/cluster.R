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

getConditionalProbability <- function(x_i, params) {
    prop <- params$Prop
    means <- params$Mean
    variances <- params$Variance
    K <- length(prop)
    probs <- {}
    for(k in 1:K) {
        probs[k] <- prop[k]*dmvnorm(x_i,
                                     mean=means[, k],
                                     sigma=variances[, , k]
                                    ) 
    }
    t_ik <- probs/sum(probs)
    return(t_ik)
}

getMAP <- function(t_ik) {
    z_i <- t_ik %>% which.max
    return(z_i)
}

getTs <- function(dat, params, K) {
    ts <- {}
    for(i in 1:nrow(dat)) {
        ts[[i]] <- getConditionalProbability(dat[i, ],
                                             params)
    }
    return(ts)
}

getZs <- function(ts) {
    zs <- ts %>% sapply(getMAP)
    return(zs)
}

BIC <- function(log_lik, nu, n) {
    bic <- log_lik - nu*log(n)/2
    return(bic)
}

getEntropy <- function(ts) {
    t_vec <- ts %>% unlist
    ent <- -sum(t_vec*log(t_vec))
    return(ent)
}

d1 <- readData("Data/4.1.csv")
d2 <- readData("Data/4.2.csv")
d3 <- readData("Data/4.3.csv")
d41 <- readData("Data/4.4.1.csv")
d42 <- readData("Data/4.4.2.csv")

d51_raw <- readData("Data/GvHD+.csv", change_names=FALSE)
names(d51_raw) <- c("CD4", "CD8beta", "CD3", "CD8")

d52_raw <- readData("Data/GvHD-.csv", change_names=FALSE)
names(d52_raw) <- c("CD4", "CD8beta", "CD3", "CD8")

K_min <- 1
K_max <- 10

mc_BIC <- mclustBIC(d1)
mc <- Mclust(d1, x=mc_BIC)
mc_params <- mc %>% getModelParameters
ts <- getTs(d1, mc_params)
zs <- getZs(ts)
plot(d1, col=zs, pch=19)
ent <- getEntropy(ts)

getDensities <- function(params) {
    K <- length(params$Prop)
    fs <- vector(mode="list",
                 length=K)
    for(k in 1:K) {
        fs[[k]] <- list(Props=list(params$Prop[k]),
                        Mean=list(params$Mean[, k]),
                        Var=list(params$Variance[, , k])
                        )
    }
    return(fs)
}

combineDensities <- function(f1, f2) {
    f <- list(Props=c(f1$Prop,
                      f2$Prop),
              Mean=c(f1$Mean,
                     f2$Mean),
              Var=c(f1$Var,
                    f2$Var)
              )
    return(f)
}

getDensityFunction <- function(f) {
    num_components <- length(f$Props)
    d <- function(x) {
        res <- 0
        for(i in 1:num_components) {
            res <- res + f$Props[[i]]*dmvnorm(x, 
                                              mean=f$Mean[[i]],
                                              sigma=f$Var[[i]])
        }
        return(res)
    }
    return(d)
}

f <- combineDensities(fs[[1]], fs[[2]])
f <- fs[[1]] %>% combineDensities(fs[[2]]) %>%
    combineDensities(fs[[3]]) %>%
    combineDensities(fs[[4]]) %>%
    combineDensities(fs[[5]]) %>%
    combineDensities(fs[[6]])
d <- getDensityFunction(f)
m <- 100
xs <- seq(-4, 10, length.out=m)
ys <- seq(-4, 10, length.out=m)
z <- matrix(0, nrow=m, ncol=m)
for(i in 1:m) for(j in 1:m) z[i, j] <- d(c(xs[i], ys[j]))

surface3d(xs, ys, z, col="lightgreen")
aspect3d(1, 1, 1)
bbox3d(back="lines")
