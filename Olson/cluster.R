library(data.table)
library(dplyr)
library(Rmixmod)

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

getModelParameters <- function(mixmod_object) {
    means <- mixmod_object %>%
        getSlot("mean")
    variances <- mixmod_object %>%
        getSlot("variance")
    return(list(Means=means,
                Variances=variances))
}

getModelLikelihood <- function(mixmod_object) {
    lik <- mixmod_object %>%
        slot("bestResult") %>%
        slot("likelihood")
    return(lik)
}

getBICFromModel <- function(mixmod_object) {
    bic <- mixmod_object %>%
        slot("bestResult") %>%
        slot("criterionValue")
    return(bic)
}

BIC <- function(log_lik, nu, n) {
    bic <- log_lik - nu*log(n)/2
    return(bic)
}

getKBIC <- function(Ks, models) {
    K_count <- length(Ks)
    bics <- {}
    for(i in 1:K_count) {
        K <- Ks[i]
        mod <- models[[i]]
        bics[i] <- getBICFromModel(mod)
    }  
    K_BIC <- Ks[which.min(bics)]
    print(bics)
    return(K_BIC)
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
K_max <- 6
strategy <- mixmodStrategy(
                           algo="EM",
                           nbIterationInInit=10,
                           nbIterationInAlgo=5000,
                           initMethod="random"
                          )
mixmod <- mixmodCluster(d1, 
                        K_min:K_max, 
                        strategy=strategy,
                        criterion="BIC"
                        )

K_BIC <- getKBIC(Ks, mixmod)
