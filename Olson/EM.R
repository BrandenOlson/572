source("cluster.R")

createComponents <- function(p_init,
                             mu_init,
                             Sigma_init) {

    components <- {}
    for(k in 1:K) {
        components[[k]] <- list(
                                Props=list(p_init[k]),
                                Mean=list(mu_init[[k]]),
                                Var=list(Sigma_init[[k]])
                               ) 

    }
    return(components)
}

generateRandomMeans <- function(K, dat) {
    mu_lower <- apply(dat, 2, min)
    mu_upper <- apply(dat, 2, max)
    means <- replicate(K,
                            mapply(function(x, y) { 
                                       runif(1, min=x, max=y) 
                                   },
                                   mu_lower,
                                   mu_upper),
                       simplify=FALSE
                      )
    return(means)
}

initializeEM <- function(K,
                         dat,
                         trial_count=1
                         ) {
    print("Initializing the parameters for EM...")
    liks <- {}
    d <- ncol(dat)
    components <- {}
    for(i in 1:trial_count) {
        print(i)
        p_init <- rep(1/K, K)
        mu_init <- rep(NA, K) %>% as.list
        Sigma_init <- {}
        mu_init <- kMeans(K, dat)
        for(k in 1:K) {
            Sigma_init[[k]] <- diag(rep(1, d))
        }
        components[[i]] <- createComponents(p_init,
                                            mu_init,
                                            Sigma_init)
        liks[i] <- getModelLikelihood(dat,
                                  components[[i]])
    }
    return(components[[which.max(liks)]])
}

runEM <- function(K, 
                  dat,
                  components,
                  tol=1e-4,
                  max_iterations=50
                  ) {
    print("Running the EM algorithm...")
    dat <- dat %>% as.matrix

    n <- nrow(dat)
    likelihood <- getModelLikelihood(dat, components)
    error <- Inf
    iteration <- 1
    while(error > tol && iteration <= max_iterations) {
        cat(error, likelihood, '\n')
        likelihood_prev <- likelihood
        tik <- getTs(dat, 
                     components,
                     K
                     ) %>%
             as.data.frame
        ps <- tik %>% apply(1, mean)
        ti_sums <- tik %>% apply(1, sum)

        for(k in 1:K) {
            components[[k]]$Props <- list(ps[k])
            tk <- tik[k, ] %>% as.matrix
            mu_k <- 0
            mu_k <- tk %*% dat
            components[[k]]$Mean <- list(mu_k/ti_sums[k])

            Sigma_k <- 0
            Sigma_k <- dat %>% 
                sweep(2, components[[k]]$Mean[[1]]) %>% # Subtract mu_k fro dat
                plyr::alply(1, tcrossprod) %>% # Compute outer product of each row
                Map("*", ., tk) %>% # Multiply each matrix the t_ik
                Reduce("+", .) # Sum matrices

            components[[k]]$Var <- list(Sigma_k/ti_sums[k])
        }

        likelihood <- getModelLikelihood(dat, components)
        error <- (likelihood - likelihood_prev)/abs(likelihood_prev)
        components %>% plotDensities(dat=dat)
        iteration <- iteration + 1
    }
    return(components)
}

kMeans <- function(K, 
                   dat,
                   tol=1e-3
                   ) {
    means <- generateRandomMeans(K, dat)
    theta_norm <- norm(means %>%
                       unlist %>%
                       as.matrix,
                   "F")
    error <- Inf
    while(error > tol) {
        theta_norm_prev <- theta_norm
        rs <- dat %>% apply(1,
                      function(x) {
                          distances <- {}
                          for(k in 1:K) {
                              diff <- as.matrix(x - means[[k]])
                               distances[k] <- norm( diff^2, "F")
                          }   
                          return(which.min(distances))
                      }
                )

        for(k in 1:K) {
            means[[k]] <- dat[rs == k, ] %>% apply(2, mean)
        }
        theta_norm <- means %>%
            unlist %>%
            as.matrix %>%
            norm("F")
        error <- abs(theta_norm - theta_norm_prev)/abs(theta_norm_prev)
    }
    return(means)
}


K <- 6
p_init <- rep(1/K, K)
mu_init <- list(c(0, 0), c(-1, 5), c(1, 5), c(8, 1), c(8, -1), c(8,5))
Sigma_init <- list(diag(rep(1, 2)),
                   diag(rep(1, 2)),
                   diag(rep(1, 2)),
                   diag(rep(1, 2)),
                   diag(rep(1, 2)),
                   diag(rep(1, 2)))

kmeans <- kMeans(K, d1)

comp_init <- initializeEM(K,
                          d1)

comp_init %>% plotDensities(dat=d1)

theta_mle <- runEM(K=K,
                   dat=d1,
                   components=comp_init
                  )

theta_mle %>% plotDensities(dat=d1)
