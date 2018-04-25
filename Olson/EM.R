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

initializeEM <- function(K,
                         dat,
                         trial_count=20
                         ) {
    print("Initializing the parameters for EM...")
    liks <- {}
    mu_lower <- apply(d1, 2, min)
    mu_upper <- apply(d1, 2, max)
    d <- ncol(dat)
    components <- {}
    for(i in 1:trial_count) {
        print(i)
        p_init <- rep(1/K, K)
        mu_init <- rep(NA, K) %>% as.list
        Sigma_init <- {}
        for(k in 1:K) {
            mu_k <- mapply(function(x, y) { runif(1, min=x, max=y) },
                           mu_lower,
                           mu_upper)
            mu_init[[k]] <- mu_k 
            Sigma_init[[k]] <- diag(rep(1, d))
        }
        components[[i]] <- createComponents(p_init,
                                            mu_init,
                                            Sigma_init)
        liks[i] <- getModelLikelihood(dat,
                                  components[[i]])
    }
    print(liks)
    return(components[[which.max(liks)]])
}

runEM <- function(K, 
                  dat,
                  components,
                  tol=1e-4
                  ) {
    print("Running the EM algorithm...")
    dat <- dat %>% as.matrix

    n <- nrow(dat)
    theta <- getModelLikelihood(dat, components)
    error <- Inf
    while(error > tol) {
        print(error)
        theta_prev <- theta
        tik <- getTs(dat, 
                     components,
                     K
                     ) %>%
             as.data.frame
        ps <- tik %>% apply(1, mean)
        ti_sums <- tik %>% apply(1, sum)

        for(k in 1:K) {
            components[[k]]$Props <- list(ps[k])
            mu_k <- 0
            for(i in 1:n) {
                mu_k <- mu_k + dat[i, ]*tik[k, i]
            }
            components[[k]]$Mean <- list(mu_k/ti_sums[k])

            Sigma_k <- 0
            for(i in 1:n) {
                y_center <- (dat[i, ] - components[[k]]$Mean[[1]])
                Sigma_k <- Sigma_k + 
                    tik[k, i] * (y_center %*% t(y_center))
            }
            components[[k]]$Var <- list(Sigma_k/ti_sums[k])
        }

        theta <- getModelLikelihood(dat, components)
        error <- abs(theta_prev - theta)/abs(theta_prev)
    }
    return(components)
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

comp_init <- initializeEM(K,
                          d1)

theta_mle <- runEM(K=K,
                   dat=d1,
                   components=comp_init
                  )

comp_init %>% plotDensities(dat=d1,
                            output_prefix="random")
theta_mle %>% plotDensities(dat=d1,
                            output_prefix="test")
