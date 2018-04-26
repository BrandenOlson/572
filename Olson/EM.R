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
                         trial_count=10
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
        error <- abs(likelihood_prev - likelihood)/abs(likelihood_prev)
        components %>% plotDensities(dat=dat)
        iteration <- iteration + 1
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

comp_init %>% plotDensities(dat=d1)

theta_mle <- runEM(K=K,
                   dat=d1,
                   components=comp_init
                  )

theta_mle %>% plotDensities(dat=d1)
