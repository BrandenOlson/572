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
    means <- dat %>%
        sample_n(K) %>%
        plyr::alply(1, as.numeric)
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
        mu_init <- rep(NA, K) %>% as.list
        kmeans_fit <- kMeans(K, dat)
        p_init <- rep(1/K, K) 
            # For some reason, using kmeans_fit$ps doesn't work as well
        mu_init <- kmeans_fit$means
        cov_init <- kmeans_fit$covariances

        components[[i]] <- createComponents(p_init,
                                            mu_init,
                                            cov_init)
        liks[i] <- getModelLikelihood(dat,
                                  components[[i]])
    }
    component_best <- components[[which.max(liks)]]
    zs <- getTs(dat, component_best, K) %>% getZs
    plotDensities(component_best, dat=dat, zs=zs)
    return(component_best)
}

runEM <- function(K, 
                  dat,
                  tol=1e-5,
                  max_iterations=50
                  ) {
    print("Running the EM algorithm...")
    components <- initializeEM(K, dat)
    dat <- dat %>% as.matrix

    likelihood <- getModelLikelihood(dat, components)
    error <- Inf
    iteration <- 1
    while((error > tol || iteration <= 3) && (iteration <= max_iterations)) {
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
        zs <- tik %>% getZs
        components %>% plotDensities(dat=dat, zs=zs)
        iteration <- iteration + 1
    }
    return(components)
}

kMeans <- function(K, 
                   dat,
                   tol=1e-1
                   ) {
    dat <- data.frame(dat)
    means <- dat %>%
        sample_n(K) %>%
        plyr::alply(1, as.numeric)
    covariances <- replicate(K,
                             diag(1, ncol(dat)),
                             simplify=F)
    ps <- rep(1/K, K)
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

        hcs <- dat %>%
            hc %>% 
            hclass(K)
        for(k in 1:K) {
            if(nrow(dat[hcs==k,])) {
                dat_k <- dat[hcs == k, ]
                means[[k]] <- dat_k %>% apply(2, mean)
                covariances[[k]] <- cov(dat_k)
                ps[k] <- nrow(dat_k)/nrow(dat)
                print(ps[k])
            }
        }

        theta_norm <- means %>%
            unlist %>%
            as.matrix %>%
            norm("F")
        error <- abs(theta_norm - theta_norm_prev)/abs(theta_norm_prev)
    }
    return(list(ps=ps, means=means, covariances=covariances))
}

classLikSolo <- function(xi, mean, cov) {
    if(!is.positive.definite(cov)) {
        cov <- diag(0.5, length(xi))
    }
    dmvn(xi,
         mean,
         cov,
         log=TRUE)
}

classLik <- function(dat, means, covs) {
    liks <- dat %>% 
        plyr::alply(1,
                    as.numeric) %>%
        mapply(classLikSolo,
               .,
               means,
               covs
               ) %>%
        unname
    return(liks)
}

getTotalLik <- function(class_liks) {
    class_liks %>% sum
}

mbhc <- function(K,
                 dat) {
    n <- nrow(dat)
    d <- ncol(dat)
    zs <- 1:n
    means <- dat %>%
        plyr::alply(1, as.numeric)
    covs <- replicate(n,
                      diag(1, d),
                      simplify=FALSE
                      )
    class_liks <- classLik(dat,
                           means,
                           covs)
    argmax <- 0
    for(K_curr in n:K) {
        clusters <- zs %>% unique
        means_current <- 0
        cov_current <- 0
        max_diff <- 0
        for(ci in clusters) {
            cat(ci, '\n')
            for(cj in clusters) {
                if(ci < cj) {
                    dat_ij <- dat[zs %in% c(ci, cj), ]
                    mu_ij <- dat_ij %>% apply(2, mean)
                    cov_ij <- cov(dat_ij)
                    lik_ij <- classLik(dat_ij, 
                                       replicate(nrow(dat_ij), mu_ij, simplify=FALSE),
                                       replicate(nrow(dat_ij), cov_ij, simplify=FALSE)
                                       ) %>% sum
                    lik_difference <- lik_ij - class_liks[ci] - class_liks[cj]
                    if(lik_difference > max_diff) {
                        argmax <- c(ci, cj)
                        means_current <- mu_ij
                        cov_current <- cov_ij
                    }
                } 
            }
        }
        zs[zs == argmax[2]] <- argmax[1]
        means[[argmax[1]]] <- means_current
        covs[[argmax[1]]] <- cov_current
    }

    clusters <- zs %>% unique
    means <- means[zs %>% unique]
    covs <- covs[zs %>% unique]
    return(list(means, covs))
}

# Try a model-based hierarchical clustering initialization

if(FALSE) {
K <- 6
# theta_init <- mbhc(K, dat=d1)
theta_mle <- runEM(K=K, dat=d1)
theta_mle %>% plotDensities(dat=d1)
}
