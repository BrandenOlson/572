library(data.table)
library(dplyr)
library(magrittr)
library(mclust)
library(mvnfast)

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

getModel <- function(dat,
                     G,
                     model_names=NULL) {
    mc_BIC <- mclustBIC(dat,
                        G=G,
                        modelNames=model_names)
    mc <- Mclust(dat, x=mc_BIC)
    return(mc)
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

getModelLikelihood <- function(dat, components) {
    ell <- 0
    log_sum <- dat %>% 
        apply(1,
              function(x) {
                  res <- components %>% 
                      sapply(function(y) { 
                                y$Prop[[1]]*dmvn(
                                    x, 
                                    mu=y$Mean[[1]],
                                    sigma=y$Var[[1]],
                                )
                             }
                            ) %>%
                      sum
                  return(res)
               }
        ) %>%
        log %>%
        sum 
    return(log_sum)
}

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

getDensityFunction <- function(f) {
    num_components <- length(f$Props)
    d <- function(x) {
        res <- 0
        for(i in 1:num_components) {
            res <- res + f$Props[[i]]*dmvn(x %>% as.numeric, 
                                              mu=f$Mean[[i]],
                                              sigma=f$Var[[i]])
        }
        return(res)
    }
    return(d)
}

getComponentProbability <- function(x_i, density_object) {
    d <- density_object %>% getDensityFunction
    prob <- d(x_i)
    return(prob)
}

getConditionalProbability <- function(x_i, components) {
    K <- length(components)
    probs <- components %>%
        sapply(function(c) { getComponentProbability(x_i,
                                                     c)
                      })
    t_ik <- probs/sum(probs)
    return(t_ik)
}

getMAP <- function(t_ik) {
    z_i <- t_ik %>% which.max
    return(z_i)
}

getTs <- function(dat, components, K) {
    ts <- dat %>% 
        plyr::alply(1, getConditionalProbability, components=components)
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

combineDensities <- function(f1, f2) {
    f <- list(Props=c(f1$Prop,
                      f2$Prop),
              Mean=c(f1$Mean,
                     f2$Mean),
              Var=c(f1$Var,
                    f2$Var)
              )
    return(list(f))
}

plotDensities <- function(density_list,
                          output_prefix,
                          dat,
                          num_points=100
                          ) {
    color <- 2 # Skip black (1) since points are black 
    xrange <- range(dat[, 1])
    yrange <- range(dat[, 2])
    xrange <- c(xrange[1] - length(xrange)*0.2,
              xrange[2] + length(xrange)*0.2)
    yrange <- c(yrange[1] - length(yrange)*0.2,
              yrange[2] + length(yrange)*0.2)
    xs <- seq(xrange[1], xrange[2], length.out=num_points)
    ys <- seq(yrange[1], yrange[2], length.out=num_points)
    if(!missing(output_prefix)) {
        pdf(paste0(output_prefix, "_contour.pdf"), width=10, height=6)
    }
    plot(dat, xlim=xrange, ylim=yrange, pch=19)
    for(density_object in density_list) {
        d <- density_object %>%
            getDensityFunction
        z <- matrix(0, nrow=num_points, ncol=num_points)
        for(i in 1:num_points) {
            for(j in 1:num_points) {
                z[i, j] <- d(c(xs[i], ys[j]))
            }
        }

        if(FALSE) {
            pdf(paste0(output_prefix, "_density.pdf"), width=10, height=6)
            surface3d(xs, ys, z, col="lightgreen")
            aspect3d(1, 1, 1)
            bbox3d(back="lines")
            dev.off()
        }

        contour(x=xs, y=ys, z, 
                col=color, 
                levels=0.01, 
                add=TRUE,
                drawlabels=FALSE)
        color <- color + 1
    }
    if(!missing(output_prefix)) {
        dev.off()
    }
}

# Helper function to check if t*log(t) is a valid number, and if not,
# set it to zero, since lim_{t -> 0} t*log(t) = 0
getPointEntropy <- function(t) {
    entropy <- ifelse(!is.nan(t*log(t)),
                      t*log(t),
                      0)
    return(entropy)
}

deltaEntropy <- function(j,
                         k,
                         ts
                         ){
    delta <- 0
    for(i in 1:length(ts)) {
        t_ij <- ts[[i]][j]
        t_ik <- ts[[i]][k]
        t_union <- t_ij + t_ik
        entropy_union <- t_union %>% getPointEntropy
        entropy_ij <- t_ij %>% getPointEntropy
        entropy_ik <- t_ik %>% getPointEntropy

        delta <- delta + entropy_union - entropy_ij - entropy_ik
    } 
    return(delta)
}

argmaxDelta <- function(ts, K) {
    argmax <- {}
    max_delta <- -Inf
    for(j in 1:K) {
        for(k in 1:K) {
            if(j != k) {
                delta <- deltaEntropy(j, k, ts)
                if(delta > max_delta) {
                    max_delta <- delta
                    argmax <- c(j, k) 
                }
            }
        }
    }
    return(argmax)
}

mergeClusters <- function(components, ts) {
    K <- length(components)
    argmax <- argmaxDelta(ts, K)
    new_components <- list()
    merged <- FALSE
    for(i in 1:K) {
        if(i %in% argmax) {
            if(!merged) {
                new_components <- c(new_components,
                                    combineDensities(components[[argmax[1]]],
                                                     components[[argmax[2]]])
                                    )
                merged <- TRUE
            }
        }
        else {
            new_components <- c(new_components,
                                list(components[[i]]))
        }
    }
    return(new_components)
}

getClusterSequence <- function(dat,
                               components) {
    cluster_list <- list(components)
    K <- length(components)
    while(K > 1) {
        ts <- getTs(dat, components)
        components <- mergeClusters(components, ts)
        cluster_list <- c(cluster_list,
                          list(components)
                          )
        K <- K - 1
    }
    return(cluster_list)
}

