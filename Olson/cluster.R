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

getComponentProbability <- function(x_i, density_object) {
    d <- density_object %>% getDensityFunction
    prob <- d(x_i)
    return(prob)
}

getConditionalProbability <- function(x_i, components) {
    K <- length(components)
    probs <- {}
    for(k in 1:K) {
        component_density <- components[[k]]
        probs[k] <- getComponentProbability(x_i,
                                            component_density)
    }
    t_ik <- probs/sum(probs)
    return(t_ik)
}

getMAP <- function(t_ik) {
    z_i <- t_ik %>% which.max
    return(z_i)
}

getTs <- function(dat, components, K) {
    ts <- {}
    for(i in 1:nrow(dat)) {
        ts[[i]] <- getConditionalProbability(dat[i, ],
                                             components)
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
    return(list(f))
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

plotDensities <- function(density_list,
                          output_prefix,
                          xrange,
                          yrange,
                          dat,
                          num_points=100
                          ) {
    color <- 2 # Skip black (1) since points are black 
    pdf(paste0(output_prefix, "_contour.pdf"), width=10, height=6)
    plot(dat, pch=19)
    for(density_object in density_list) {
        d <- density_object %>%
            getDensityFunction
        xs <- seq(xrange[1], xrange[2], length.out=num_points)
        ys <- seq(yrange[1], yrange[2], length.out=num_points)
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
    dev.off()
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
        delta <- delta + t_union*log(t_union) - t_ij*log(t_ij) - t_ik*log(t_ik)
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
                    cat(delta, max_delta, '\n')
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

if(FALSE) {
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

    working_dat <- d1
    mc_BIC <- mclustBIC(working_dat)
    mc <- Mclust(working_dat, x=mc_BIC)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    ts <- getTs(working_dat, densities)
    zs <- getZs(ts)
    plot(working_dat, col=zs, pch=19)
    ent <- getEntropy(ts)

    fs <- mc_params %>% 
        getDensities
    f <- combineDensities(fs[[1]], fs[[2]])

    plotDensities(list(f),
                output_prefix="Figures/d1_combined",
                xrange=range(working_dat[, 1]),
                yrange=range(working_dat[, 2]),
                dat=working_dat
                )
    
    plotDensities(fs,
                output_prefix="Figures/d1_BIC",
                xrange=range(working_dat[, 1]),
                yrange=range(working_dat[, 2]),
                dat=working_dat
                )
}

merged <- mergeClusters(fs, ts)
t_merged <- getTs(working_dat, merged)
plotDensities(merged,
              output_prefix="Figures/merged",
              xrange=range(working_dat[, 1]),
              yrange=range(working_dat[, 2]),
              dat=working_dat)

merged_2 <- mergeClusters(merged, t_merged)
plotDensities(merged_2,
              output_prefix="Figures/merged_2",
              xrange=range(working_dat[, 1]),
              yrange=range(working_dat[, 2]),
              dat=working_dat)
