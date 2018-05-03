library(data.table)
library(dplyr)
library(magrittr)
library(mclust)
library(mvnfast)
library(plyr)

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
                     model_names=NULL,
                     criterion="BIC"
                     ) {
    if(criterion == "BIC") {
        mc_raw <- mclustBIC(dat,
                            G=G,
                            modelNames=model_names)
        mc <- Mclust(dat, x=mc_raw)
    } else if(criterion == "ICL") {
        mc_raw <- mclustICL(dat,
                            G=G,
                            modelNames=model_names
                            )

        # Mclust tries REALLY hard to 
        K <- mc_raw %>%
            summary %>%
            attributes %$%
            names %>%
            first %>%
            strsplit(",") %>%
            unlist %>%
            last %>%
            as.numeric
        mc <- Mclust(dat, G=K, modelNames=model_names)
    }
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
    log_sum <- dat %>% 
        apply(1,
              function(x) {
                  res <- components %>% 
                      sapply(function(y) { 
                                y$Prop[[1]]*dmvn(
                                    x, 
                                    mu=y$Mean[[1]],
                                    sigma=y$Var[[1]]
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
                          zs,
                          num_points=100,
                          contour_levels=0.01
                          ) {
    colors <- rainbow(length(density_list))
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
        pdf(paste0(output_prefix, "_contour.pdf"), width=10, height=10)
    }
    
    color <- 1
    plot(dat, xlim=xrange, ylim=yrange, col=colors[zs], pch=19)
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
                col=colors[color], 
                levels=contour_levels,
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
    return(list(argmax=argmax,
                delta_ent=max_delta
           )
    )
}

mergeClusters <- function(components, 
                          ts,
                          dat
                          ) {
    K <- length(components)
    delta_ent_object <- argmaxDelta(ts, K)
    argmax <- delta_ent_object$argmax
    delta_ent <- delta_ent_object$delta_ent
    zs <- getTs(dat=dat, components=components, K=k) %>%
        getZs
    num_merged <- dat[zs %in% argmax, ] %>% nrow
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
    total_ent <- 
    return(list(components=new_components,
                delta_ent=delta_ent,
                num_merged=num_merged
                ))
}

# Get cluster sequence along with the changes in entropy for each K-value
getClusterSequence <- function(dat,
                               components
                               ) {
    cluster_list <- list(components)
    delta_ents <- NA
    num_merged <- NA
    K <- length(components)
    total_ents <- getTs(dat=dat,
                        components=components,
                        K=K) %>%
        getEntropy
    k <- K
    while(k > 1) {
        ts <- getTs(dat=dat, components=components, K=k)
        merge_clusters_object <- mergeClusters(components, ts, dat=dat)
        components <- merge_clusters_object$components
        delta_ent <- merge_clusters_object$delta_ent
        num_merged_k <- merge_clusters_object$num_merged

        delta_ents <- c(delta_ents, delta_ent)
        total_ent <- getTs(dat=dat,
                           components=components,
                           K=k) %>%
            getEntropy
        total_ents <- c(total_ents, total_ent)
        num_merged <- c(num_merged, num_merged_k)
        cluster_list <- c(cluster_list,
                          list(components)
                          )
        k <- k - 1
    }

    # Baudry plots the cum sums in reverse for some reason...
    num_merged_cumsum <- num_merged[!is.na(num_merged)] %>%
        cumsum %>%
        append(0, .) %>%
        rev


    return(
        list(cluster_list=cluster_list,
            entropy=data.frame(K=K:1,
                               DeltaEntropy=delta_ents,
                               TotalEntropy=total_ents,
                               NumMerged=num_merged,
                               NumMergedCumSum=num_merged_cumsum,
                               NormalizedDiff=delta_ents/num_merged
                              )
                )
    )
}

