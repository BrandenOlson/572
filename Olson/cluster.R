library(data.table)
library(dplyr)
library(magrittr)
library(mclust)
library(mvnfast)
library(plyr)
library(scatterplot3d)
library(strucchange)
library(xtable)

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
                        Var=list(params$Variance[, , k]),
                        Label=k # Avoid k = 1, which plots black points
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
    argmax <- t_ik %>% which.max
    return(argmax)
}

getTs <- function(dat, components, K) {
    ts <- dat %>% 
        plyr::alply(1, getConditionalProbability, components=components)
    return(ts)
}

getLabelFromIndex <- function(argmax, components) {
    label <- components[[argmax]]$Label
    return(label)
}

getZs <- function(ts, components) {
    zs <- ts %>% 
        sapply(getMAP) %>%
        sapply(getLabelFromIndex, components=components)

    return(zs)
}

BIC <- function(log_lik, nu, n) {
    bic <- log_lik - nu*log(n)/2
    return(bic)
}

getEntropy <- function(ts) {
    t_vec <- ts %>% unlist
    ent <- -sum(t_vec*log(t_vec)) %>% abs
    return(ent)
}

combineDensities <- function(f1, f2) {
    f <- list(Props=c(f1$Prop,
                      f2$Prop),
              Mean=c(f1$Mean,
                     f2$Mean),
              Var=c(f1$Var,
                    f2$Var),
              Label=min(f1$Label, f2$Label)
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
    xrange <- range(dat[, 1])
    yrange <- range(dat[, 2])
    xrange <- c(xrange[1] - length(xrange)*0.2,
              xrange[2] + length(xrange)*0.2)
    yrange <- c(yrange[1] - length(yrange)*0.2,
              yrange[2] + length(yrange)*0.2)
    xs <- seq(xrange[1], xrange[2], length.out=num_points)
    ys <- seq(yrange[1], yrange[2], length.out=num_points)
    if(!missing(output_prefix)) {
        pdf(paste0(output_prefix, "_contour.pdf"), width=6, height=6)
    }
    
    plot(dat, xlim=xrange, ylim=yrange, col=zs + 1, pch=19, asp=1,
         xlab='x', ylab='y',
         cex=0.5, cex.axis=2, cex.lab=2)
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
            pdf(paste0(output_prefix, "_density.pdf"), width=10, height=10)
            surface3d(xs, ys, z, col="lightgreen")
            aspect3d(1, 1, 1)
            bbox3d(back="lines")
            dev.off()
        }

        contour(x=xs, y=ys, z, 
                col=density_object$Label + 1, 
                levels=contour_levels,
                add=TRUE,
                drawlabels=FALSE)
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

argmaxDelta <- function(components, ts, K) {
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
    delta_ent_object <- argmaxDelta(components, ts, K)
    argmax <- delta_ent_object$argmax
    argmax_labels <- argmax %>% sapply(getLabelFromIndex,
                                       components=components)
    delta_ent <- delta_ent_object$delta_ent
    zs <- getTs(dat=dat, components=components, K=k) %>%
        getZs(components=components)


    num_merged <- zs[zs %in% argmax_labels] %>% length
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
    delta_ents <- {}
    num_merged <- {}
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

    delta_ents <- c(delta_ents, NA)
    num_merged <- c(num_merged, 0)

    num_merged_cumsum <- num_merged %>%
        rev %>%
        cumsum %>%
        rev

    normalized_diff <- c((delta_ents/num_merged)[1:(K - 1)], NA)

    return(
        list(cluster_list=cluster_list,
            entropy=data.frame(K=K:1,
                               DeltaEntropy=delta_ents,
                               TotalEntropy=total_ents,
                               NumMerged=num_merged,
                               NumMergedCumSum=num_merged_cumsum,
                               NormalizedDiff=normalized_diff
                              )
                )
    )
}

# Basically specialized for the "3D uniform cross" example
plotDensities3D <- function(dat, prefix, zs) {
    pdf(paste0(prefix, "_scatter.pdf"))
    scatterplot3d(dat, pch=19, asp=1, color=zs)
    dev.off()

    pdf(paste0(prefix, "_top_2d.pdf"))
    plot(dat[, c(1, 2)], pch=19, asp=1, ylim=c(1, 0), col=zs)
    dev.off()
}

getComponentMean <- function(component) {
    mean <- mapply(function(x, y) { x*y },
                   component$Prop,
                   component$Mean
                   ) %>%
        c
    mean <- mean/sum(component$Prop %>% unlist)
    return(mean)
}

isCD3Positive <- function(component,
                          threshold=280
                          ) {
    component_mean <- getComponentMean(component)
    CD3_index <- 3
    is_cd3_pos <- component_mean[CD3_index] > threshold
    return(is_cd3_pos)
}

plotCD3Clusters <- function(dat,
                            components,
                            zs,
                            prefix
                            ) {
    cd3_pos <- components %>%
        sapply(isCD3Positive) %>%
        which %>%
        sapply(getLabelFromIndex, components=components)

    if(length(cd3_pos) > 0) {
        cd3_pos_dat <- dat[zs %in% cd3_pos, ]
        cd3_zs <- zs[zs %in% cd3_pos]
        pdf(paste0(prefix, "_CD3_CD8beta.pdf"), width=6, height=6)
        par(mar=c(5, 5, 2, 2))
        plot(CD8beta ~ CD4, 
             data=cd3_pos_dat, 
             col=cd3_zs %>% as.factor,
             pch=19, 
             cex=0.3, cex.axis=2, cex.lab=2,
             asp=1)
        abline(h=280)
        abline(v=280)
        dev.off()
    }
}

runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_min=1,
                        K_max=12,
                        contour_levels=0.01,
                        plot_types="2d",
                        init_params=NULL
                       ) {
    # Plot unlabeled data
    point_size <- ifelse(plot_types == "4d", 0.1, 0.5)
    pdf(paste0(output_dir, "unclustered.pdf"), width=6, height=6)
    plot(working_dat, pch=19, asp=1, cex=point_size, cex.axis=2, cex.lab=2, xlab="x", ylab="y")
    dev.off()

    if(missing(init_params)) {
        mc <- getModel(working_dat,
                       G=K_min:K_max,
                       model_names)
        init_params <- mc %>% getModelParameters
    }
    densities <- init_params %>% getDensities
    K_BIC <- densities %>% length
    ts <- getTs(dat=working_dat, components=densities, K=K_BIC)
    zs <- getZs(ts, densities)
    entropy_bic <- ts %>% getEntropy
    
    cluster_seq_object <- getClusterSequence(working_dat,
                                      densities)
    cluster_seq <- cluster_seq_object$cluster_list
    entropy <- cluster_seq_object$entropy

    if(!file.exists(output_dir)) {
        dir.create(output_dir)
    }

    mc_ICL <- getModel(working_dat,
                       G=K_min:K_max,
                       model_names,
                       criterion="ICL"
                       )
    mc_icl_params <- mc_ICL %>% getModelParameters
    icl_densities <- mc_icl_params %>%
        getDensities
    K_ICL <- icl_densities %>% length
    ts_icl <- getTs(dat=working_dat,
                    components=icl_densities,
                    K=K_ICL)

    zs_icl <- ts_icl %>%
        getZs(components=icl_densities)
    entropy_icl <- ts_icl %>% getEntropy

    K_ent_df <- data.table(
                           Method=c("BIC", "ICL"),
                           K=c(K_BIC, K_ICL),
                           Entropy=c(entropy_bic, entropy_icl)
                          )

    ICL_prefix <- paste0(output_dir, "ICL", K_ICL)

    # Plot ICL solutions
    if(plot_types=="2d") {
        plotDensities(icl_densities,
                  output_prefix=ICL_prefix,
                  dat=working_dat,
                  zs=zs_icl,
                  contour_levels=contour_levels
                  )
    } else if(plot_types=="3d") {
        plotDensities3D(dat=working_dat,
                        prefix=ICL_prefix,
                        zs=zs_icl)
    } else if(plot_types=="4d") {
        plotCD3Clusters(dat=working_dat,
                        components=icl_densities,
                        zs=zs_icl,
                        prefix=ICL_prefix
                        )
    }

    # Plot all merged solutions
    k_ent_dfs <- Map(function(x, k) {
                         plot_name <- paste0(output_dir, "merged_", k)
                         ts <- getTs(dat=working_dat,
                                     components=x,
                                     K=K)
                         zs <- ts %>%
                             getZs(components=x)
                         entropy_k <- ts %>% getEntropy
                         k_ent_df <- data.table(Method="Merged",
                                                K=k,
                                                Entropy=entropy_k)

                         if(plot_types=="2d") {
                             plotDensities(x,
                                       output_prefix=plot_name,
                                       dat=working_dat,
                                       zs=zs,
                                       contour_levels=contour_levels
                                       )
                         } else if(plot_types=="3d") {
                              plotDensities3D(dat=working_dat,
                                              prefix=plot_name,
                                              zs=zs)
                         } else if(plot_types=="4d") {
                              plotCD3Clusters(dat=working_dat,
                                              components=x,
                                              zs=zs,
                                              prefix=plot_name)
                         }

                         return(k_ent_df)
                     },
                     cluster_seq,
                     K_BIC:1
                 ) %>% Reduce(rbind, .)

    K_ent_df <- rbind(K_ent_df,
                      k_ent_dfs)

    print(dim(K_ent_df))
    printTable(K_ent_df,
               paste0(output_dir, "K_ent.tex")
               )

    pdf(paste0(output_dir, "total_entropy.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$K,
                           xlab="K",
                           ylab="Total entropy",
                           do_regression=TRUE
                          )
    dev.off()
    pdf(paste0(output_dir, "diff_entropy.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$K - 1,
                           ylab="Difference in entropy from K + 1",
                           xlab="K")
    dev.off()
    pdf(paste0(output_dir, "cum_count.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$NumMergedCumSum,
                           xlab="Cumulative count of merged observations",
                           ylab="Total entropy",
                           do_regression=TRUE
                           )
    dev.off()
    pdf(paste0(output_dir, "norm_diff.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$NormalizedDiff,
                           entropy$K,
                           xlab="K",
                           ylab="Normalized difference in entropy")
    dev.off()
}

#' Follow the strategy of Baudry by iterating through each possible
#' change point in 2:(K_max - 1), computing a bi-piecewise regression,
#' and choosing the one that minimzes the total sum of squares
plotPiecewiseFitToData <- function(
                                   y_col,
                                   x_col,
                                   xlab,
                                   ylab,
                                   do_regression=FALSE,
                                   num_changepoints=1
                                   ) {
    valid_indices <- which(!is.na(y_col))
    ys <- y_col[valid_indices]
    xs <- x_col[valid_indices]
    n <- length(xs)
    yrange <- c(0, max(ys))
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    plot(ys ~ xs, pch=19, xlab=xlab, ylab=ylab, ylim=yrange, 
         cex=1.5, cex.axis=2, cex.lab=2)
    a1 <- b1 <- a2 <- b2 <- NA
    lse <- NA
    if(n > 2 && do_regression) {
        for(candidate in 2:(n - 1)) {
            x1 <- xs[1:candidate]
            x2 <- xs[candidate:n]
            y1 <- ys[1:candidate]
            y2 <- ys[candidate:n]

            a1[candidate] <- a1_c <- (sum(x1*y1) - sum(x1)*mean(y1))/
                (sum(x1^2) - sum(x1)^2/length(x1))
            b1[candidate] <- b1_c <- mean(y1) - a1[candidate]*mean(x1)

            a2[candidate] <- a2_c <- (sum(x2*y2) - sum(x2)*mean(y2))/
                (sum(x2^2) - sum(x2)^2/length(x2))
            b2[candidate] <- b2_c <- mean(y2) - a2[candidate]*mean(x2)

            lse[candidate] <- sum( (a1_c*x1 + b1_c - y1)^2 ) + 
                sum( (a2_c*x2 + b2_c - y2)^2 )
        }

        bp <- which.min(lse)
        lines(a1[bp]*xs[1:bp] + b1[bp] ~ xs[1:bp], col="red", lty=2)
        lines(a2[bp]*xs[bp:n] + b2[bp] ~ xs[bp:n], col="red", lty=2)
    }
    
}

printTable <- function(dat, filename, digits=3, hline_pos=0) {
    sink(filename)
    cat("\\begin{center}", '\n')
    dat %>%
        xtable::xtable(digits=digits) %>%
        print(
              floating=FALSE,
              include.rownames=FALSE,
              latex.environments="center",
              hline.after=c(0, hline_pos),
              sanitize.text.function=function(x){x}
             )
    cat("\\end{center}", '\n')
    sink()
}
