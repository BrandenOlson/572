library(strucchange)

source("cluster.R")
source("EM.R")
source("processData.R")

runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_min=1,
                        K_max=12,
                        contour_levels=0.01,
                        plot_types="2d"
                       ) {


    mc <- getModel(working_dat,
                   G=K_min:K_max,
                   model_names)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    K_BIC <- densities %>% length
    ts <- getTs(dat=working_dat, components=densities, K=K_BIC)
    zs <- getZs(ts)
    
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
    zs_icl <- getTs(dat=working_dat,
                    components=icl_densities,
                    K=K_ICL
                    ) %>%
        getZs

    # Plot unlabeled data
    pdf(paste0(output_dir, "unclustered.pdf"))
    plot(working_dat, pch=19, asp=1, cex.axis=1.5, cex.lab=2)
    dev.off()

    ICL_prefix <- paste0(output_dir, "ICL")

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
    plot_names <- paste0(output_dir, "merged_", K_BIC:1)
    mapply(function(x, y) {
               zs <- getTs(dat=working_dat,
                           components=x,
                           K=K) %>%
                   getZs
               if(plot_types=="2d") {
                   plotDensities(x,
                             output_prefix=y,
                             dat=working_dat,
                             zs=zs,
                             contour_levels=contour_levels
                             )
               } else if(plot_types=="3d") {
                    plotDensities3D(dat=working_dat,
                                    prefix=y,
                                    zs=zs)
               } else if(plot_types=="4d") {
                    plotCD3Clusters(dat=working_dat,
                                    components=x,
                                    zs=zs,
                                    prefix=y)
               }
           },
           cluster_seq,
           plot_names
    )
    
    pdf(paste0(output_dir, "entropy.pdf"))
    par(mfrow=c(2, 2))
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$K - 1,
                           ylab="Difference in entropy from K + 1",
                           xlab="K")
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$K,
                           xlab="K",
                           ylab="Total entropy"
                          )
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$NumMergedCumSum,
                           xlab="Cumulative count of merged observations",
                           ylab="Difference in entropy")
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
                                   h=3
                                   ) {
    valid_indices <- which(!is.na(y_col))
    ys <- y_col[valid_indices]
    xs <- x_col[valid_indices]
    n <- length(xs)
    yrange <- c(0, max(ys))
    plot(ys ~ xs, pch=19, xlab=xlab, ylab=ylab, ylim=yrange)
    axis(1, at=xs, labels=xs)
    a1 <- b1 <- a2 <- b2 <- NA
    lse <- NA
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
                                  
runAnalysis(d1,
            output_dir="Figures/d1/")
runAnalysis(d2,
            K_max=15,
            model_names=c("EEI", "VEI", "EVI", "VVI"),
            output_dir="Figures/d2/")
runAnalysis(d3,
            model_names=c("VII"),
            output_dir="Figures/d3/")
runAnalysis(d41,
            output_dir="Figures/d41/",
            contour_levels=0.001
            )
runAnalysis(unif_dat[, 1:2],
            model_names="VII",
            output_dir="Figures/unif/",
            contour_levels=1/10)
stop()
runAnalysis(d51_raw,
            output_dir="Figures/d51/",
            plot_types="4d",
            K_max=15,
            model_names="VVV"
            )
runAnalysis(d52_raw,
            output_dir="Figures/d52/",
            model_names="VVV",
            K_max=15,
            plot_types="4d"
            )
runAnalysis(d42,
            output_dir="Figures/d42/",
            plot_types="3d"
            )

# Need to make a plot just for the cluster labels,
# since it's difficult to plot 3D contours

stop()
