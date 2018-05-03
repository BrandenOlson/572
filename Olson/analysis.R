library(strucchange)

source("cluster.R")
source("EM.R")
source("processData.R")

runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_min=1,
                        K_max=12,
                        contour_levels=0.01
                        ) {

    pdf(paste0(output_dir, "unclustered.pdf"))
    plot(working_dat, pch=19, asp=1)
    dev.off()

    mc <- getModel(working_dat,
                   G=K_min:K_max,
                   model_names)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    # densities <- runEM(K, working_dat)
    K_BIC <- densities %>% length
    ts <- getTs(dat=working_dat, components=densities, K=K_BIC)
    zs <- getZs(ts)
    plot(working_dat, col=zs, pch=19)
    
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
    plotDensities(icl_densities,
                  output_prefix=paste0(output_dir, "ICL"),
                  dat=working_dat,
                  zs=zs_icl,
                  contour_levels=contour_levels
                  )

    plot_names <- paste0(output_dir, "merged_", K_BIC:1)
    mapply(function(x, y) {
               zs <- getTs(dat=working_dat,
                           components=x,
                           K=K) %>%
                   getZs
               plotDensities(x,
                             output_prefix=y,
                             dat=working_dat,
                             zs=zs,
                             contour_levels=contour_levels
                             )
           },
           cluster_seq,
           plot_names
    )
    
    
    pdf(paste0(output_dir, "entropy.pdf"))
    par(mfrow=c(2, 2))
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$K,
                           ylab="Difference in entropy",
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

plotPiecewiseFitToData <- function(
                                   y_col,
                                   x_col,
                                   xlab,
                                   ylab,
                                   h=3
                                   ) {
    plot(y_col ~ x_col, pch=19, xlab=xlab, ylab=ylab)
    if(FALSE) {
        # Need to find a way to fit piecewise linear regression with few data points
        bp <- breakpoints(y_col ~ x_col, h=h)
        # seg_fit <- lm(y_col ~ x_col*(x_col < bp) + x_col*(x_col >= bp))
        lines(fitted(bp, breaks=1) ~ xcol, lty=2, col="red")
    }
}
                                  
runAnalysis(unif_dat[, 1:2],
            model_names="VII",
            output_dir="Figures/unif/",
            contour_levels=1/10)
runAnalysis(d1,
            output_dir="Figures/d1/")
runAnalysis(d2,
            model_names=c("EEI", "VEI", "EVI", "VVI"),
            output_dir="Figures/d2/")
runAnalysis(d3,
            model_names=c("EII", "VII"),
            output_dir="Figures/d3/")
runAnalysis(d41,
            output_dir="Figures/d41/",
            contour_levels=0.001
            )
runAnalysis(d_square_circle,
            model_names="VII",
            output_dir="Figures/square_circle/",
            contour_levels=1/10)
# Need to make a plot just for the cluster labels,
# since it's difficult to plot 3D contours
runAnalysis(d42,
            output_dir="Figures/d42/"
            )
runAnalysis(d51_raw,
            output_dir="Figures/d51_raw/"
            )
runAnalysis(d52,
            output_dir="Figures/d52_raw/"
            )

stop()

# Proof that the paper and mclust give different clusterings
mc <- mclustBIC(d1, G=6, modelNames="VVV") %>% Mclust(data=d1, x=.)
