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
    
    
    pdf(paste0(output_dir, "entropy.pdf"), width=10, height=6)
    par(mfrow=c(2, 2))
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$K)
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$K)
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$NumMergedCumSum)
    plotPiecewiseFitToData(entropy$NormalizedDiff,
                           entropy$K)
    dev.off()
}

plotPiecewiseFitToData <- function(
                                   x_col,
                                   y_col,
                                   breakpoint
                                   ) {
    bp <- breakpoints(y_col ~ x_col, h=3)
    seg_fit <- lm(y_col ~ x_col*(x_col < bp) + x_col*(x_col >= bp))
    plot(y_col ~ x_col, data=dat, pch=19)
    lines(fitted(bp, breaks=1) ~ xcol, lty=2, col="red")
}
                                  

unif_dat <- # rbind(cbind(runif(400, 0, 1), runif(200, 0, 1)),
            #       rmvn(200, mu=c(1, 2), sigma=diag(1/5, 2))
            #       )
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
            contour_levels=1/10000)
stop()
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
