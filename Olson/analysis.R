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
    ts <- getTs(working_dat, densities)
    zs <- getZs(ts)
    plot(working_dat, col=zs, pch=19)
    
    cluster_seq <- getClusterSequence(working_dat,
                                      densities)

    if(!file.exists(output_dir)) {
        dir.create(output_dir)
    }

    plot_names <- paste0(output_dir, "merged_", K_BIC:1)
    mapply(function(x, y) {
               zs <- getTs(working_dat,
                           x,
                           K) %>%
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
}

runAnalysis(d42,
            output_dir="Figures/d42/"
            )
runAnalysis(d1,
            output_dir="Figures/d1/")
runAnalysis(d2,
            model_names=c("EEI", "VEI", "EVI", "VVI"),
            output_dir="Figures/d2/")
runAnalysis(d3,
            model_names=c("EII", "VII"),
            output_dir="Figures/d3/")

stop()
# Need to make a plot just for the cluster labels,
# since it's difficult to plot 3D contours
runAnalysis(d41,
            output_dir="Figures/d41/",
            contour_levels=0.001
            )
runAnalysis(d51_raw,
            output_dir="Figures/d51_raw/"
            )
runAnalysis(d52,
            output_dir="Figures/d52_raw/"
            )

