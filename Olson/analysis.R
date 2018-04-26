source("cluster.R")
source("EM.R")
source("processData.R")

runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_min=1,
                        K_max=12
                        ) {

    mc <- getModel(working_dat,
                   G=K_min:K_max,
                   model_names)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    # mc_params <- runEM(K, working_dat)
    K_BIC <- densities %>% length
    ts <- getTs(working_dat, densities)
    zs <- getZs(ts)
    plot(working_dat, col=zs, pch=19)
    
    fs <- mc_params %>% 
        getDensities
    f <- combineDensities(fs[[1]], fs[[2]])
    
    cluster_seq <- getClusterSequence(working_dat,
                                      fs)
    
    if(!file.exists(output_dir)) {
        dir.create(output_dir)
    }

    plot_names <- paste0(output_dir, "merged_", K_BIC:1)
    mapply(function(x, y) {
               plotDensities(x,
                             output_prefix=y,
                             dat=working_dat)
           },
           cluster_seq,
           plot_names
    )
}

runAnalysis(d1,
            output_dir="Figures/d1/")
stop()
runAnalysis(d2,
            model_names="VVI",
            output_dir="Figures/d2/")
runAnalysis(d3,
            model_names="VII",
            output_dir="Figures/d3/",
            K_min=4)
