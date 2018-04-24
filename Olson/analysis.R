runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_max=12
                        ) {

    mc_BIC <- mclustBIC(working_dat,
                        G=1:K_max,
                        modelNames=model_names)
    mc <- Mclust(working_dat, x=mc_BIC)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    K_BIC <- densities %>% length
    ts <- getTs(working_dat, densities)
    zs <- getZs(ts)
    plot(working_dat, col=zs, pch=19)
    ent <- getEntropy(ts)
    
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
runAnalysis(d2,
            model_names="VVI",
            output_dir="Figures/d2/")
