source("cluster.R")

# Analysis using initialization parameters from Baudry et al.
# Only run on the dataset in Example 4.3

getInitParameters <- function(param_dir,
                              K,
                              d=2
                              ) {
    bic_p <- fread(paste0(param_dir, "BIC_p"))
    
    bic_mu <- fread(paste0(param_dir, "BIC_mu"))
    
    K <- 5
    bic_Sigma <- array(NA, c(d, d, K))
    for(k in 1:K) {
        bic_Sigma[, , k] <- fread(paste0(param_dir, "BIC_Sigma_", k)) %>% 
            as.matrix
    }
    
    params <- list(Prop=bic_p %>% as.numeric,
                   Mean=bic_mu %>% t %>% as.matrix,
                   Variance=bic_Sigma)
                   
    return(params)
}

init_params <- getInitParameters("../Baudry/",
                                 K=5)
runAnalysis(d3,
            model_names=c("VII"),
            output_dir="Figures/d3_B/",
            init_params=params
            )
