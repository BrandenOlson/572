source("cluster.R")

shift <- 0.1
unif_dat <- rbind(cbind(runif(200, 0, 1),
                             runif(200, 0, 1)
                             ),
                  rmvnorm(300, mean=c(0.5,1.7), sigma=0.1*diag(2))
                       )

runAnalysis(unif_dat[, 1:2],
            output_dir="Figures/unif2/",
            contour_levels=1/10
            )
