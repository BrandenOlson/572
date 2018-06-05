source("cluster.R")
source("EM.R")
source("processData.R")


runAnalysis(d51_raw,
            output_dir="Figures/d51/",
            plot_types="4d",
            K_max=16,
            model_names="VVV"
            )
stop()
runAnalysis(d52_raw,
            output_dir="Figures/d52/",
            K_max=15,
            model_names="VVV",
            plot_types="4d"
            )
runAnalysis(unif_dat[, 1:2],
            model_names="VII",
            output_dir="Figures/unif/",
            contour_levels=1/10)
d_tmp <- rbind(
               x=cbind(runif(200, 0, 2),
                     runif(200, 1.1, 2)),
               y=cbind(runif(200, 0, 2),
                     runif(200, 0, 0.9)),
               rmvn(400, mu=c(1, 3), sigma=diag(1/5, 2))
               )

runAnalysis(d_tmp,
            model_names="VII",
            output_dir="Figures/square_circle/",
            contour_levels=0.1
            )
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
runAnalysis(d42,
            output_dir="Figures/d42/",
            plot_types="3d"
            )
