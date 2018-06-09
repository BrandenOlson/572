source("cluster.R")
source("EM.R")
source("processData.R")


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
runAnalysis(d51_raw,
            output_dir="Figures/d51/",
            plot_types="4d",
            K_max=16,
            model_names="VVV"
            )
runAnalysis(d52_raw,
            output_dir="Figures/d52/",
            K_max=15,
            model_names="VVV",
            plot_types="4d"
            )
