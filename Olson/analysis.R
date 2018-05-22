library(strucchange)
library(xtable)

source("cluster.R")
source("EM.R")
source("processData.R")

printTable <- function(dat, filename, digits=3, hline_pos=0) {
    sink(filename)
    cat("\\begin{center}", '\n')
    dat %>%
        xtable::xtable(digits=digits) %>%
        print(
              floating=FALSE,
              include.rownames=FALSE,
              latex.environments="center",
              hline.after=c(0, hline_pos),
              sanitize.text.function=function(x){x}
             )
    cat("\\end{center}", '\n')
    sink()
}

runAnalysis <- function(working_dat,
                        model_names=NULL,
                        output_dir,
                        K_min=1,
                        K_max=12,
                        contour_levels=0.01,
                        plot_types="2d"
                       ) {
    # Plot unlabeled data
    point_size <- ifelse(plot_types == "4d", 0.1, 0.5)
    pdf(paste0(output_dir, "unclustered.pdf"), width=6, height=6)
    plot(working_dat, pch=19, asp=1, cex=point_size, cex.axis=2, cex.lab=2)
    dev.off()

    mc <- getModel(working_dat,
                   G=K_min:K_max,
                   model_names)
    mc_params <- mc %>% getModelParameters
    densities <- mc_params %>% getDensities
    K_BIC <- densities %>% length
    ts <- getTs(dat=working_dat, components=densities, K=K_BIC)
    zs <- getZs(ts, densities)
    entropy_bic <- ts %>% getEntropy
    
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
    ts_icl <- getTs(dat=working_dat,
                    components=icl_densities,
                    K=K_ICL)

    zs_icl <- ts_icl %>%
        getZs(components=icl_densities)
    entropy_icl <- ts_icl %>% getEntropy

    K_ent_df <- data.table(
                           Method=c("BIC", "ICL"),
                           K=c(K_BIC, K_ICL),
                           Entropy=c(entropy_bic, entropy_icl)
                          )


    ICL_prefix <- paste0(output_dir, "ICL", K_ICL)

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
    k_ent_dfs <- Map(function(x, k) {
                         plot_name <- paste0(output_dir, "merged_", k)
                         ts <- getTs(dat=working_dat,
                                     components=x,
                                     K=K)
                         zs <- ts %>%
                             getZs(components=x)
                         entropy_k <- ts %>% getEntropy
                         k_ent_df <- data.table(Method="Merged",
                                                K=k,
                                                Entropy=entropy_k)

                         if(plot_types=="2d") {
                             plotDensities(x,
                                       output_prefix=plot_name,
                                       dat=working_dat,
                                       zs=zs,
                                       contour_levels=contour_levels
                                       )
                         } else if(plot_types=="3d") {
                              plotDensities3D(dat=working_dat,
                                              prefix=plot_name,
                                              zs=zs)
                         } else if(plot_types=="4d") {
                              plotCD3Clusters(dat=working_dat,
                                              components=x,
                                              zs=zs,
                                              prefix=plot_name)
                         }

                         return(k_ent_df)
                     },
                     cluster_seq,
                     K_BIC:1
                 ) %>% Reduce(rbind, .)

    K_ent_df <- rbind(K_ent_df,
                      k_ent_dfs)

    print(dim(K_ent_df))
    printTable(K_ent_df,
               paste0(output_dir, "K_ent.tex")
               )

    pdf(paste0(output_dir, "total_entropy.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$K,
                           xlab="K",
                           ylab="Total entropy",
                           do_regression=TRUE
                          )
    dev.off()
    pdf(paste0(output_dir, "diff_entropy.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$DeltaEntropy,
                           entropy$K - 1,
                           ylab="Difference in entropy from K + 1",
                           xlab="K")
    dev.off()
    pdf(paste0(output_dir, "cum_count.pdf"), width=10, height=6)
    plotPiecewiseFitToData(entropy$TotalEntropy,
                           entropy$NumMergedCumSum,
                           xlab="Cumulative count of merged observations",
                           ylab="Total entropy",
                           do_regression=TRUE
                           )
    dev.off()
    pdf(paste0(output_dir, "norm_diff.pdf"), width=10, height=6)
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
                                   do_regression=FALSE
                                   ) {
    valid_indices <- which(!is.na(y_col))
    ys <- y_col[valid_indices]
    xs <- x_col[valid_indices]
    n <- length(xs)
    yrange <- c(0, max(ys))
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    plot(ys ~ xs, pch=19, xlab=xlab, ylab=ylab, ylim=yrange, 
         cex=1.5, cex.axis=2, cex.lab=2)
    a1 <- b1 <- a2 <- b2 <- NA
    lse <- NA
    if(n > 2 && do_regression) {
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
    
}
                                  
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
stop()
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
runAnalysis(d42,
            output_dir="Figures/d42/",
            plot_types="3d"
            )

# Need to make a plot just for the cluster labels,
# since it's difficult to plot 3D contours

stop()
