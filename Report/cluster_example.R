library(mvtnorm)

xs <- rbind(
            rmvnorm(200, mean=c(0, 0), sigma=diag(0.5, 2)),
            rmvnorm(300, mean=c(3, 3), sigma=matrix(c(0.5, 0.3, 0.3, 0.5), 
                                                    nrow=2, ncol=2)),
            rmvnorm(100, mean=c(4, -1), sigma=matrix(c(0.6, -0.4, -0.4, 0.6),
                                                     nrow=2, ncol=2))
            )

labels <- kmeans(xs, 3)$cluster

pdf("unclustered.pdf")
plot(xs, asp=1, pch=19, xlab="x", ylab="y", cex.lab=1.5, cex.axis=1.5)
dev.off()

pdf("clustered.pdf")
plot(xs, asp=1, pch=19, col=labels + 1, xlab="x", ylab="y", cex.lab=1.5, cex.axis=1.5)
dev.off()
