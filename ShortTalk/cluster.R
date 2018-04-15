library(MASS)

sample_a <- mvrnorm(100, mu=c(0, 0), Sigma=0.2*matrix(c(1, 0.5, 0.5, 1), 2, 2))
sample_b <- mvrnorm(100, mu=c(3, 0), Sigma=matrix(c(1, 0.5, 0.5, 3), 2, 2))
sample_c <- mvrnorm(100, mu=c(1.5, -1.75), Sigma=0.5*matrix(c(2, 0.1, 0.1, 1), 2, 2))
true_labels <- c(rep(1, 100), rep(2, 100), rep(3, 100))

full_sample <- rbind(sample_a, sample_b, sample_c)

kmeans_object <- kmeans(full_sample, 3)
labels <- kmeans_object$cluster

pdf("clustered.pdf", width=10, height=6)
par(mfrow=c(1, 2))
plot(full_sample, pch=19, xlab="x", ylab="y", main="Observed data")
plot(full_sample, col=labels + 1, pch=19, xlab="x", ylab="y",
     main="Clustered data")
dev.off()
