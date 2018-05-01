d1 <- readData("Data/4.1.csv")
d2 <- readData("Data/4.2.csv")
d3 <- readData("Data/4.3.csv")
d41 <- readData("Data/4.4.1.csv")
d42 <- readData("Data/4.4.2.csv")

d_square_circle <- readData("Data/square_circle.csv")

d51_raw <- readData("Data/GvHD+.csv", change_names=FALSE)
names(d51_raw) <- c("CD4", "CD8beta", "CD3", "CD8")

d52_raw <- readData("Data/GvHD-.csv", change_names=FALSE)
names(d52_raw) <- c("CD4", "CD8beta", "CD3", "CD8")

