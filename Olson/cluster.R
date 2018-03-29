library(data.table)
library(dplyr)

readData <- function(filename, change_names=TRUE) {
    dat <- fread(filename)
    if(change_names) {
        dat_names <- c("x", "y")
        if(ncol(dat) == 3) {
            dat_names <- c(dat_names, "z")
        } 
        names(dat) <- dat_names
    }
    return(dat)
}

d1 <- readData("Data/4.1.csv")
d2 <- readData("Data/4.2.csv")
d3 <- readData("Data/4.3.csv")
d41 <- readData("Data/4.4.1.csv")
d42 <- readData("Data/4.4.2.csv")

d51_raw <- readData("Data/GvHD+.csv", change_names=FALSE)
names(d51_raw) <- c("CD4", "CD8beta", "CD3", "CD8")
d51 <- d51_raw[d51_raw$CD3 > 280, ]



d52_raw <- readData("Data/GvHD-.csv", change_names=FALSE)
names(d52_raw) <- c("CD4", "CD8beta", "CD3", "CD8")
d52 <- d52_raw[d52_raw$CD3 > 280, ]
