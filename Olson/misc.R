
# Proof that the paper and mclust give different clusterings
mc <- mclustBIC(d1, G=6, modelNames="VVV") %>% Mclust(data=d1, x=.)
