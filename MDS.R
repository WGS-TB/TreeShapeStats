source("FinalStats.R")


MDS = function(tip,met) {
  nstat <-2
  outputdir <- "./trees"
  inputdir <- "./trees"
  inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
  trees <- read.tree(inputfile, keep.multi = TRUE)
  # fill the Stat Matrix 
  n <- length(trees)
  metrics <- DistanceMatrix(trees, tip, metric = met)
  Ds <- metrics^2
  H=diag(rep(1,n))-1/n
  XD <- -(1/2) * H %*% Ds %*% H
  eye = eigen(XD , symmetric = TRUE, EISPACK = FALSE)
  values <- eye$values
  vectors <- eye$vectors
  treeLocation = matrix(0, nrow = n , ncol = 2)
  for (i in (1:n)){
    treeLocation[i, 1] <- (values[1]^(1/2)) * vectors[i, 1]
    treeLocation[i, 2] <- (values[2]^(1/2)) * vectors[i, 2]
  }
  
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  hist(1)
  hist(2)
  hist(3)
  
  title = paste("Distibution of Trees with ", tip , " Leaves")
  plot(treeLocation, main = list (title, col= "black" ), xlab = "X", ylab = "Y")
  idx = which(treeLocation[,1]==max(treeLocation[,1]))
  plot (trees[[idx[1]]], type="c", main = list("Outlier with maximum X value" , col = "black") )
  idx = which(treeLocation[,1]==min(treeLocation[,1]))
  plot (trees[[idx[1]]],  type="c",main = list("Outlier with minimum X value" , col = "black") )
}  


