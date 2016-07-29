source("FinalStats.R")


statDistanceMatrix = function(tip) {
  nstat <-2
  outputdir <- "./trees"
  inputdir <- "./trees"
  inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
  trees <- read.tree(inputfile, keep.multi = TRUE)
 
  n <- length(trees)
  statMatrix <- matrix(0, nrow = n, ncol= nstat) 
  for (j in (1:n)) {
    root <- trees[[j]]
    tree <- setTree(root, tip)
    parents <- setParents(tree, tip)
    Nis <- setNis (root)
    MHat <- setMaxDist4node2tips(parents, tip)
    NTips <- extractNTips(tree, tip)
    cherries <- extractCherries(tree, tip)
    pitchforks <- extractPitchforks(tree, tip)
    Width <- extractWidth(Nis,MHat, tip)
    
    # # calculate Ic
  statMatrix[j,1] <- statIc(root)
    # 
    # # calculate  Sackin
  statMatrix[j,2] <- statSackin(root)
    # 
    
  }
  DM = matrix(0, nrow= n, ncol = n)
  DM = as.matrix(dist(statMatrix, method= "euclidean", diag = FALSE, upper = FALSE, p = 2))
  Ds <- DM^2
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





nstat <- 2
met = 'NNI'
main = function(tip,met,nstat){
  ntip=tip
  inputdir <- "./trees"
  inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
    trees <- read.tree(inputfile, keep.multi = TRUE)
    # fill the Stat Matrix 
    n <- length(trees)
    statMatrix <- matrix(0, nrow = n, ncol= nstat) 
    for (j in (1:n)) {
      root <- trees[[j]]
      Nis <- setNis (root)
      statMatrix[j,1] <- statIc(root)
      statMatrix[j,2] <- statSackin(root)
      
    }
    }
Colless=statMatrix[,1]
Sakin=statMatrix[,2]
plot(Colless,Sakin)
plot(trees[[which(Sakin==max(Sakin))]])
plot(trees[[which(SaKin==max(Sakin))]])

