require(combinat)
require(rARPACK)
require(bigmemory)
library(igraph)
library(apTreeshape)
library(phangorn)
library(ape)
library(hash)


# returns a matrix(ntip-1,3), each row includes (#internal node, #child1, #child2) 
setTree = function (root , ntip){
  edges <- root$edge
  tree <- matrix(0, nrow= (ntip-1), ncol =3)
  for (i in ((ntip+1):(2*ntip -1))) {
    tree[i-ntip,1] <- i
    indx <- which(edges[,1] %in% i)
    tree[i-ntip,2] = edges[indx[1],2]
    tree[i-ntip,3] = edges[indx[2],2]
  }
  return (tree)
}

# returns a vector that represents parent of each node in tree. 
setParents = function (tree, ntip){
  #tree <- setTree(root, ntip)
  parent <- numeric(2*ntip-1)
  for (i in 1:(2*ntip-1)){
    idx <- which(tree == i)
    if (i <= ntip +1 ){
      if (i == ntip +1){
        #  i is root
        parent[i] <- -1 
      }else {
        # i is tip
        if (idx[1]%%(ntip-1) == 0) {row <- (ntip-1)}
        else { row = idx[1]%%(ntip-1) }
        parent[i] <- tree[row, 1]
      }
    }else {
      # i is an internal node, therefore it exists 2 times in table
      # the first value is less than ntip
      if (idx[2]%%(ntip-1) == 0) {row <- (ntip-1)}
      else { row = idx[2]%%(ntip-1) }
      parent[i] <- tree[row, 1]
    }
  }
  return (parent)
}

# computes the number of nodes between external node to root (root is included)

setNis = function(root) {
  Nis<-node.depth.edgelength(root)
  return (Nis)
}


# returns the number of left node and right node of each internal node in a tree rooted by root which has ntip tips 
statI2 = function (tree, NTips, ntip, n){
  #tree is already defined as global variable
  #tree = setTree(root, ntip)
  # calculate I2
  I2 <- 0
  for (i in ((ntip+1):(2*ntip-1))){
    ch1 <- tree[i-ntip,2]
    if (ch1 > ntip){
      ri <- NTips[ch1]
    } else {
      ri <- 1
    }
    ch2 <- tree[i-ntip,3]
    if (ch2 > ntip){
      si <- NTips[ch2]
    } else {
      si <- 1
    }
    if ((ri + si) > 2) {
      I2 <- I2 + (abs(ri -si)/abs (ri +si -2)) 
    }        
  }
  I2 <- I2 / (n -2)
  return (I2)
}

# Helper function that returns number of tips for all nodes in tree
# This value for tips is zero 
extractNTips = function (tree, ntip) {
  NTips <<-rep(1, 2*ntip -1)
  for (i in ((ntip+1): (2*ntip -1))){
    NTips[i] <- findNtip(tree, i, ntip, NTips)
  }
  return (NTips)
}
findNtip = function (tree, node, ntip, Ntip) {
  if (node <= ntip ) {
    return (1)
  }
  return (findNtip(tree, tree[node-ntip, 2], ntip, Ntip) + findNtip(tree, tree[node-ntip, 3], ntip, Ntip))
}

# calculates cobination of Sackin and Colless stat
statSaless = function (root){
  root=as.treeshape(root)
  Ns=sackin(root, norm = NULL)
  Ic=colless(root)
  
  return (1.300745*Ns+Ic)
}


# calculates Sackin stat
statSackin = function (root){
  root=as.treeshape(root)
  Ns=sackin(root, norm = NULL)
  return (Ns)
}


# calculate Variance stat
statVarianc = function (Nis,ntip){
  
  Nvar=var(Nis[1:ntip]) 
  return (Nvar)
}

# calculates Ic stat
statIc = function (root){
  root=as.treeshape(root)
  Ic=colless(root)
  return (Ic)
}

# calculates B2
statB2 = function(Nis, ntip){
  B2 <- 0 
  for (i in 1:ntip){
    B2 <- B2 + (Nis[i] / (2^ Nis[i]))
  }
  return (B2)
}

# calculates B1 
statB1 = function(MHat, ntip){
  #MHat <- setMaxDist4node2tips(root, ntip)
  B1 <- 0
  # number of internal nodes excluding root
  for (i in ((ntip+2):(2*ntip-1))){
    B1 <- B1 + (1/MHat[i])
  }
  return (B1)
}

# helper function 
setMaxDist4node2tips = function (parents, ntip){
  # node is an internal node excluding root
  Mdist <- numeric(2*ntip-1)
  for (i in ((2*ntip -1):(ntip+1))) {
    list <- which(parents==i)
    for (j in (1:2)){
      if ( Mdist[i] <= Mdist[list[j]]) {
        Mdist[i] <- Mdist[list[j]] + 1 
      }
    }
  }
  return (Mdist)
}

# a helper function that returns number of cheries for a given internal node (which is an integer number in [ntip+1 .. 2*ntip -1]) in tree
# for nodes that are tips it returns zero.   

extractCherries = function (tree, ntip) {
  cherrys <<-rep(0, 2*ntip -1)
  for (i in ((ntip+1): (2*ntip -1))){
    cherrys[i] <- findCherry(tree, i, ntip, cherrys)
  }
  return (cherrys)
}
findCherry = function (tree, node, ntip, cherry) {
  if (node <= ntip ) {
    return (0)
  }
  idx <- which(tree[,1] %in% node)
  if (tree[idx,2] <= ntip & tree[idx,3] <= ntip) {
    return (1)
  }   
  return (findCherry(tree, tree[idx, 2], ntip, cherry) + findCherry(tree, tree[idx, 3], ntip, cherry))
}

# a helper function that returns number of pitchforks for a given internal node (which is an integer number in [ntip+1 .. 2*ntip -1]) in tree
# for nodes that are tips it returns zero.   

extractPitchforks = function (tree, ntip) {
  pitchforks <<-rep(0, 2*ntip -1)
  for (i in ((ntip+1): (2*ntip -1))){
    pitchforks[i] <- findPitchforks(tree, i, ntip, pitchforks, NTips)
  }
  return (pitchforks)
}
findPitchforks = function (tree, node, ntip, pitchforks, NTips) {
  if (node <= ntip ) {
    return (0)
  }
  if (NTips[node] ==3){
    return (1) 
  }
  idx <- which(tree[,1] %in% node)
  return (findPitchforks(tree, tree[idx, 2], ntip, pitchforks, NTips) + findPitchforks(tree, tree[idx, 3], ntip, pitchforks, NTips))
}


# extracts the width of each level in tree. Width of level zero is 1 which includes root of tree.
extractWidth = function (Nis, MHat, ntip){
  Mroot <- MHat[ntip+1] 
  Nwidth <- numeric(Mroot)
  for (i in (1:Mroot)){
    list <- which(Nis==i)
    Nwidth[i] <- length(list)
  }
  return (Nwidth)
}



statTP = function (tree, NTips, pitchforks, ntip){
  sum <- 0
  for (i in ((ntip+1):(2*ntip -1))){ 
    ch1 <- tree[i-ntip,2]
    if (ch1 > ntip){
      Tri <- NTips[ch1]
      Pri <- pitchforks[ch1]
    }else{
      Tri <- 1
      Pri <- 0
    }  
    ch2 <- tree[i-ntip,3]
    if (ch2 > ntip){
      Tsi <- NTips[ch2]
      Psi <- pitchforks[ch2]
    } else {
      Tsi <- 1
      Psi <- 0
    }  
    sum = sum + (Tri - Tsi)^2+(Pri - Psi)^2
    
  } 
  return (sum)
}
getCanonical = function (root, ntip){
  tree1 <- setTree(root, ntip)
  NTip <- extractNTips(tree1, ntip)
  for (i in (1:(ntip-1))){
    nodel <- tree1[i,2]
    noder <- tree1[i,3]
    if (NTip[nodel] < NTip[noder]){
      tree1[i,2] <- noder
      tree1[i,3] <- nodel
    }
  } 
  return (tree1)
}

# returns a matrix with 4 columns and ntip-1 rows. 
# It contains tree in 3 first column and the last columns is 1 if there is rotation, and 0 otherwise. 
findEquals = function (root, ntip){
  list = array(0,dim = c(ntip-1))
  tree1 <- setTree(root, ntip)
  NTip <- extractNTips(tree1, ntip)
  for (i in (1:(ntip-1))){
    nodel <- tree1[i,2]
    noder <- tree1[i,3]  
    # Save different shapes of tree by turning tree around a node
    if ((NTip[nodel] == NTip[noder]) & (NTip[nodel] > 1)){
      list[i] <- 1
    }
  } 
  Ctree <- matrix(c(tree1, list), nrow=ntip-1 , ncol= 4)
  return (Ctree)
}

getNNI = function (root, ntip) {
  proned <<- list()
  neighbours <- nni(root)
  for (i in (1: length(neighbours))){
    cantree <- getCanonical(neighbours[[i]], ntip)
    strNew <- extractNewick(cantree, ntip+1, ntip)
    flag <- TRUE
    for (t in proned) {
      if (strNew == t) {
        flag <- FALSE
        break
      }
    }
    if (flag) {
      .GlobalEnv$proned[[(length(.GlobalEnv$proned)+1)]] <- strNew
    }   
  }
  return (proned)
}

extractNewick = function(tree, node, ntip){
  if (node <= ntip){
    return ("")
  }
  i <- which(tree[,1]==node)
  return (paste ("(", extractNewick(tree, tree[i,2], ntip), "," ,extractNewick(tree, tree[i,3], ntip), ")", sep="" ))
}

# extracts different sahpe of a tree based on rotation on 4th column of Ctree, and it returns as a list of trees
extractEquals = function (Ctree, ntip) {
  ls <- list()
  ls[[length(ls) +1]] <- Ctree[,1:3]
  for (j in (1: (ntip-1))) {
    if (Ctree[j,4] == 1) {
      hls = list()
      for ( p in ls ){
        hls [[length(hls) +1]] <- p
        t <- p[j,2] 
        p[j,2] <- p[j,3]
        p[j,3] <- t 
        hls[[length(hls) +1]] <- p
      }
      ls <- hls
    }
  }
  return (ls)
} 

# returns a distance matrix (n*n) for all the trees in the input file by using the NNI metric
DistanceMatrix = function (trees, ntip, metric) {
  n <- length(trees)
  if (metric=="NNI"){
    # creates a hash table for (newickformat of each phylogenic tree with given leaves, an index in [1..n]) 
    ha <- new.env()
    for (i in (1:n)){
      tree1 <- getCanonical(trees[[i]],ntip)
      Ctree <- findEquals(trees[[i]], ntip)
      listTree <- extractEquals(Ctree, ntip)
      ll <-length(listTree)
      listTree[[ll+1]] <- tree1
      for (k in listTree) {
        key <- extractNewick (k, ntip+1, ntip)
        ha[[key]] <- i
      }    
    }
    # creates the adjacency matrix
    adj <- matrix(0, nrow= n, ncol =n)
    for (i in (1:n)){ 
      row <- gsub("\\;","", gsub("(\\d*)(\\:)(\\d+)", "", write.tree(trees[[i]])))
      proned <- getNNI(trees[[i]], ntip)
      for (j in (1:length(proned))){
        col <- proned[[j]]
        if (col != row) {
          
          adj[ha[[row]], ha[[col]]] <- 1
        }
      }
      
    }
    ig <- graph.adjacency(adj, mode="undirected")
    distMatrix <- shortest.paths(ig, v=V(ig), to=V(ig))
    
    return (distMatrix)
  } 
  if (metric=="SPR"){
    spr = read.table(paste(inputdir, "/SPR/spr", ntip, ".txt", sep=""))
    return (as.matrix(spr))
  }
}  

# There are 8 stat considered to be evaluated by NNI metric on phylognic trees of #tip leaves. These stats are stored in statMatrix:
# 1- Ic, 2- Sackin, 3- Variance, 4- I2,  5- B1, 6- B2
# 7- TP, 8- Saless  
main = function(low,high,met,nstat){
  
  outputdir <- "C:\\Users\\Maryam\\Desktop\\Courses\\Spring 2016\\CMPT829\\Project\\trees"
  filename <- paste(outputdir,"/resolutionspr.csv",sep='')
  cat(", Ic, Sackin, Variance, I2, B1, B2,TP, Saless, TCP\n", file = filename, append = TRUE)
  inputdir <- "C:\\Users\\Maryam\\Desktop\\Courses\\Spring 2016\\CMPT829\\Project\\trees"
  resolution <- matrix (0, nrow = high -low +1 , ncol= nstat)
  for (tip in (low:high)) {
    # reading all phylogenic trees for a given number of tips from a file
    inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
    trees <- read.tree(inputfile, keep.multi = TRUE)
    # fill the Stat Matrix 
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
      # calculate Ic
      statMatrix[j,1] <- statIc(root)
      
      # calculate  Sackin
      statMatrix[j,2] <- statSackin(root)
      
      # calculate stat Variance
      statMatrix[j,3] <- statVarianc(Nis, tip)
      
      # calculate stat I2
      statMatrix[j,4] <- statI2(tree, NTips, tip, n)
      
      # calculate stat B1
      statMatrix[j,5] <- statB1(MHat, tip)
      
      # calculate stat B2
      statMatrix[j,6] <- statB2(Nis, tip)
      
      # calculate stat Saless: combination of sackin and colless on each node for all nodes in tree
      statMatrix[j,7] <- statSaless(root)
      
      # calculate stat TCP: combination of number of Tips, Cherries and Pitchforks for all internal nodes. 
      statMatrix[j,8] <- statTP(tree, NTips,  pitchforks, tip)
      
      
    }
    
    
    # Calcolates vector XF for all stats. each column corresponds to one stat 
    H <- diag(rep(1,n))-1/n
    
    XF <- matrix(0, nrow = n, ncol= nstat)
    for (s in (1:nstat)){
      YF <- statMatrix[,s]
      HYF <- H %*% YF
      XF[,s] <- HYF / norm(HYF, type="2") 
    }
    
    # Calculate the NNI resolution for all stats
    metrics <- DistanceMatrix(trees, tip, metric = met)
    Ds <- metrics^2
    tmp <- -(1/2) * H %*% Ds %*% H
    eye = eigen(tmp , symmetric = TRUE ,only.values = TRUE, EISPACK = FALSE)
    values <- eye$values
    dist <- values[1] - values[n]
    cat(tip, file = filename, append=TRUE)
    for (s in (1:nstat)){
      result <- (-0.5) * ((t(XF[,s]) %*% Ds) %*% XF[,s])
      # scale the result
      resolution[tip-(low-1),s] <- (result - values[n])/ dist 
      cat(", ", file = filename, append=TRUE)
      cat(resolution[tip-(low-1),s], file = filename, append=TRUE)
    }
    cat("\n" , file = filename, append=TRUE)
  }  
}

statDistanceMatrix = function(tip) {
  nstat <- 2
  outputdir <- "C:\\Users\\Maryam\\Desktop\\Courses\\Spring 2016\\CMPT829\\Project\\trees"
  inputdir <- "C:\\Users\\Maryam\\Desktop\\Courses\\Spring 2016\\CMPT829\\Project\\trees"
  inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
  trees <- read.tree(inputfile, keep.multi = TRUE)
  # fill the Stat Matrix 
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
    
    # calculate Ic
    statMatrix[j,1] <- statIc(root)
    
    # calculate  Sackin
    statMatrix[j,2] <- statSackin(root)
    
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
  
  layout(matrix(c(1,1,1,2,3,4),2,3, byrow = TRUE))
  hist(1)
  hist(2)
  hist(3)
  hist(4)
  title = paste("Distibution of Trees with ", tip , " Leaves")
  plot(treeLocation, main = list (title, col= "blue" ), xlab = "X", ylab = "Y")
  idx = which(treeLocation[,2]==min(treeLocation[,2]))
  plot (trees[[idx]], main = list("Outlier with minimum Y value" , col = "red") )
  idx = which(treeLocation[,1]==max(treeLocation[,1]))
  plot (trees[[idx]], main = list("Outlier with maximum X value" , col = "red") )
  idx = which(treeLocation[,1]==min(treeLocation[,1]))
  plot (trees[[idx]], main = list("Outlier with minimum X value" , col = "red") )
}  





