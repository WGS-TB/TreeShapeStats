

require(combinat)
require(rARPACK)
library(igraph)
library(apTreeshape)
library(phangorn)
library(ape)
library(hash)
library(Matrix)

# returns a matrix(tip-1,3), each row includes (#internal node, #child1, #child2) 
setTree = function (root , tip){
  edges <- root$edge
  tree <- matrix(0, nrow= (tip-1), ncol =3)
  for (i in ((tip+1):(2*tip -1))) {
    tree[i-tip,1] <- i
    indx <- which(edges[,1] %in% i)
    tree[i-tip,2] = edges[indx[1],2]
    tree[i-tip,3] = edges[indx[2],2]
  }
  return (tree)
}

# returns a vector that represents parent of each node in tree. 
setParents = function (tree, tip){
  #tree <- setTree(root, tip)
  parent <- numeric(2*tip-1)
  for (i in 1:(2*tip-1)){
    idx <- which(tree == i)
    if (i <= tip +1 ){
      if (i == tip +1){
        #  i is root
        parent[i] <- -1 
      }else {
        # i is tip
        if (idx[1]%%(tip-1) == 0) {row <- (tip-1)}
        else { row = idx[1]%%(tip-1) }
        parent[i] <- tree[row, 1]
      }
    }else {
      # i is an internal node, therefore it exists 2 times in table
      # the first value is less than tip
      if (idx[2]%%(tip-1) == 0) {row <- (tip-1)}
      else { row = idx[2]%%(tip-1) }
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


# returns the number of left node and right node of each internal node in a tree rooted by root which has tip tips 
statI2 = function (tree, NTips, tip, n){
  #tree is already defined as global variable
  #tree = setTree(root, tip)
  # calculate I2
  I2 <- 0
  for (i in ((tip+1):(2*tip-1))){
    ch1 <- tree[i-tip,2]
    if (ch1 > tip){
      ri <- NTips[ch1]
    } else {
      ri <- 1
    }
    ch2 <- tree[i-tip,3]
    if (ch2 > tip){
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
extractNTips = function (tree, tip) {
  NTips <<-rep(1, 2*tip -1)
  for (i in ((tip+1): (2*tip -1))){
    NTips[i] <- findNtip(tree, i, tip, NTips)
  }
  return (NTips)
}
findNtip = function (tree, node, tip, Ntip) {
  if (node <= tip ) {
    return (1)
  }
  return (findNtip(tree, tree[node-tip, 2], tip, Ntip) + findNtip(tree, tree[node-tip, 3], tip, Ntip))
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
statVarianc = function (Nis,tip){
  
  Nvar=var(Nis[1:tip]) 
  return (Nvar)
}

# calculates Ic stat
statIc = function (root){
  root=as.treeshape(root)
  Ic=colless(root)
  return (Ic)
}

# calculates B2
statB2 = function(Nis, tip){
  B2 <- 0 
  for (i in 1:tip){
    B2 <- B2 + (Nis[i] / (2^ Nis[i]))
  }
  return (B2)
}

# calculates B1 
statB1 = function(MHat, tip){
  #MHat <- setMaxDist4node2tips(root, tip)
  B1 <- 0
  # number of internal nodes excluding root
  for (i in ((tip+2):(2*tip-1))){
    B1 <- B1 + (1/MHat[i])
  }
  return (B1)
}

# helper function 
setMaxDist4node2tips = function (parents, tip){
  # node is an internal node excluding root
  Mdist <- numeric(2*tip-1)
  for (i in ((2*tip -1):(tip+1))) {
    list <- which(parents==i)
    for (j in (1:2)){
      if ( Mdist[i] <= Mdist[list[j]]) {
        Mdist[i] <- Mdist[list[j]] + 1 
      }
    }
  }
  return (Mdist)
}

# a helper function that returns number of cheries for a given internal node (which is an integer number in [tip+1 .. 2*tip -1]) in tree
# for nodes that are tips it returns zero.   

extractCherries = function (tree, tip) {
  cherrys <<-rep(0, 2*tip -1)
  for (i in ((tip+1): (2*tip -1))){
    cherrys[i] <- findCherry(tree, i, tip, cherrys)
  }
  return (cherrys)
}
findCherry = function (tree, node, tip, cherry) {
  if (node <= tip ) {
    return (0)
  }
  idx <- which(tree[,1] %in% node)
  if (tree[idx,2] <= tip & tree[idx,3] <= tip) {
    return (1)
  }   
  return (findCherry(tree, tree[idx, 2], tip, cherry) + findCherry(tree, tree[idx, 3], tip, cherry))
}

# a helper function that returns number of pitchforks for a given internal node (which is an integer number in [tip+1 .. 2*tip -1]) in tree
# for nodes that are tips it returns zero.   

extractPitchforks = function (tree, tip) {
  pitchforks <<-rep(0, 2*tip -1)
  for (i in ((tip+1): (2*tip -1))){
    pitchforks[i] <- findPitchforks(tree, i, tip, pitchforks, NTips)
  }
  return (pitchforks)
}
findPitchforks = function (tree, node, tip, pitchforks, NTips) {
  if (node <= tip ) {
    return (0)
  }
  if (NTips[node] ==3){
    return (1) 
  }
  idx <- which(tree[,1] %in% node)
  return (findPitchforks(tree, tree[idx, 2], tip, pitchforks, NTips) + findPitchforks(tree, tree[idx, 3], tip, pitchforks, NTips))
}


# extracts the width of each level in tree. Width of level zero is 1 which includes root of tree.
extractWidth = function (Nis, MHat, tip){
  Mroot <- MHat[tip+1] 
  Nwidth <- numeric(Mroot)
  for (i in (1:Mroot)){
    list <- which(Nis==i)
    Nwidth[i] <- length(list)
  }
  return (Nwidth)
}



statTP = function (tree, NTips, pitchforks, tip){
  sum <- 0
  for (i in ((tip+1):(2*tip -1))){ 
    ch1 <- tree[i-tip,2]
    if (ch1 > tip){
      Tri <- NTips[ch1]
      Pri <- pitchforks[ch1]
    }else{
      Tri <- 1
      Pri <- 0
    }  
    ch2 <- tree[i-tip,3]
    if (ch2 > tip){
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
getCanonical = function (root, tip){
  tree1 <- setTree(root, tip)
  NTip <- extractNTips(tree1, tip)
  for (i in (1:(tip-1))){
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
findEquals = function (root, tip){
  list = array(0,dim = c(tip-1))
  tree1 <- setTree(root, tip)
  NTip <- extractNTips(tree1, tip)
  for (i in (1:(tip-1))){
    nodel <- tree1[i,2]
    noder <- tree1[i,3]  
    # Save different shapes of tree by turning tree around a node
    if ((NTip[nodel] == NTip[noder]) & (NTip[nodel] > 1)){
      list[i] <- 1
    }
  } 
  Ctree <- matrix(c(tree1, list), nrow=tip-1 , ncol= 4)
  return (Ctree)
}
#return all trees wich are equal to a canonical tree
findEqualsCanonical = function (tree, tip){
  list = array(0,dim = c(tip-1))
  #tree1 <- setTree(root, tip)
  NTip <- extractNTips(tree, tip)
  for (i in (1:(tip-1))){
    nodel <- tree[i,2]
    noder <- tree[i,3]  
    # Save different shapes of tree by turning tree around a node
    if ((NTip[nodel] == NTip[noder]) & (NTip[nodel] > 1)){
      list[i] <- 1
    }
  } 
  Ctree <- matrix(c(tree, list), nrow=tip-1 , ncol= 4)
  return (Ctree)
}

getNNI = function (root, tip) {
  proned <<- list()
  neighbours <- nni(root)
  for (i in (1: length(neighbours))){
    cantree <- getCanonical(neighbours[[i]], tip)
    strNew <- extractNewick(cantree, tip+1, tip)
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

extractNewick = function(tree, node, tip){
  if (node <= tip){
    return ("")
  }
  i <- which(tree[,1]==node)
  return (paste ("(", extractNewick(tree, tree[i,2], tip), "," ,extractNewick(tree, tree[i,3], tip), ")", sep="" ))
}

# extracts different sahpe of a tree based on rotation on 4th column of Ctree, and it returns as a list of trees
extractEquals = function (Ctree, tip) {
  ls <- list()
  ls[[length(ls) +1]] <- Ctree[,1:3]
  for (j in (1: (tip-1))) {
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


main = function(low,high,nstat){
  outputdir <- "/Users/maryam/Google Drive/Research/Tree Stat/Codes/New"
  filename <- paste(outputdir,"/resolution.csv",sep='')
  cat(", Ic, Sackin, Variance, I2, B1, B2,TP, Saless\n", file = filename, append = TRUE)
  inputdir <- "/Users/maryam/Google Drive/Research/Tree Stat/Codes/Trees"
  resolutionNew <- matrix (0, nrow = high -low +1 , ncol= nstat)
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
    
    NSt <- matrix(0, nrow = n, ncol= nstat)
    for (s in (1:nstat)){
      St <- statMatrix[,s]
      m=mean(St)
      HSt=St-m
      HSt=as.matrix(HSt)
      NSt[,s] <- HSt / norm(HSt, type="f")
    }
    
    
    adj=adjacencygraph(trees, tip)
    ig <- graph.adjacency(adj, mode="undirected")
    LM=laplacian_matrix(ig, normalized = FALSE, weights = NULL,
                        sparse = TRUE)
    
    
    eyeL=eigs_sym(LM,1,"LM",opts = list(retvec = FALSE), lower=TRUE)
    eyeS=eigs_sym(LM,2,"SM",opts = list(retvec = FALSE),lower=TRUE)
    VL=eyeL$values[1]
    VS=eyeS$values[1]
    dist <- VL - VS
    cat(tip, file = filename, append=TRUE)
    for (s in (1:nstat)){
      result <- (t(NSt[,s]) %*% LM) %*% NSt[,s]
      result=as.numeric(result)
      # scale the result
      resolutionNew[tip-(low-1),s] <- (result - VS)/ dist 
      cat(", ", file = filename, append=TRUE)
      cat(resolutionNew[tip-(low-1),s], file = filename, append=TRUE)
    }
    cat("\n" , file = filename, append=TRUE)
  }  
}


#======================================================================

adjacencygraph = function (trees, tip) {
  n <- length(trees)
  
  # creates a hash table for (newickformat of each phylogenic tree with given leaves, an index in [1..n]) 
  ha <- new.env()
  for (i in (1:n)){
    tree1 <- getCanonical(trees[[i]],tip)
    Ctree2 <- findEqualsCanonical(tree1, tip)
    listTree2 <- extractEquals(Ctree2, tip)
    Ctree <- findEquals(trees[[i]], tip)
    listTree <- extractEquals(Ctree, tip)
    ll <-length(listTree)
    listTree[[ll+1]] <- tree1
    listTree=c(listTree,listTree2)
    for (k in listTree) {
      key <- extractNewick (k, tip+1, tip)
      ha[[key]] <- i
    }    
  }
  # creates the adjacency matrix
  I=numeric(2*n*(tip-2))
  J=numeric(2*n*(tip-2))
  K=numeric(2*n*(tip-2))
  t=1
 
  for (i in (1:n)){ 
    Spos=t
    row <- gsub("\\;","", gsub("(\\d*)(\\:)(\\d+)", "", write.tree(trees[[i]])))
    proned <- getNNI(trees[[i]], tip)
    for (j in (1:length(proned))){
      col <- proned[[j]]
      
      if (ha[[row]] != ha[[col]] && is.null(ha[[col]])==FALSE) {
        if(sum(J[Spos:t]==ha[[col]])==0 ){
          I[t]=ha[[row]]
          J[t]=ha[[col]]
          K[t]=1
          t=t+1
         
        }
      }
    }
   
  }
  
  ind=length(which(K==0))
  I=head(I,-ind)
  J=head(J,-ind)
  K=head(K,-ind)
  adj=sparseMatrix(I,J,x=K)
  return(adj)
}



