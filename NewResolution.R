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
    indx <- which(edges[,1] == i)
    tree[i-tip,2] = edges[indx[1],2]
    tree[i-tip,3] = edges[indx[2],2]
  }
  return (tree)
}

#compute the number of tips of each subtree rooted at each internal node.
extractNTips = function (root, tip) {
  NTips <<-rep(1, 2*tip -1)
  for (i in ((tip+1): (2*tip -1))){
    NTips[i] <- length(Descendants(root, i, type=c("tips"))[[1]])
  }
  return (NTips)
}

#compute the number of Pitchforks of each subtree rooted at each internal node.
extractPitchforks = function (tree,  NTips, tip) {
  pitchforks <<-rep(0, 2*tip -1)
  for (i in (2*tip -1):(tip+1)){
    if(NTips[i]==3){ 
        pitchforks[i] <- 1
    }
    if(NTips[i]>3){ 
      pitchforks[i] <- pitchforks[tree[i-tip,2]]+pitchforks[tree[i-tip,3]]
    }
  }
  return (pitchforks)
}

#function to compute the canonical form of each tree
getCanonical = function (tree, NTips,tip){
  for (i in (1:(tip-1))){
    nodel <- tree[i,2]
    noder <- tree[i,3]
    if (NTips[nodel] < NTips[noder]){
      tree[i,2] <- noder
      tree[i,3] <- nodel
    }
  } 
  return (tree)
}

# extracts different sahpe of a tree based on rotation on 4th column of Ctree, and it returns as a list of trees
extractEquals = function (tree, NTips,tip){
  ls <- list()
  ls[[length(ls) +1]] <- tree
  for (i in (1:(tip-1))){
    nodel <- tree[i,2]
    noder <- tree[i,3]  
      if ((NTips[nodel] == NTips[noder]) & (NTips[nodel] > 2)){
        Als=ls
        for ( p in ls ){
          t <- p[i,2] 
          p[i,2] <- p[i,3]
          p[i,3] <- t 
         Als[[length(Als) +1]] <- p
        }
        ls=Als
      }
   
  } 
 
  return (ls)
}

#find the set of canonical trees that are in distance 1 from a specific tree
getNNI = function (root, tip) {
  NTips=extractNTips(root, tip)
  neighbours <- nni(root)
  proned=lapply(neighbours,function(x,tip) setTree(x,tip), tip=tip)
  ntips=lapply(neighbours,function(x,tip) extractNTips(x, tip), tip=tip)
  proned=mapply(function(x,y,tip) getCanonical(x,y,tip), proned,ntips, tip=tip, SIMPLIFY = FALSE)
  proned=lapply(proned,function(x, y, tip) extractNewick(x, y, tip), y=tip+1,tip=tip)
  proned=unique(proned, incomparables = FALSE)
  return (proned)
}

#find the newick format of a tree
extractNewick = function(tree, node, tip){
  if (node <= tip){
    return ("")
  }
  i <- which(tree[,1]==node)
  return (paste ("(", extractNewick(tree, tree[i,2], tip), "," ,extractNewick(tree, tree[i,3], tip), ")", sep="" ))
}


#compute the I2 statistics 
statI2 = function (tree, NTips,tip, n){
  
  I2 <- 0
  for (i in ((tip+1):(2*tip-1))){
    
    ri=NTips[tree[i-tip,2]]
    si=NTips[tree[i-tip,3]]
    
    if ((ri + si) > 2) {
      I2 <- I2 + (abs(ri -si)/abs (ri +si -2)) 
    }        
  }
  I2 <- I2 / (n -2)
  return (I2)
}

#compute the TP statistics
statTP=function (tree,NTips, pitchforks,tip){
  sum <- 0
  for (i in ((tip+1):(2*tip -1))){ 
    ch1 <- tree[i-tip,2]
    ch2 <- tree[i-tip,3]
    
    sum = sum + (NTips[ch1] - NTips[ch2])^2+(pitchforks[ch1] - pitchforks[ch2])^2
    
  } 
  return (sum)
  
}

#The main function to compute the scaled resolution
main = function(low,high,nstat){
  outputdir <- "/Users/maryam/Google Drive/Research/Tree Stat/Codes/New"
  filename <- paste(outputdir,"/resolution.csv",sep='')
  cat(", Ic, Sackin, Variance, I2, B1, B2, Saless,TP\n", file = filename, append = TRUE)
  inputdir <-  "/Users/maryam/Google Drive/Research/Tree Stat/Codes/Trees"
  resolutionNew <- matrix (0, nrow = high -low +1 , ncol= nstat)
  for (tip in (low:high)) {
    
    # reading all phylogenic trees for a given number of tips from a file
    inputfile <- paste(inputdir,"/tr",tip,".tre",sep="")
    trees <- read.tree(inputfile, keep.multi = TRUE)
    # fill the Stat Matrix 
    n <- length(trees)
    ha <- new.env()
    statMatrix <- matrix(0, nrow = n, ncol= nstat) 
    for (j in (1:n)) {
      root <- trees[[j]]
      tree <- setTree(root, tip)
      parents <- Ancestors(root, 1:(2*tip-1), type=c("parent"))
      Nis<-node.depth.edgelength(root)
      MHat <- node.depth(root,2)-1
      NTips <- extractNTips(root, tip)
      pitchforks <- extractPitchforks(tree, NTips,tip ) 
     
      # calculate Ic
      statMatrix[j,1] <- colless(as.treeshape(root))
      
      # calculate  Sackin
      statMatrix[j,2] <- sackin(as.treeshape(root),norm = NULL)
      
      # calculate stat Variance
      statMatrix[j,3] <- var(Nis[1:tip]) 
      
      # calculate stat I2
      statMatrix[j,4] <- statI2(tree, NTips,tip, n)
      
      # calculate stat B1
      statMatrix[j,5] <- sum(apply(as.matrix(MHat[(tip+2):(2*tip-1)]),1,function(x) 1/x))
      
      # calculate stat B2
      statMatrix[j,6] <- sum(apply(as.matrix(Nis[1:tip]),1,function(x) x/2^x))
      
      # calculate stat Saless: combination of sackin and colless on each node for all nodes in tree
      statMatrix[j,7] <- 1.300745*statMatrix[j,2]+statMatrix[j,1]
      
      # calculate stat TCP: combination of number of Tips and Pitchforks for all internal nodes. 
      statMatrix[j,8] <- statTP(tree,NTips, pitchforks,tip)
  
    # creates a hash table for (newickformat of each phylogenic tree with given leaves, an index in [1..n])
    tree1 <- getCanonical(tree, NTips,tip)
    listTree1 <- extractEquals(tree, NTips,tip)
    listTree1[[length(listTree1)+1]]=tree1
    listTree2 <- extractEquals(tree1, NTips,tip)
    listTree=c(listTree1,listTree2)
    for (k in listTree) {
      key <- extractNewick (k, tip+1, tip)
      ha[[key]] <- j
    }    
   }
  
    NSt <- matrix(0, nrow = n, ncol= nstat)
    for (s in (1:nstat)){
      St <- statMatrix[,s]
      m=mean(St)
      HSt=St-m
      HSt=as.matrix(HSt)
      NSt[,s] <- HSt / norm(HSt, type="f")
    }
    
    
    adj=adjacencygraph(trees, ha,tip)
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

adjacencygraph = function (trees, ha,tip) {
  n <- length(trees)
  
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



