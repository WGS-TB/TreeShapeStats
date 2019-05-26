require(reshape2)

load("~/Trees/allTrees9.RData")

ntip=9
trees=allTree
n <- length(trees)
metric="NNI"
Ds=DistanceMatrix(trees,ntip,metric)
Ds <- Ds^2
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

Colless_ind=sapply(allTree, function(x) colless(as.treeshape(x),"pda"))
Sackin_ind=sapply(allTree, function(x) sackin(as.treeshape(x),"pda"))
#==============================================================================================
Colless_ind=sapply(allTree, function(x) colless(as.treeshape(x),"pda"))
df=cbind(treeLocation,Colless_ind)
df=as.data.frame(df)
names(df)=c("x","y","index")
Colless_dis=distribution.f(Colless_ind)
tree_col=character(n) 
tree_col[Colless_dis$avg]="red"
tree_col[Colless_dis$lower]="blue"
tree_col[Colless_dis$upper]="green"

#==============================================================================================
Saless_ind=0.11*Sackin_ind+Colless_ind
df=cbind(treeLocation,Saless_ind)
df=as.data.frame(df)
names(df)=c("x","y","index")
Saless_dis=distribution.f(Saless_ind)
tree_col=character(n) 
tree_col[Saless_dis$avg]="red"
tree_col[Saless_dis$lower]="blue"
tree_col[Saless_dis$upper]="green"

#==============================================================================================
Sackin_ind=sapply(trees, function(x) sackin(as.treeshape(x),"pda"))
df=cbind(treeLocation,Sackin_ind)
df=as.data.frame(df)
names(df)=c("x","y","index")
Sackin_dis=distribution.f(Sackin_ind)
tree_col=character(n) 
tree_col[Sackin_dis$avg]="Average"
tree_col[Sackin_dis$lower]="Low"
tree_col[Sackin_dis$upper]="High"

#==============================================================================================
I2_ind=sapply(trees, function(x) computeI2(x))
I2_dis=distribution.f(I2_ind)
tree_col=character(n) 
tree_col[I2_dis$avg]="Average"
tree_col[I2_dis$lower]="Low"
tree_col[I2_dis$upper]="High"

df=cbind(treeLocation,I2_ind)
df=as.data.frame(df)
names(df)=c("x","y","index")

#==============================================================================================
distribution.f <- function(x){
  low=as.numeric(quantile(x)[2])
  high=as.numeric(quantile(x)[4])
  list(avg=which(x>=low & x<=high),lower=which(x<low), upper=which(x>high))
}
#==============================================================================================
p=ggplot(data=df, aes(x = x, y = y)) +
  geom_point(stat="identity", position="identity")+geom_point(aes(colour = factor(tree_col)),show.legend = TRUE)+
  ggtitle("Embedding the set of trees with 9 tips using I2 statistics")+
  theme(plot.title = element_text(hjust = 0.5,size = 10),axis.text=element_text(size=6),axis.title=element_text(size=8))  
p+theme(legend.position="bottom", legend.title = 
         element_blank())


