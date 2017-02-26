source("SpectralComparison.R")
source("SpectralTrees.R")
source("PlottingScript.R")
source("Multiplot.R")
Dir = ""
if (getwd() == "/Users/admin/Documents/HSPHwork/Project with Caroline") {
  Dir = "/Users/admin/Documents/HSPHwork/TedProject/"
}
source(paste0(Dir, "NeighborJoin.R"))
require(RPANDA)
require(ggplot2)
RANGE = 5:20

driver = function() {
  for (SIZE in RANGE) {
    print(SIZE)
    curSize = allSizes[SIZE]
    allTree <<- vector("list", curSize)
    print(paste("There are", curSize, "trees to process"))
    setupLinspace(SIZE - 2)
    generateTrees(p = 1, s = 0, cL = 0, n = SIZE - 2)
    save(allTree, file = paste0("allTrees", SIZE, ".RData"))
    allDiams = rep(NA, curSize)
    allWieners = rep(NA, curSize)
    allBetweenMax = rep(NA, curSize)
    allBetweenMin = rep(NA, curSize)
    allClosenessMax = rep(NA, curSize)
    allClosenessMin = rep(NA, curSize)
    allEigenMax = rep(NA, curSize)
    allEigenMin = rep(NA, curSize)
    allEigenvalues = rep(NA, curSize)
    for (ind in 1:curSize) {
      if (ind %% 10000 == 0) {
        print(ind)
      }
      curTree = allTree[[ind]]
      allDiams[ind] = computeDiameter(curTree)
      allWieners[ind] = computeWienerIndex(curTree, double = FALSE)
      curBetween = computeBetweenness(curTree)
      allBetweenMin[ind] = min(curBetween[(SIZE+1):(length(curBetween))])
      allBetweenMax[ind] = max(curBetween)
      curCloseness = computeFarness(curTree)
      allClosenessMax[ind] = min(curCloseness)
      allClosenessMin[ind] = max(curCloseness)
      curEigen = computeEigenvector(curTree)
      allEigenvalues[ind] = curEigen[[2]]
      curEigenvector = curEigen[[1]]
      allEigenMax[ind] = max(curEigenvector)
      allEigenMin[ind] = min(curEigenvector)
    }
    save(allDiams, allWieners, allBetweenMin, allBetweenMax, allClosenessMin, allClosenessMax, 
         allEigenMin, allEigenMax, allEigenvalues, file = paste0("allStats", SIZE, ".RData"))
  }
}

# minDiams = vector("list", length(RANGE))
# maxDiams = vector("list", length(RANGE))
# minWieners = vector("list", length(RANGE))
# maxWieners = vector("list", length(RANGE))
# minBetweenMax = vector("list", length(RANGE))
# maxBetweenMax = vector("list", length(RANGE))
# minClosenessMax = vector("list", length(RANGE))
# maxClosenessMax = vector("list", length(RANGE))
# minEigenMax = vector("list", length(RANGE))
# maxEigenMax = vector("list", length(RANGE))
# minEigenMaxNorm = vector("list", length(RANGE))
# maxEigenMaxNorm = vector("list", length(RANGE))
# minEigenvalues = vector("list", length(RANGE))
# maxEigenvalues = vector("list", length(RANGE))
# rangeDiams = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeWieners = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeBetweenMax = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeClosenessMax = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeEigenMax = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeEigenMaxNorm = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# rangeEigenvalues = matrix(NA, length(RANGE), 2, dimnames = list(RANGE, c("min", "max")))
# for (SIZE in RANGE) {
#   print(SIZE)
#   load(paste0("allTrees", SIZE, ".RData"))
#   load(paste0("allStats", SIZE, ".RData"))
#   rangeDiams[as.character(SIZE),] = range(allDiams)
#   minWieners[[SIZE]] = allTree[[which.min(allWieners)]]
#   maxWieners[[SIZE]] = allTree[[which.max(allWieners)]]
#   rangeWieners[as.character(SIZE),] = range(allWieners)
#   if (allDiams[which.min(allWieners)] != min(allDiams)) {
#     print(paste("Found a counterexample", SIZE))
#   }
#   minDiams[[SIZE]] = allTree[[which.min(allWieners)]]
#   maxDiams[[SIZE]] = allTree[[which.max(allDiams)]]
#   rangeBetweenMax[as.character(SIZE),] = range(allBetweenMax)
#   minBetweenMax[[SIZE]] = allTree[[which.min(allBetweenMax)]]
#   maxBetweenMax[[SIZE]] = allTree[[which.max(allBetweenMax)]]
#   rangeClosenessMax[as.character(SIZE),] = range(allClosenessMax)
#   minClosenessMax[[SIZE]] = allTree[[which.min(allClosenessMax)]]
#   maxClosenessMax[[SIZE]] = allTree[[which.max(allClosenessMax)]]
#   rangeEigenMax[as.character(SIZE),] = range(allEigenMax)
#   rangeEigenMaxNorm[as.character(SIZE),] = range(allEigenMax/allEigenMin)
#   minEigenMax[[SIZE]] = allTree[[which.min(allEigenMax)]]
#   maxEigenMax[[SIZE]] = allTree[[which.max(allEigenMax)]]
#   minEigenMaxNorm[[SIZE]] = allTree[[which.min(allEigenMax/allEigenMin)]]
#   maxEigenMaxNorm[[SIZE]] = allTree[[which.max(allEigenMax/allEigenMin)]]
#   rangeEigenvalues[as.character(SIZE),] = range(allEigenvalues)
#   minEigenvalues[[SIZE]] = allTree[[which.min(allEigenvalues)]]
#   maxEigenvalues[[SIZE]] = allTree[[which.max(allEigenvalues)]]
# }

findSameEigenvector = function(SIZE = 18, tol = 1e-9) {
  load(paste0("allTrees", SIZE, ".RData"))
  load(paste0("allStats", SIZE, ".RData"))
  allE=cbind(allEigenMin,allEigenMax, 1:allSizes[SIZE])
  orderE=order(allE[,1], allE[,2])
  allES=allE[orderE,]
  goodInds=which(diff(allES[,1]) < tol & diff(allES[,2]) < tol)
  for (ind in goodInds) {
    cur1=allES[ind,3]
    cur2=allES[ind+1,3]
    tr1=allTree[[cur1]]
    tr2=allTree[[cur2]]
    spec1=computeEigenvector(tr1)[[1]]
    spec2=computeEigenvector(tr2)[[1]]
    u=sum(abs(sort(spec1)-sort(spec2)))
    if(u < tol) {
      print(paste(cur1,cur2,u))
      return(list(tr1, tr2))
    }
  }
}

findSameVector = function(SIZE = 9, metrics = "Farness") {
  load(paste0("allTrees", SIZE, ".RData"))
  load(paste0("allStats", SIZE, ".RData"))
  curSize = allSizes[SIZE]
  Lengths = ifelse(metrics == "Betweenness", SIZE, 2 * SIZE - 1)
  Length = sum(Lengths)
  Mat = matrix(NA, curSize, Length)
  for (ind in 1:curSize) {
    curTree = allTree[[ind]]
    pos = 0
    for (index in 1:length(metrics)) {
      curFun = get(paste0("compute", metrics[index]))
      curVec = curFun(curTree)
      Mat[ind, pos + (1:Lengths[index])] = tail(curVec, Lengths[index])
      pos = pos + Lengths[index]
    }
  }
  pos = 0
  for (index in 1:length(metrics)) {
    curRange = pos + (1:Lengths[index])
    Mat[,curRange] = t(apply(Mat[,curRange], 1, sort))
    pos = pos + Lengths[index]
  }
  dup = which(duplicated(Mat))
  if (length(dup) > 0) {
    dup = lapply(dup, function(x) {
      allTree[which(rowSums(Mat==matrix(Mat[x,], curSize, Length, byrow = TRUE))==Length)]
    })
  }
  dup
}

### This function produces the figure containing the three counterexamples
makeCounterexampleFigure = function() {
  CounterB = findSameVector(9, "Betweenness")
  B1 = CounterB[[1]][[1]]
  B2 = CounterB[[1]][[2]]
  B2r = rotate(rotate(rotate(rotate(rotate(B2,16),15),10),15),16)
  BVec1 = computeBetweenness(B1)
  BVec2 = computeBetweenness(B2)
  BMat = cbind(as.character(c(9:6,1:3,5:4)), as.character(c(5:1,6:7,9:8)))
  objectB = cophylo(B1, B2r, assoc = BMat)
  pdf("CounterB.pdf")
  plot(objectB, ftype = "off")
  myInds = c(10, 15, 16, 17, 11, 12, 13, 14)
  nodelabels.cophylo(BVec1[-(1:9)], which = "left" , cex = 0.9, bg = "white")
  nodelabels.cophylo(BVec2[myInds], which = "right", cex = 0.9, bg = "white")
  tiplabels.cophylo(BVec1[(1:9)], which = "left" , cex = 0.9, bg = "white")
  tiplabels.cophylo(BVec2[(1:9)], which = "right", cex = 0.9, bg = "white")
  dev.off()
  CounterC = findSameVector(13, "Farness")
  C1 = CounterC[[1]][[1]]
  C2 = CounterC[[1]][[2]]
  C1r = rotate(rotate(rotate(C1,21),20),19)
  CVec1 = computeFarness(C1)
  CVec2 = computeFarness(C2)
  CInds1 = c(13:8,1:3,7:4)
  CMat = cbind(as.character(CInds1), as.character(c(10:11,13:12,9:8,1:3,7:4))) 
  objectC = cophylo(C1r, C2, assoc = CMat)
  tiplabels.cophylo(CVec1[1:13], which = "left" , cex = 0.9, bg = "white")
  tiplabels.cophylo(CVec2[1:13], which = "right", cex = 0.9, bg = "white")
  plot(objectC, ftype = "off")
  CounterE = findSameEigenvector(18)
  
}

checkWienerIndices = function(N = 30, offset = 2) {
  Mins  = rep(NA,N)
  Maxes = rep(NA,N)
  Means = rep(NA,N)
  for (ind in 1:N) {
    cur = ind + offset
    load(paste0("TracesAllNodes/TracesAllNodes",cur,".Rdata"))
    curValues = as.numeric(names(curTracesRed))/2
    curFreqs  = as.vector(curTracesRed)
    Mins[ind]  = min(curValues)
    Maxes[ind] = max(curValues)
    Means[ind] = sum(curValues * curFreqs)/sum(curFreqs)
    if (cur %in% 2^(1:5)) {
      allValues = c(curValues, max(curValues) + 1)
      myhist <-list(breaks = allValues, counts = curFreqs, density = curFreqs/diff(allValues),
      xname = paste("Wiener indices for phylogenetic trees on", cur, "leaves"))
      class(myhist) <- "histogram"
      pdf(paste0("WienerHistogram",cur,".pdf"))
      plot(myhist)
      dev.off()
    }
  }
  output = list(Mins, Maxes, Means)
  output
}

checkCounterexamples = function() {
  Tree = mergeSubtrees(Cherry, Singleton)
  for (ind in 1:3) {
    baseLengths = rbind(c(4.124905584084558159332936, 1.777414646626160824397943, 
                          3.087661446451583894502595, 0.947565530795418107055579),
                        c(1.812363501995566708609640, 2.785812690566482823514912, 
                          4.282526048434207590995765, 1.213116008005954294179023),
                        c(3.513639916017646947222228, 0.1675441776785849445448657,
                          1.316446674855575760995505, 0.245549273439368873692737))
    altLengths = rbind(c(4, 2, 3, 1), c(2, 3, 4, 1), c(2, 2, 1, 1))
    Tree1 = Tree
    Tree1$edge.length = baseLengths[ind,]
    Vec1 = spectR(Tree1)$eigenvalues
    Tree2 = Tree
    Tree2$edge.length = altLengths[ind,]
    Vec2 = spectR(Tree2)$eigenvalues
    print(paste("Example", ind, "the maximum difference between spectra is", max(abs(Vec1 - Vec2))))
  }
}

### This function creates the first two figures in the paper - the trees are the same as above
createFigures = function() {
  Tr = createEmptyTree(numNodes = 13)
  Tr=createNewNode(Tr, c(4,5), c(1,1), 8)	
  Tr=createNewNode(Tr, c(6,7), c(1,1), 9)
  Tr=createNewNode(Tr, c(2,3),  c(2,2), 10)
  Tr=createNewNode(Tr, c(8,9), c(1,1), 11)
  Tr=createNewNode(Tr, c(10,11), c(2,2), 12)
  Tr=createNewNode(Tr, c(1,12), c(8,4), 13)
  Tree=convertPhylo(Tr)
  Tree$tip.label=c("a","f","g","d","e","b","c")
  Tree$node.label=LETTERS[1:6]
  pdf("fig1a.pdf")
  plot.phylo(Tree, show.node.label = TRUE)
  dev.off()
  Tree1 = Tree
  Tree1$edge.length = rep(1,12)
  pdf("fig1b.pdf")
  plot.phylo(Tree1, show.node.label = TRUE)
  dev.off()
  output = list(Tree, Tree1)
}

produceResultFigures = function() {
  Res00 = metaDriverOld("hivdenmeas.Rdata", decreasing = FALSE)
  save(Res00, file = "ThreeVirusesLowestEigs.Rdata")
  Res01 = metaDriverOld("hivdenmeas.Rdata", decreasing = TRUE)
  save(Res01, file = "ThreeVirusesHighestEigs.Rdata")
  Res10 = metaDriverNew("hivdenmeas.Rdata", decreasing = FALSE)
  save(Res10, file = "ThreeVirusesLowestMetrics.Rdata")
  Res11 = metaDriverNew("hivdenmeas.Rdata", decreasing = TRUE)
  save(Res11, file = "ThreeVirusesHighestMetrics.Rdata")
  plotResults("ThreeVirusesLowestEigs.Rdata", "Res00", metric = FALSE, unitary = FALSE, test = FALSE)
  plotResults("ThreeVirusesHighestEigs.Rdata", "Res01", metric = FALSE, unitary = FALSE, test = FALSE)
  plotResults("ThreeVirusesLowestMetrics.Rdata", "Res10", metric = TRUE, unitary = FALSE, test = FALSE)
  plotResults("ThreeVirusesHighestMetrics.Rdata", "Res11", metric = TRUE, unitary = FALSE, test = FALSE)
  plotResults("ThreeVirusesLowestMetrics.Rdata", "Res10", metric = TRUE, unitary = TRUE, test = FALSE)
  plotResults("ThreeVirusesHighestMetrics.Rdata", "Res11", metric = TRUE, unitary = TRUE, test = FALSE)
  e = new.env()
  load("globtrees120.Rdata", e)
  load("fiveyrtrees120.Rdata", e)
  load("usatrees120.Rdata" , e)
  setwd("chronic")
  abbrs = c("HEPCV1a", "HEPCV1b", "HEPCV5a", "HDELG1", "HIVB1", "HIVC1")
  prepareFile(list.files(), abbrs, "ChronicTrees.Rdata")
  load("ChronicTrees.Rdata")
  setwd("../")
  fiveyrflu = get("fiveyrtrees", envir = e)
  globalflu = get("globtrees", envir = e)
  usaflu = get("usatrees", envir = e)
  fullAbbrs = c("fiveyrflu", "globalflu", abbrs, "usaflu")
  save(list=fullAbbrs, file = "NineViruses.Rdata")
  Res00 = metaDriverOld("NineViruses.Rdata", decreasing = FALSE)
  save(Res00, file = "NineVirusesLowestEigs.Rdata")
  Res01 = metaDriverOld("NineViruses.Rdata", decreasing = TRUE)
  save(Res01, file = "NineVirusesHighestEigs.Rdata")
  Res10 = metaDriverNew("NineViruses.Rdata", decreasing = FALSE)
  save(Res10, file = "NineVirusesLowestMetrics.Rdata")
  Res11 = metaDriverNew("NineViruses.Rdata", decreasing = TRUE)
  save(Res11, file = "NineVirusesHighestMetrics.Rdata")
  plotResults("NineVirusesLowestEigs.Rdata", "Res00", metric = FALSE, unitary = FALSE, test = FALSE)
  plotResults("NineVirusesHighestEigs.Rdata", "Res01", metric = FALSE, unitary = FALSE, test = FALSE)
  plotResults("NineVirusesLowestMetrics.Rdata", "Res10", metric = TRUE, unitary = FALSE, test = FALSE)
  plotResults("NineVirusesHighestMetrics.Rdata", "Res11", metric = TRUE, unitary = FALSE, test = FALSE)
  plotResults("NineVirusesLowestMetrics.Rdata", "Res10", metric = TRUE, unitary = TRUE, test = FALSE)
  plotResults("NineVirusesHighestMetrics.Rdata", "Res11", metric = TRUE, unitary = TRUE, test = FALSE)
}