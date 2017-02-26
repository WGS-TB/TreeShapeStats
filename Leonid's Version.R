initDir = getwd()
setwd("/Users/admin/Documents/HSPHwork/Project with Caroline/")
source("ArtPaperScript.R")
require(Matrix)
require(apTreeshape)
SIZE = 20
statNames = outer(c("Between", "Closeness", "Eigen"), c("Min", "Max"), paste0)
statNames = c(statNames, c("Diams", "Wieners", "Eigenvalues"))
statNames = c(statNames, c("colless", "sackin", "variance", "I2", "B1", "B2"))
load(paste0("allTrees", SIZE, ".RData")) # load the trees, assumed to be precomputed
load(paste0("allStats", SIZE, ".RData")) # assumes the first 9 stats are precomputed
setwd(initDir)

precomputeGrid = function(size = SIZE) {
  A = matrix(NA, size, size)
  diag(A) = allSizes[1:size]
  B = matrix(0, size, size)
  for (ind1 in 2:size) {
    curMin = ceiling(ind1/2)
    curRange = (ind1 - 1):curMin
    for (ind2 in curRange) {
      if (ind2 * 2 == ind1) {
          A[ind1, ind2] = choose(allSizes[ind2] + 1, 2)
      }
      else {
        A[ind1, ind2] = allSizes[ind2] * allSizes[ind1 - ind2]
      }
    }
    B[ind1, curRange] = cumsum(A[ind1, curRange])
  }
  B
}

treeNumber = function(tree) {
  n = Ntip(tree)
  edges = tree$edge
  numbering = matrix(NA, 2 * n - 1, 2)
  numbering[1:n, ] = cbind(rep(1, n), rep(0, n))
  for (index in Nedge(tree):1) {
    curEdge = edges[index,]
    curFather = curEdge[1]
    curSon = curEdge[2]
    if (is.na(numbering[curFather, 1])) {
      numbering[curFather, ] = numbering[curSon, ] # first child
    }
    else {
      first = numbering[curFather, ]
      second = numbering[curSon, ]
      numbering[curFather, ] = combineNumberings(first, second)
    }
  }
  numbering[n + 1, 2]
}

combineNumberings = function(first, second) {
  if (first[1] < second[1] || (first[1] == second[1] && first[2] < second[2])) { # switch
    temp = first
    first = second
    second = temp
  }
  x1 = first[1]
  x2 = second[1]
  y1 = first[2]
  y2 = second[2]
  comboX = x1 + x2
  if (x1 > x2) {
    comboY = Grid[comboX, x1 + 1] + allSizes[x2] * y1 + y2
  }
  else {
    comboY = Grid[comboX, x1 + 1] + choose(allSizes[x2] + 1, 2) - choose(y1 + 2, 2) + y2
  }
  c(comboX, comboY)
}

createCayleyGraphLaplacian = function(size = SIZE) { # assumes Grid, allTree are global!
  totalLength = numTrees * 2 * (size - 2)
  I = c(0:(numTrees - 1), rep(NA, totalLength))
  J = c(0:(numTrees - 1), rep(NA, totalLength))
  Order = rep(NA, numTrees)
  Deg = rep(NA, numTrees)
  pos = numTrees
  print(paste("There are", numTrees, "trees to process"))
  for (ind in 1:numTrees) {
    if (ind %% 10000 == 0) {  print(ind)  }
    curTree = allTree[[ind]]
    curIndex = treeNumber(curTree)
    curNeighbors = nni(curTree)
    curNeighbors = lapply(curNeighbors, reorder)
    curIndices = sapply(curNeighbors, treeNumber)
    curIndices = setdiff(unique(curIndices), curIndex)
    curLen = length(curIndices)
    Order[curIndex + 1] = ind
    Deg[curIndex + 1] = curLen
    I[pos + (1:curLen)] = curIndex
    J[pos + (1:curLen)] = curIndices
    pos = pos + curLen
  }
  K = c(Deg, rep(-1, totalLength))
  Lap = sparseMatrix(I[1:pos], J[1:pos], x = K[1:pos], dims = c(numTrees, numTrees), index1 = FALSE)
  output = list(Lap, Order)
  output
}

computeStats = function(size = SIZE) {
  allStats = matrix(NA, numTrees, 15, dimnames = list(NULL, statNames))
  for (ind in 1:9) {
    allStats[, ind] = get(paste0("all", statNames[ind]))
  }
  allTreeShapes = lapply(allTree, as.treeshape)
  allStats[,"colless"]  = sapply(allTreeShapes, colless)
  allStats[,"sackin" ]  = sapply(allTreeShapes, sackin )
  allStats[,"variance"] = sapply(allTree, function(x) {
    var(node.depth.edgelength(x)[1:size])
  })
  allStats[,"I2"] = sapply(allTree, computeI2)
  allStats[,"B1"] = sapply(allTree, computeB1)
  allStats[,"B2"] = sapply(allTree, computeB2)
  allStats
}

computeResolution = function(size = SIZE) {
  a0 = proc.time()
  numTrees <<- allSizes[size]
  Grid <<- precomputeGrid(size)
  a1 = proc.time()
  print(a1 - a0)
  print("Computing the Laplacian")
  prep = createCayleyGraphLaplacian(size)
  LapMatrix = prep[[1]]
  print(nnzero(LapMatrix))
  print(object.size(LapMatrix))
  a2 = proc.time()
  print(a2 - a1)
  print("Computing the statistics")
  allStats = computeStats(size)
  allStats = allStats[prep[[2]], ]
  numStats = ncol(allStats)
  allRes = rep(NA, numStats)
  names(allRes) = colnames(allStats)
  a3 = proc.time()
  print(a3 - a2)
  print("Computing the eigenvalues")
  Max = eigs_sym(LapMatrix, 1, "LM", opts = list(retvec = FALSE), lower = TRUE)$values[1]
  Min = eigs_sym(LapMatrix, 2, "SM", opts = list(retvec = FALSE), lower = TRUE)$values[1]
  a4 = proc.time()
  print(a4 - a3)
  print("Computing the resolutions")
  for (ind in 1:numStats) {
    curStat = allStats[,ind, drop = FALSE]
    curStat = curStat - mean(curStat)
    curStat = curStat/norm(curStat, "F")
    allRes[ind] = ((t(curStat) %*% (LapMatrix %*% curStat)) - Min) / (Max - Min)
  }
  a5 = proc.time()
  print(a5 - a4)
  allRes
}