CORENUM = 7     # Core used here

Dir = ""
if (getwd() == "/Users/admin/Documents/HSPHwork/Project with Caroline") {
	Dir = "/Users/admin/Documents/HSPHwork/TedProject/"
}

source(paste0(Dir, "Alternative.R"))
source(paste0(Dir, "NeighborJoin.R"))
require(combinat)
require(rARPACK)

Singleton = list(Nnode = 0, tip.label = "1", edge = matrix(NA, nrow = 0, ncol = 2))
class(Singleton) = "phylo"
Cherry = list(Nnode = 1, tip.label = c("1", "2"), edge = rbind(c(3,1), c(3,2)), edge.length = rep(1, 2))
class(Cherry) = "phylo"

allSizes = c(1, 1, 1, 2, 3, 6, 11, 23, 46, 98, 207, 451, 983, 2179, 4850, 10905, 24631, 56011, 127912, 
    293547, 676157, 1563372, 3626149, 8436379, 19680277, 46026618, 107890609, 253450711, 
    596572387, 1406818759, 3323236238, 7862958391) # Wedderburn - Etherington numbers, OEIS:A001190)
specialPairs = matrix(NA, nrow = 2, ncol = 0)
specialInds = NULL
specialTrees = list()

LAP  = TRUE		  # TRUE if the spectra computed are Laplacian spectra
DIST = TRUE     # TRUE if the spectra computed are distance spectra
FULL = TRUE		  # TRUE if the spectra computed are node-node spectra
WIENER = TRUE		# TRUE if the Wiener index of each tree is calculated (only makes sense if DIST = TRUE and FULL = TRUE)
DOUBLE = TRUE		# TRUE if the Wiener index includes each pair of nodes twice (only makes sense if WIENER = TRUE)
ROUND = 4			  # Number of digits to which eigenvalues are rounded
SIZE = 6 #32	  # Number of leaves in the trees being analyzed
CORES = 10      # Number of cores (which could be virtual)
GENONLY = TRUE  # TRUE if the process is only used to generate trees

if (GENONLY) { allTree <<- vector("list", allSizes[SIZE]) }

### This function is a driver that performs the search co-spectral trees of specified spectral type
driver = function(size = SIZE, num = CORENUM) {
	filename = paste0("SpecialPairs", size, "X", num, ".Rdata")
	if (!(filename %in% list.files())) {
		generateSpectraLinspace(size = size, num = num)
		postprocessAllSpectra(size = size)
	}
	load(filename)
	specialInds <<- sort(unique(as.vector(specialPairs)))
	print(paste("There are", ncol(specialPairs), "collisions;", length(specialInds), "colliding trees"))
	if (length(specialInds) > 0) {
		WIENER <<- FALSE # no longer useful at this point!
		pos <<- 1
		nextInd <<- specialInds[pos]
		generateSpectraLinspace(size = size)
		Dists = postprocessSpecialTrees(specialPairs)
		print(paste("The closest distance is", min(Dists)))
		Inds = specialPairs[,which.min(Dists)]
		for (index in 1:2) {
		  assign(paste0("Tr", index), specialTrees[[as.character(Inds[index])]]) 
		  assign(paste0("Mat", index), computeMatrix(get(paste0("Tr", index)), full = FULL, lap = LAP, dist = DIST))
		  write.table(get(paste0("Mat", index)), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",", file = paste0("Matrix", size, "X", num, "I", index, ".csv"))
		}
		save(Tr1, Tr2, file = paste0("SpecialTrees", size, "X", num, ".Rdata"))
		save(Inds, file = paste0("BestInds", size, "X", num, ".Rdata"))
	}
}

### This function is a setup function that prepares the global environment for generating trees fast
setupLinspace = function(n) {
  Par <<- rep(NA, n)
  Child <<- rep(0, n + 1)
  WIndex <<- 0
  if (WIENER) {
    Sub <<- 2 * (n:1) + 1
    WIndex <<- (1 + DOUBLE) * (sum(Sub * (2 * n + 3 - Sub)) + (n + 2) * (2 * n + 2))
    allProducts <<- (1 + DOUBLE) * (1:(2 * n + 3)) * ((2 * n + 2):0)
    increases <<- c(allProducts[3:(2 * n + 3)], 0, 0) - allProducts
    decreases <<- c(0, 0, allProducts[1:(2 * n + 1)]) - allProducts 
  }
  useInds <<- ifelse(rep(FULL, 3), c(1, 2 * n + 2, 2 * n + 3), c(1, n + 1, n + 2))
  index <<- 1
}

### This function generates the spectra of required type for all size n trees using only O(n) space
generateSpectraLinspace = function(size = SIZE, num = CORENUM) {
  setupLinspace(n = size - 2)
	numRows = allSizes[size]
	numCols = 3
	if (LAP && DIST && is.null(specialInds) && !exists("Q")) {
	  load(paste0(ifelse(FULL, "TracesAllNodes/TracesAllNodes", "TracesLeavesOnly/TracesLeavesOnly"), size, ".Rdata"))
	  fullTable = splitFrequencyMatrix(curTracesRed, CORES)[[CORENUM]]
	}
	else {
	  fullTable = ifelse(length(specialInds) > 0, length(specialInds), allSizes[size])
	  names(fullTable) = 0
	}
	goodTraces = fullTable[fullTable > 1]
	allTr <<- as.numeric(names(goodTraces))
	numTraces = length(goodTraces)
	print(paste("There are", numTraces, "traces to process"))
	allSpectra <<- vector("list", numTraces)
	allIndices <<- rep(1, numTraces)
	names(allSpectra) <<- names(goodTraces)
	names(allIndices) <<- names(goodTraces)		
	for (ind in 1:numTraces) {
	  allSpectra[[ind]] <<- matrix(0, goodTraces[ind], numCols)
	}
	print(paste("There are", numRows, "trees to process"))
	generateTrees(p = 1, s = 0, cL = 0, n = size - 2)
}

### This function postprocesses the spectra to identify any matching ones for further investigation
postprocessAllSpectra = function(size = SIZE, num = CORENUM) {
	print("Postprocessing all spectra")
  Len = length(allSpectra)
  print(paste("There are", Len, "traces to process"))
  for (ind in 1:Len) {
    if (ind %% 100 == 0) {
      print(ind)
    }
    curSpectra = allSpectra[[ind]]
    extraMat = findMatchingRows(curSpectra, byLastCol = TRUE)
    specialPairs = cbind(specialPairs, t(extraMat))	
  }
	save(specialPairs, file=paste0("SpecialPairs", size, "X", num, ".Rdata"))
}

### This function re-generates trees based on indices obtained from a search for matching spectra
postprocessSpecialTrees = function(specialPairs) {
  print("Postprocessing special trees")
  specialSpectra = lapply(specialTrees, function(x) {computeSpectrum(x, lap = LAP, full = FULL, trace = TRUE, dist = DIST, K = 0)})
  numSpecial = ncol(specialPairs)
  print(paste("There are", numSpecial, "pairs to process"))
  allDists = rep(NA, numSpecial)
  for (ind in 1:numSpecial) {
    if (ind %% 10000 == 0) {
      print(ind)
    }
    ind1 = specialPairs[1,ind]
    ind2 = specialPairs[2,ind]
    allDists[ind] = sum(abs(specialSpectra[[as.character(ind1)]] - specialSpectra[[as.character(ind2)]]))
  }
  return(allDists)
}

### This function identifies any identical spectra, optionally enforcing a match in the last column
findMatchingRows = function(curSpectra, byLastCol = TRUE) {
	extraMat = matrix(NA, nrow = 0, ncol = 2)
	if (byLastCol) {
		numCols = ncol(curSpectra)
		Z = extractUniqueRows(curSpectra[,-numCols])
	}
	else {
		Z = extractUniqueRows(curSpectra)
	}
	allGroups = Z[[2]]
	allLens = sapply(allGroups, length)
	goodGroups = allGroups[allLens > 1]
	Len = length(goodGroups)
	if (Len > 0) {
		if (byLastCol) {
			curIndices = curSpectra[,numCols]
			extraMat = lapply(goodGroups, function(x) {combn2(curIndices[x])})
		}
		else {
			extraMat = lapply(goodGroups, combn2)
		}
		extraMat = do.call(rbind, extraMat)
	}
	extraMat
}

### This function identifies all pairs based on groups, i.e. computes a graph's transitive closure
postprocessSpecialPairs = function(size = SIZE, num = CORENUM) {
	splitPairs = split(specialPairs[1,], specialPairs[2,])
	roots = as.numeric(names(splitPairs))
	fullSplitPairs = lapply(1:length(roots), function(x) {combn2(c(roots[x], splitPairs[[x]]))})
	specialPairs = t(do.call(rbind, fullSplitPairs))
	save(specialPairs, file=paste0("SpecialPairs", size, "X", num, ".Rdata"))
}

### This function, based on Gang Li's thesis, generates trees of size n in amortized constant time
generateTrees = function(p, s, cL, n) {
	if (p > n) { 
	  if (GENONLY) { allTree[[index]] <<- extractTree(Par, Child, n); index <<- index + 1 }
	  else { processCurTree(n) }
	}
	else {
		if (cL == 0) {
			Par[p] <<- p - 1
		}
		else {
			if (Par[p - cL] < s) {
				if (WIENER) { old = Par[p]; new = Par[s]; updateSubVector(old, new) }
				Par[p] <<- Par[s]
			}
			else {
				if (WIENER) { old = Par[p]; new = cL + Par[p - cL]; updateSubVector(old, new) }
				Par[p] <<- cL + Par[p - cL]
			}
		}
		Child[Par[p] + 1] <<- Child[Par[p] + 1] + 1
		if (Child[Par[p] + 1] <= 2) {
			generateTrees(p + 1, s, cL, n)
		}
		Child[Par[p] + 1] <<- Child[Par[p] + 1] - 1
		while (Par[p] > 0) {
			s = Par[p]
			if (WIENER) { WIndex <<- WIndex + decreases[Sub[s]]; Sub[s] <<- Sub[s] - 2 }
			Par[p] <<- Par[s]
			Child[Par[p] + 1] <<- Child[Par[p] + 1] + 1
			if (Child[Par[p] + 1] <= 2) {
				generateTrees(p + 1, s, p - s, n)
			}
			Child[Par[p] + 1] <<- Child[Par[p] + 1] - 1
		}
	}
}

### This is a utility function used to maintain the correct value of a tree's Wiener index
updateSubVector = function(old, new) {
	cur = old
	while (cur != new && cur > 0) {
		WIndex <<- WIndex + decreases[Sub[cur]]
		Sub[cur] <<- Sub[cur] - 2
		cur = Par[cur]
	}
	if (cur == 0) {
		cur = new
		while (cur > 0) {
			WIndex <<- WIndex + increases[Sub[cur]]
			Sub[cur] <<- Sub[cur] + 2
			cur = Par[cur]
		}
	}
}

### This function processes a given tree, either computing its spectrum or extracting it to save it
processCurTree = function(n) {
  if (!is.null(specialInds)) {
    if (index == nextInd) {
      curTree = extractTree(Par, Child, n)
      specialTrees[[as.character(index)]] <<- curTree
      pos <<- pos + 1
      if (pos <= length(specialInds)) {
        nextInd <<- specialInds[pos]
      }
      if (pos %% 1000 == 0) { print(paste("Processed", pos, "collisions")) }
    }
  }
  else {
    curTree = extractTree(Par, Child, n)
    if (curTree$Wiener %in% allTr) { # if WIENER = TRUE, the index must be eligible 
      curSpectrum = computeSpectrum(curTree, lap = LAP, full = FULL, trace = TRUE, dist = DIST, K = 0)[useInds]
      curSpectrum = round(curSpectrum, ROUND)
      curTrace = as.character(curSpectrum[length(useInds)])
      curIndex = allIndices[curTrace]
      if (!is.na(curIndex)) {
        allSpectra[[curTrace]][curIndex,] <<- c(curSpectrum[-length(useInds)], index)
        allIndices[curTrace] <<- curIndex + 1
      }
    }
  }
	if (index %% 100000 == 0) { print(index) }
	index <<- index + 1
}

### This function extracts a tree from the current state of the specified global variables
### Par is the parent array, Child stores the number of children, N is the number of internal nodes - 1
extractTree = function(Par, Child, N) {
	if (is.null(Child)) {
		Child = rep(0, N)
		for (ind in 1:N) {
			Child[Par[ind] + 1] = Child[Par[ind] + 1] + 1
		}
	}
	ParExt = c(Par, rep(0:N, 2 - Child))
	edges = cbind(ParExt + N + 3, c(N + 3 + (1:N), 1:(N + 2)))
	Tree = list(Nnode = N + 1, tip.label = as.character(1:(N + 2)), edge = edges, edge.length = rep(1, 2 * N + 2), Wiener = WIndex)
	class(Tree) = "phylo"
	Tree
}

### This function computes the adjacency/Laplacian regular/distance matrix of full tree/only leaves
computeMatrix = function(curTree, full, lap, dist, addRow = NULL, addCol = NULL) {
	if (dist) {
		if (full) {
			curMatrix = dist.nodes(curTree)
		}
		else {
			curMatrix = cophenetic.phylo(curTree)
		}
	}
	else {
		Dim = 2 * curTree$Nnode + 1
		curMatrix = matrix(0, Dim, Dim)
		curMatrix[curTree$edge] = 1
		curMatrix = curMatrix + t(curMatrix)
	}
	if (!is.null(addRow)) {
		Dim = nrow(curMatrix)
		heights = ifelse(rep(full, Dim), curMatrix[(Dim + 3)/2, ], node.depth.edgelength(curTree)[1:Dim])
		curMatrix = rbind(curMatrix,   ifelse(rep(addRow, Dim), heights, rep(1, Dim)))
		curMatrix = cbind(curMatrix, c(ifelse(rep(addCol, Dim), heights, rep(1, Dim)), 0))	
	}
	if (lap) {
		curMatrix = diag(rowSums(curMatrix)) - curMatrix
	}
	curMatrix
}

### This function computes the trace of the powers 1:K of a tree matrix based on certain parameters
computeTraces = function(curTree, full, lap, dist, K = 5) {
	curMatrix = computeMatrix(curTree, full = full, lap = lap, dist = dist)
	output = rep(NA, K)
	output[1] = sum(diag(curMatrix))
	for (ind in 2:K) {
		curMatrix = curMatrix %*% curMatrix
		output[ind] = sum(diag(curMatrix))
	}
	output
}

### This function computes the spectrum of a phylo tree based on certain parameter specifications 
computeSpectrum = function(curTree, full, lap, trace, dist, K = 0, addRow = NULL, addCol = NULL) {
	curMatrix = computeMatrix(curTree, full = full, lap = lap, dist = dist, addRow = addRow, addCol = addCol)
	if (K == 0) {
		curSpectrum = eigen(curMatrix, symmetric = TRUE, only.values = TRUE)$values
	}
	else { # Note that this option may be unstable for small matrices!
		curSpectrum = eigs_sym(curMatrix, k = K, which = "BE", opts = list(retvec = FALSE))$values
	}
	if (trace) {
		curSpectrum[length(curSpectrum)] = ifelse(WIENER, WIndex, sum(diag(curMatrix)))
	}
	curSpectrum
}

### This function merges two subtrees into one; it is similar to bind.tree with added root edges
mergeSubtrees = function(tree1, tree2) { 
	L1 = tree1$Nnode + 1
	L2 = tree2$Nnode + 1
	L = L1 + L2
	tree = list(Nnode = L - 1, tip.label = 1:L)
	Edges1 = tree1$edge
	Edges2 = tree2$edge
	Edges1[Edges1 > L1]  = Edges1[Edges1 > L1] + L2 + 1
	Edges2[Edges2 > L2]  = Edges2[Edges2 > L2] + 2 * L1
	Edges2[Edges2 <= L2] = Edges2[Edges2 <= L2] + L1
	root = L + 1
	root1 = ifelse (L1 == 1, 1, root + 1)
	root2 = ifelse (L2 == 1, L, root + L1) 
	tree$edge = rbind(c(root, root1), c(root, root2), Edges1, Edges2)
	class(tree) = "phylo"
	tree
}

### This function splits a list of items based on frequency values into approximately equal parts
splitFrequencyMatrix = function(freqList, numParts = CORES, sizes = FALSE) {
	L = length(freqList)
	freqListS = sort(freqList, decreasing = TRUE)
	M = floor(L / numParts)
	N = L %% numParts
	splitVector = c(rep(1:numParts, M), rep(1:N, (N > 0)))
	Lists = split(freqListS, splitVector)
	if (sizes) {
		partSums = sapply(Lists, sum)
		output = list(Lists, partSums)
	}
	else {
		output = Lists
	}
	output
}

computeDiameters = function(N = 30, offset = 2) {
  
}

### This function computes the diameter of a phylo tree without multifurcations
computeDiameter = function(tree) {
  Tab = computeLRDepths(tree)
  diam = max(rowSums(Tab))
  diam
}

### This function computes the Wiener index of a phylo tree without multifurcations
computeWienerIndex = function(tree, double = DOUBLE) {
  q = rowSums(computeLRSizes(tree)) + 1
  n = length(q)
  N = 2 * n + 1
  stopifnot(q[1] == N)
  W = (1 + double) * (sum(q * (N - q)) + (n + 1) * (N - 1))
  W
}

### This function computes the betweenness centrality of a phylo tree without multifurcations
computeBetweenness = function(tree) {
  Tab = computeLRSizes(tree)
  n = nrow(Tab)
  rSums = rowSums(Tab)
  Centralities = c(rep(0, n + 1), Tab[,1] * Tab[,2] + rSums * (2 * n - rSums))
  Centralities
}

### This function computes the closeness centrality of a phylo tree without multifurcations
computeCloseness = function(tree) {
  return(1/computeFarness(tree))
}

### This function computes the farness of each node of a phylo tree without multifurcations
computeFarness = function(tree) {
  sizes = rowSums(computeLRSizes(tree))
  n = Ntip(tree)
  N = 2 * n - 1
  Farness = rep(NA, N)
  Farness[n + 1] = sum(sizes)
  edges = tree$edge
  for (ind in 1:(N - 1)) {
    curRow = edges[ind,]
    kid = curRow[2]
    Farness[kid] = N + Farness[curRow[1]] - 2 * (1 + ifelse(kid <= n, 0, sizes[kid - n]))
  }
  Farness
}

computeEigenvector = function(tree, fast = TRUE) {
  if (fast) {
    adj_matrix = get.adjacency(as.igraph.phylo(tree, directed = FALSE), sparse = TRUE)
    EV = eigs_sym(adj_matrix, k = 1, which = "LM", opts = list(retvec = TRUE))
  }
  else {
    Edges <<- tree$edge
    N = nrow(Edges) + 1
    findPos = findFirstPositions(Edges[,1])
    first  <<- findPos[[1]]
    second <<- findPos[[2]]
    updateVector = function(x, args) {
      y = rep(0, length(x))
      y[Edges[ ,2]] = x[Edges[ ,1]]
      y[Edges[first, 1]] = y[Edges[first, 1]] + x[Edges[first, 2]]
      y[Edges[second,1]] = y[Edges[second,1]] + x[Edges[second,2]]
      y
    }
    EV = eigs_sym(updateVector, k = 1, which = "LM", n = N, opts = list(retvec = TRUE))
  }
  evector = abs(EV$vectors[,1])
  names(evector) = ifelse (fast, rownames(adj_matrix), 1:N)
  evalue = abs(EV$values[1])
  output = list(evector, evalue)
  output
}

### This function determines the first and last position of occurrence of a "promise" vector,
### where the promise is that it has consecutive elements, each occurs twice and min is first
findFirstPositions = function(dupVector) {
  Min = dupVector[1]
  N = length(dupVector)/2
  Max = N + Min - 1
  reps = rep(FALSE, Max)
  first = rep(NA, N)
  second = rep(NA, N)
  pos1 = 1
  pos2 = 1
  for (ind in 1:(2 * N)) {
    cur = dupVector[ind]
    if (reps[cur]) {
      second[pos2] = ind
      pos2 = pos2 + 1
    }
    else {
      first[pos1] = ind
      reps[cur] = TRUE
      pos1 = pos1 + 1
    }
  }
  output = list(first, second)
  output
}

### This function computes the sizes of the left and right subtrees rooted at internal nodes
computeLRSizes = function(tree) {
  return(computeLRValues(tree, sum))
}

### This function computes the depths of the left and right subtrees rooted at internal nodes
computeLRDepths = function(tree) {
  return(computeLRValues(tree, max))
}

### This function factory recursively computes values for subtrees rooted at internal nodes  
computeLRValues = function(tree, FUN) {
  n = Ntip(tree)
  N = 2 * n - 1
  Tab = matrix(NA, n - 1, 2)
  edges = tree$edge
  for (ind in (N - 1):1) {
    curRow = edges[ind,] - n
    pos = Tab[curRow[1], 1]
    Tab[curRow[1], 2 - is.na(pos)] = 1 + ifelse(curRow[2] <= 0, 0, FUN(Tab[curRow[2],]))
  }
  Tab
}

### This function uses a BFS to compute the depth of each node (internal or leaf) in a tree
computeDepths = function(tree) {
  n = Ntip(tree)
  myGraph = as.igraph(tree)
  depths = bfs(myGraph, root = 1, neimode = "out", dist = TRUE)$dist
  depths = c(tail(depths, n), head(depths, n - 1))
  depths
}

### This function computes the I2 statistic of a tree
computeI2 = function(tree) {
  n = Ntip(tree)
  LRMat = computeLRSizes(tree)
  values = abs(LRMat[,1] - LRMat[,2]) / (rowSums(LRMat) - 2)
  output = sum(values[is.finite(values)])
  output
}

### This function computes the B1 statistic of a tree
computeB1 = function(tree) {
  n = Ntip(tree)
  depths = node.depth(tree, method = 2)
  output = sum(1/(depths[-(1:(n+1))] - 1))
  output
}

### This function computes the B2 statistic of a tree
computeB2 = function(tree) {
  n = Ntip(tree)
  depths = node.depth.edgelength(tree)[1:n]
  output = sum(depths/2^depths)
  output
}