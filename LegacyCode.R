### From SpectralComparison.R
### When starting a new instance of SpectralComparison.R do this!
install.packages(c("combinat","rARPACK","ape","igraph","phangorn"))

checkExchangeProperty = function(Tr1, Tr2) {
  Maxes = rep(NA, 8)
  ind = 1
  for (full in c(FALSE, TRUE)) {
    for (addRow in c(FALSE, TRUE)) {
      for (addCol in c(FALSE, TRUE)) {
        Spec1 = computeSpectrum(Tr1, full = full, lap = FALSE, trace = FALSE, dist = TRUE, addRow = addRow, addCol = addCol)
        Spec2 = computeSpectrum(Tr2, full = full, lap = FALSE, trace = FALSE, dist = TRUE, addRow = addRow, addCol = addCol)
        Maxes[ind] = max(abs(Spec1 - Spec2))
        ind = ind + 1
      }
    }
  }
  Maxes
}

postDriver = function(Vec, size, tol = 1e-10) {
  # Postprocessng for LAP = FALSE, DIST = TRUE, FULL = FALSE
  # Confirms that leaf distance cospectrality <=> node distance cospectrality <=> exchange property
  # Verified for size = 17, 18, 19, 20 on all 1, 2, 4, 9 cospectral pairs, respectively
  goodInds = which(Vec < tol)
  load(paste0("SpecialPairs", size, ".Rdata"))
  goodPairs = specialPairs[,goodInds]
  for (ind in 1:length(goodInds)) {
    print(ind)
    curPair = goodPairs[,ind]
    Tr1 = specialTrees[[as.character(curPair[1])]]
    Tr2 = specialTrees[[as.character(curPair[2])]]
    MAX = max(checkExchangeProperty(Tr1, Tr2)) 
    print(MAX)
  }
}

### Wiener index computation
WI = 0
if (Wiener) {
	WI = computeWienerIndex(ParExt)
}
	
computeWienerIndex = function(ParExt) {
	N = length(ParExt) + 1
	q = rep(0, N)
	for (ind in N:1) {
		q[ind] = q[ind] + 1
		curPar = ParExt[ind - 1] + 1
		q[curPar] = q[curPar] + q[ind]
	}
	stopifnot(q[1] == N)
	stopifnot(all(q[1 + (1:length(Sub))] == Sub)) 
	W = (1 + DOUBLE) * sum(q * (N - q))
	W
}

generateSpectraOptimize = function(maxSize = SIZE, full = FULL) {
	### An attempt at optimizing trace generation using tricks ###
	baseFilename = "Traces"
	if (full) {
		baseFilename = paste0(baseFilename, "AllNodes")
	}
	else {
		baseFilename = paste0(baseFilename, "LeavesOnly")
	}
	allHeightSums = vector("list", length = maxSize - 1)
	allTraces   = vector("list", length = maxSize)
	allSizes = rep(NA, maxSize)
	allHeightSums[[1]] = 0
	allHeightSums[[2]] = 2
	allTraces[[1]] = 0
	if (full) {
		allTraces[[2]] = 8
	}
	else {
		allTraces[[2]] = 4
	}
	allSizes[1:2] = 1
	for (ind in 3:maxSize) {
		halfSize = ceiling((ind + 1)/2)
		range = halfSize:(ind - 1)
		allSizes[ind] = sum(allSizes[range] * allSizes[ind - range])
		if (ind %% 2 == 0) {
			allSizes[ind] = allSizes[ind] + choose(allSizes[ind/2] + 1, 2)
		}
	}
	for (ind in 3:maxSize) {
		print(paste("Processing the", allSizes[ind], "trees of size", ind))
		if (full) {
			addHeight = 2 * (ind - 1)
		}
		else {
			addHeight = ind
		}
		halfSize = ceiling((ind + 1)/2)
		range = halfSize:(ind - 1)
		if (ind %% 2 == 0) {
			extraRange = 1:(choose(allSizes[ind/2] + 1, 2))
			curStart = tail(extraRange, 1)
			extraInds = combn2(1:(allSizes[ind/2] + 1))
			extraInds = extraInds[order(extraInds[,2]), 2:1, drop = FALSE]
			extraInds[,1] = extraInds[,1] - 1
		}
		else {
			curStart = 0
		}
		if (ind < maxSize) {
			curHeights = rep(NA, allSizes[ind])	
			if (ind %% 2 == 0) {
				extraHeights = outer(allHeightSums[[ind/2]], allHeightSums[[ind/2]], "+")
				extraHeights = extraHeights[extraInds]
				curHeights[extraRange] = extraHeights
			}
			for (indexL in range) {
				indexR = ind - indexL
				curSize = allSizes[indexL] * allSizes[indexR]
				curHeights[curStart + (1:curSize)] = outer(allHeightSums[[indexR]], allHeightSums[[indexL]], "+")
				curStart = curStart + curSize
			}
			allHeightSums[[ind]] = curHeights + addHeight
		}
		curTraces = rep(NA, allSizes[ind])
		if (ind %% 2 == 0) {
			if (full) {
				extraSummand = 2 * (allHeightSums[[ind/2]] + ind - 1) * ind
			}
			else {
				extraSummand = (allHeightSums[[ind/2]] + ind/2) * ind
			}
			extraTraces = outer(extraSummand, extraSummand, "+") + outer(allTraces[[ind/2]], allTraces[[ind/2]], "+")
			curTraces[extraRange] = extraTraces[extraInds]
			curStart = tail(extraRange, 1)
		}
		else {
			curStart = 0
		}
		for (indexL in range) {
			print(indexL)
			indexR = ind - indexL
			curSize = allSizes[indexL] * allSizes[indexR]
			if (full) {
				summandR = 4 * (allHeightSums[[indexR]] + 2 * indexR - 1) * indexL
				summandL = 4 * (allHeightSums[[indexL]] + 2 * indexL - 1) * indexR
			}
			else {
				summandR = 2 * (allHeightSums[[indexR]] + indexR) * indexL
				summandL = 2 * (allHeightSums[[indexL]] + indexL) * indexR
			}
			curTraces[curStart + (1:curSize)] = outer(summandR, summandL, "+") + outer(allTraces[[indexR]], allTraces[[indexL]], "+")
			curStart = curStart + curSize
		}
		if (ind < maxSize) {
			allTraces[[ind]] = curTraces
		}
		save(curTraces, file = paste0(baseFilename, ind, ".Rdata"))
	}
	# output = list(allSizes, allHeightSums, allTraces)
	# output
}

generateSpectra = function(maxSize = 20, full = TRUE, lap = TRUE) {
	allMatrices = vector("list", length = maxSize - 1)
	allHeights  = vector("list", length = maxSize - 1)
	allSizes = rep(NA, maxSize)
	allDists = rep(NA, maxSize)
	allMatrices[[1]] = list(matrix(0))
	if (full) {
		allMatrices[[2]] = list(rbind(c(0, 1, 1), c(1, 0, 2), c(1, 2, 0)))
	}
	else {
		allMatrices[[2]] = list(rbind(c(0, 2), c(2, 0)))
	}
	allHeights[[1]]  = matrix(0)
	if (full) {
		allHeights[[2]]  = matrix(c(0, 1, 1), nrow = 1)
	}
	else {
		allHeights[[2]] = matrix(c(1, 1), nrow = 1)
	}
	allSizes[1:2] = 1
	allDists[1:2] = 0
	for (ind in 3:maxSize) {
		halfSize = ceiling((ind + 1)/2)
		range = halfSize:(ind - 1)
		allSizes[ind] = sum(allSizes[range] * allSizes[ind - range])
		if (ind %% 2 == 0) {
			allSizes[ind] = allSizes[ind] + choose(allSizes[ind/2] + 1, 2)
			range = c(ind/2, range)
		}
		curSize = allSizes[ind]
		print(paste("Processing the", curSize, "trees of size", ind))
		if (full) {
			numNodes = 2 * ind - 1
		}
		else {
			numNodes = ind
		}
		if (ind < maxSize) {
			curMatrices = vector("list", length = curSize)
			curHeights  = matrix(0, nrow = curSize, ncol = numNodes)
		}
		curSpectra  = matrix(0, nrow = curSize, ncol = numNodes)
		curMatrix = matrix(0, nrow = numNodes, ncol = numNodes)
		curHeight = rep(0, numNodes)
		fullIndex = 1
		for (left in range) {
			print(paste("Processing the left subtrees of size", left))
			right = ind - left
			if (full) {
				subRange1 = 2:(2 * left)
				subRange2 = (2 * left + 1):numNodes
			}
			else {
				subRange1 = 1:left
				subRange2 = (left + 1):numNodes
			}
			for (indL in 1:allSizes[left]) {
				curMatL = allMatrices[[left]][[indL]]
				curHeightL = allHeights[[left]][indL,]
				curHeight[subRange1] = curHeightL + 1
				curMatrix[subRange1, subRange1] = curMatL
				upperBound = allSizes[right]
				if (ind %% 2 == 0 && left == right) {
					upperBound = indL # forces the left side to be no smaller than the right!
				}
				for (indR in 1:upperBound) {
					curMatR = allMatrices[[right]][[indR]]
					curHeightR = allHeights[[right]][indR,]
					curHeight[subRange2] = curHeightR + 1
					curOuter = 2 + outer(curHeightL, curHeightR, "+")
					if (full) {
						curMatrix[1,] = curHeight
						curMatrix[,1] = curHeight
					}
					curMatrix[subRange1, subRange2] = curOuter
					curMatrix[subRange2, subRange1] = t(curOuter)
					curMatrix[subRange2, subRange2] = curMatR
					if (lap) {
						curSpectrum = eigen(diag(rowSums(curMatrix)) - curMatrix, symmetric = TRUE, only.values = TRUE)$values
					}
					else {
						curSpectrum = eigen(curMatrix, symmetric = TRUE, only.values = TRUE)$values
					}
					if (ind < maxSize) {
						curMatrices[[fullIndex]] = curMatrix
						curHeights[fullIndex,] = curHeight
					}
					curSpectra[fullIndex,] = curSpectrum
					fullIndex = fullIndex + 1
				}
			}
		}
		print("Computing pairwise distances")
		curDist = findClosePairs(curSpectra)
		if (length(curDist) > 0) {
			uniqueClose = unique(curDist[-3,])
			output = list(best = curDist[3,], bestSpectra = curSpectra[uniqueClose,])
			if (ind < maxSize) {
				output = c(output, list(bestMatrices = curMatrices[uniqueClose]))
			}
			return(output)
		}
		if (ind < maxSize) {
			allMatrices[[ind]] = curMatrices
			allHeights[[ind]] = curHeights
		}
		allDists[ind] = min(c(0, curDist[3,]))
	}
	output = list(allMatrices, allHeights, curSpectra, allDists)
}

generateSpectraPhylo = function(maxSize = 20, full = TRUE, lap = TRUE) {
	allTrees = vector("list", length = maxSize)
	allTrees[[1]] = list(Singleton)
	allTrees[[2]] = list(Cherry)
	for (ind in 3:maxSize) {
		halfSize = ceiling((ind + 1)/2)
		range = halfSize:(ind - 1)
		if (ind %% 2 == 0) {
			range = c(ind/2, range)
		}
		curSize = allSizes[ind]
		curTrees = vector("list", length = curSize)
		print(paste("Processing the", curSize, "trees of size", ind))
		if (full) {
			numNodes = 2 * ind - 1
		}
		else {
			numNodes = ind
		}
		curSpectra  = matrix(0, nrow = curSize, ncol = numNodes)
		fullIndex = 1
		for (left in range) {
			print(paste("Processing the left subtrees of size", left))
			right = ind - left
			for (indL in 1:allSizes[left]) {
				curTreeL = allTrees[[left]][[indL]]
				upperBound = allSizes[right]
				if (ind %% 2 == 0 && left == right) {
					upperBound = indL # forces the left side to be no smaller than the right!
				}
				for (indR in 1:upperBound) {
					curTreeR = allTrees[[right]][[indR]]
					curTree = mergeSubtrees(curTreeL, curTreeR)
					curSpectrum = computeSpectrum(curTree, full = full, lap = lap)
					curSpectra[fullIndex,] = curSpectrum
					curTrees[[fullIndex]] = curTree
					fullIndex = fullIndex + 1
				}
			}
		}
		print("Computing pairwise distances")
		curDist = findClosePairs(curSpectra)
		if (length(curDist) > 0) {
			uniqueClose = unique(curDist[-3,])
			output = list(best = curDist[3,], bestSpectra = curSpectra[uniqueClose,], bestPairs = curDist[-3,], bestTrees = curTrees[uniqueClose])
			return(output)
		}
		allTrees[[ind]] = curTrees
	}
	allTrees
}

reconstructTree = function(sequence, hash = NULL, minSize = NULL, sep = ",") {
	L = length(sequence)
	if (!is.null(minSize) && L <= minSize) {
		hashKey = paste(sequence, sep = sep)
		return(get(hashKey, hash))
	}
	if (L == 1) {
		return(Cherry)
	}
	else if (L == 2) {
		return(Triplet)
	}
	leftLength = max(which(sequence - 1:L == 0))
	tree1 = reconstructTree(sequence[1:(leftLength - 1)], hash = hash, minSize = minSize, sep = sep)
	if (leftLength < L) {
		tree2 = reconstructTree(sequence[(leftLength + 1):L], hash = hash, minSize = minSize, sep = sep)
	}
	else {
		tree2 = Singleton
	}
	return(mergeSubtrees(tree1, tree2))
}

findClosePairs = function(Matrix, method = "euclidean", tol = 1e-4) {
	roundedMatrix = round(Matrix)
	allGroups = extractUniqueRows(roundedMatrix, repeats = TRUE)[[2]]
	allLens = sapply(allGroups, length)
	goodGroups = allGroups[allLens > 1]
	goodLens = allLens[allLens > 1]
	allClose = matrix(NA, nrow = 3, ncol = 0)
	if (length(goodGroups) > 0) {
		for (ind in 1:length(goodGroups)) {
			curGroup = goodGroups[[ind]]
			curLen = allLens[ind]
			curDists = dist(Matrix[curGroup,], method = method)
			curClose = which(curDists < tol)
			if (length(curClose) > 0) {
				bestDists = curDists[curClose]
				curClose = convertSinglesToPairs(curClose, curLen)
				allClose = cbind(allClose, rbind(curGroup[curClose[1,]], curGroup[curClose[2,]], bestDists))	
			}
		}
	}
	output = allClose
}

convertSinglesToPairs = function(indices, Dim) {
	### Translates the positions in the vectorized strict lower triangle in a matrix
	### into the positions in the actual matrix; useful for, eg, the output of dist.
	L = length(indices)
	output = matrix(0, nrow = 2, ncol = L)
	boundaries = c(0, cumsum((Dim - 1):1))
	for (ind in 1:Dim) {
		curBoundary = boundaries[ind]
		curExceeds  = (indices > curBoundary)
		output[1, curExceeds] = output[1, curExceeds] + 1 
	}
	output[2,] = indices - boundaries[output[1,]]
	output[2,] = output[2,] + output[1,]
	output 
}

postprocessHeights = function(DMatrix) {
  Vars = apply(curHeights, 1, var)
  N = length(Vars)
  stopifnot(nrow(DMatrix) == N)
  stopifnot(ncol(DMatrix) == N)
  H = diag(rep(1, N)) - 1/N
  SMatrix = -(1/2) * t(H) %*% DMatrix^2 %*% H
  Range = range(eigen(SMatrix)$values)
  centVars = Vars - mean(Vars)
  normVars = centVars/sqrt(sum(centVars^2))
  raw = - (1/2) * t(normVars) %*% DMatrix^2 %*% normVars
  value = (raw - Range[1])/(Range[2] - Range[1])
  print(value)
  value
}

### alternative to circular = TRUE in splitFrequencyMatrix
if (circular) {
  
} else {
  partSum = ceiling(Sum/numParts)
  breaks = rep(NA, numParts)
  breaks[numParts] = L
  pos = 1
  ind = 1
  curSum = 0
  while (ind < numParts) {
    curSum = curSum + freqListS[pos]
    if (curSum >= partSum) {
      breaks[ind] = pos
      ind = ind + 1
      curSum = 0
    }
    pos = pos + 1
  }
  Lists = split(freqListS, rep(1:numParts, diff(c(0, breaks)))) 
}

convertToLists = function(Pairs) {
  heads = setdiff(Pairs[1,], Pairs[2,])
  goodPairs = Pairs[, Pairs[1,] %in% heads]
  goodPairs = goodPairs[, order(goodPairs[1,])]
  splitPairs = split(goodPairs[2,], goodPairs[1,])
  splitPairs
}

computeMaxSize = function(Lists, tail = FALSE) { 
  # upper bound on hash size needed to process all the lists
  # if tail = TRUE, each list is headed by its LAST element!
  FACTOR = ifelse(tail, -1, 1)
  L = length(Lists)
  heads = as.numeric(names(Lists))
  if (!tail) {
    tails = sapply(Lists, function(x) {tail(x, 1)})
  }
  else {
    tails = sapply(Lists, function(x) {x[1]})
  }
  order = order(tails)
  tails = tails[order]
  sizesH = sapply(Lists, length)
  names(sizesH) = NULL
  sizesT = sizesH[order]
  count = 1
  maxCount = 1
  # Perform a pointer walk
  ind1 = 1
  ind2 = 1
  while ((!tail && ind1 <= L) || (tail && ind2 <= L)) {
    if (heads[ind1] < tails[ind2]) { # "opening" event
      count = count + FACTOR * sizesH[ind1]
      if (!tail && count > maxCount) {
        maxCount = count
      }
      ind1 = ind1 + 1
    }
    else { # "closing" event
      count = count - FACTOR * sizesT[ind2]
      if (tail && count > maxCount) {
        maxCount = count
      }
      ind2 = ind2 + 1
    }
  }
  maxCount
}

### from processCurTree
if (height) {
  curTree = extractTree(Par, Child, n)
  curHeights[index, ] <<- node.depth.edgelength(curTree)[1:(n+2)]
}

### from driver
if (HASH) {
  postprocessSpecialPairs(size = size)
} else {
  
}
if (TREE) {
  
} else {
  specialLists = convertToLists(specialPairs[2:1, ])
  maxSize = computeMaxSize(specialLists, tail = TRUE)
  specialTails = as.numeric(names(specialLists))
  posT = 1
  nextTail = specialTails[posT]
  spectralHash = new.env(hash = TRUE, size = maxSize)
  Dists = vector("list", length(specialLists))
  names(Dists) = specialTails
  print(paste("There are", length(specialTails), "groups to process"))
  generateSpectraLinspace(size = size)
  minDists = sapply(Dists, min)
  bestDist = which.min(minDists)
  print(paste("The closest distance is", minDists[bestDist]))
  write.table(as.matrix(Dists[[bestDist]]), row.names = FALSE, col.names = FALSE, sep = ",", file = paste0("BestDists", size, "X", num, ".csv"))
  bestTail = specialTails[bestDist]
  bestInds = c(specialLists[[as.character(bestTail)]], bestTail)
  save(bestInds, file = paste0("BestInds", size, "X", num, ".Rdata"))
}

### from generateSpectraLinspace
if (TRACES) {
  useInds <<- 1:K
}	else {

}

if (HASH) {
  allSpectra = new.env(hash = TRUE, size = numRows)
} else {

}

if (!TREE) {
  dimension <<- ifelse(FULL, 2 * SIZE - 1, SIZE)
  specialTrees <<- matrix(NA, length(specialInds), dimension, dimnames = list(specialInds, NULL))
}

if (HEIGHT) {
  curHeights <<- matrix(NA, allSizes[size], size)
}

if (SPLIT)	{

} else {
  allSpectra = matrix(0, numRows, numCols)
}

generateSpectraOptimize = function(maxSize = SIZE, full = FULL) {
  ### A new attempt at optimizing trace generation using tricks ###
  baseFilename = "Traces"
  if (full) {
    baseFilename = paste0(baseFilename, "AllNodes")
  }
  else {
    baseFilename = paste0(baseFilename, "LeavesOnly")
  }
  allTraces   = vector("list", length = maxSize)
  allTraces[[1]] = cbind(0, 0, 1)
  allTraces[[2]] = cbind(2, ifelse(full, 8, 4), 1)
  for (ind in 3:maxSize) {
    print(paste("Processing the", allSizes[ind], "trees of size", ind))
    if (full) {
      addHeight = 2 * (ind - 1)
    }
    else {
      addHeight = ind
    }
    halfSize = ceiling((ind + 1)/2)
    range = halfSize:(ind - 1)
    curTracesFull = matrix(NA, 0, 3)
    if (ind %% 2 == 0) {
      curHalf = allTraces[[ind/2]]
      halfHeights = curHalf[, 1]
      halfHeightsMat = outer(halfHeights, halfHeights, "+") + addHeight
      halfTraces = curHalf[, 2]
      if (full) {
        extraSummand = 2 * (halfHeights + ind - 1) * ind
      }
      else {
        extraSummand = (halfHeights + ind/2) * ind
      }
      halfTracesMat = outer(extraSummand, extraSummand, "+") + outer(halfTraces, halfTraces, "+")
      halfFreqs = curHalf[, 3]
      halfFreqsMat = outer(halfFreqs, halfFreqs)
      diag(halfFreqsMat) = diag(halfFreqsMat) + halfFreqs # to ensure that n goes to choose(n + 1, 2)
      halfTriplets = cbind(as.vector(halfHeightsMat), as.vector(halfTracesMat), as.vector(halfFreqsMat))
      curTracesFull = groupTriplets(halfTriplets, FALSE)
      curTracesFull[, 3] = curTracesFull[, 3] / 2
    }
    for (indexL in range) {
      indexR = ind - indexL
      curL = allTraces[[indexL]]
      curR = allTraces[[indexR]]
      curHeightsL = curL[, 1]
      curHeightsR = curR[, 1]
      curHeightsMat = outer(curHeightsR, curHeightsL, "+") + addHeight
      curTracesL = curL[, 2]
      curTracesR = curR[, 2]
      if (full) {
        summandR = 4 * (curHeightsR + 2 * indexR - 1) * indexL
        summandL = 4 * (curHeightsL + 2 * indexL - 1) * indexR
      }
      else {
        summandR = 2 * (curHeightsR + indexR) * indexL
        summandL = 2 * (curHeightsL + indexL) * indexR
      }
      curTracesMat = outer(summandR, summandL, "+") + outer(curTracesR, curTracesL, "+")
      curFreqsL = curL[, 3]
      curFreqsR = curR[, 3]
      curFreqsMat = outer(curFreqsR, curFreqsL)
      curTriplets = cbind(as.vector(curHeightsMat), as.vector(curTracesMat), as.vector(curFreqsMat))
      curTracesFull = groupTriplets(rbind(curTracesFull, curTriplets), FALSE)
    }
    if (ind < maxSize) {
      allTraces[[ind]] = curTracesFull 
    }
    curTracesRed = groupTriplets(curTracesFull, TRUE)
    save(curTracesRed, file = paste0(baseFilename, ind, ".Rdata"))
  }
  allTraces
}

groupTriplets = function(triplets, ignoreFirst = FALSE) {
  if (ignoreFirst) {
    myList = split(triplets[,3], triplets[,2])
    redList = sapply(myList, sum)
  }
  else {
    myList = split(triplets[,3], paste(triplets[,1], triplets[,2]))
    redList = sapply(myList, sum)
    Names = t(matrix(as.numeric(unlist(strsplit(names(redList), " "))), nrow = 2))
    names(redList) = NULL
    redList = cbind(Names, redList, deparse.level = 0)
  }
  redList
}

### from postprocessSpecialTrees
if (TREE) {
  
} else {
  allDists[ind] = sum(abs(specialTrees[as.character(ind1), ] - specialTrees[as.character(ind2), ]))
}

### from processCurTree
if (SAVE) {
  save(curTree, file = paste0("Tree", index, ".RData"))
}

if (TREE) {
  
} else {
  curSpectrum = computeSpectrum(curTree, lap = LAP, full = FULL, trace = SPLIT, dist = DIST, K = 0)
  assign(as.character(index), curSpectrum, envir = spectralHash)
  if (index == nextTail) {
    curList = c(specialLists[[as.character(index)]], index)
    toCompare = mget(as.character(curList), spectralHash)
    toCompare = matrix(unlist(toCompare), nrow = length(curList), byrow = TRUE)
    Dists[[posT]] <<- dist(toCompare)
    rm(list = as.character(curList), envir = spectralHash)
    posT <<- posT + 1
    if (posT <= length(specialTails)) {
      nextTail <<- specialTails[posT]
    }
  }
}

if (TRACES) {
  curSpectrum = computeTraces(curTree, lap = LAP, full = FULL, dist = DIST, K = K)
} else {

}

if (TRACES) {
  curTrace = curSpectrum[1]
  curVector = curSpectrum[-1]
} else {
  
}

if (SPLIT) {

} else {
  allSpectra[index,] = curSpectrum
}

if (HASH) {
  curKey = paste(curSpectrum, collapse = ",")
  if (exists(curKey, envir = allSpectra)) {
    index1 = get(curKey, envir = allSpectra)
    specialPairs = cbind(specialPairs, c(index, index1))
  }
  else {
    assign(curKey, index, envir = allSpectra)
  }
}

### from postprocessAllSpectra
if (SPLIT) {

} else {
  extraMat = findMatchingRows(allSpectra, byLastCol = FALSE)
  specialPairs = cbind(specialPairs, t(extraMat))
}

Q0 = generateSpectraOptimize(full = FALSE, maxSize = 32)
save(Q0, file = "TracesLeavesOnly.Rdata")
Q1 = generateSpectraOptimize(full = TRUE, maxSize = 32)
save(Q1, file = "TracesAllNodes.Rdata")

### From SpectralTrees.R
metaDriver = function(numIter = 100, numNodes = 25, scale = 100, frac = 0.5, K = 5, sep = ",") {
  Sample = vector("list", numIter)
  for (ind in 1:numIter) {
    curTree = rtree(numNodes)
    curTree$tip.label = paste("t", 1:numNodes, sep = "")
    Sample[[ind]] = curTree
  }
  Trees = lapply(Sample, function(x) {x$edge})
  initW = lapply(Sample, function(x) {Lens = x$edge.length; Lens/sum(Lens)*length(Lens)})
  scaledUpW = lapply(Sample, function(x) {Lens = x$edge.length; y = Lens/sum(Lens); Tips = which(x$edge[,2] <= numNodes); 
  y[Tips] = y[Tips] * scale; y/sum(y)*length(y)})	
  someTipW  = lapply(Sample, function(x) {Lens = x$edge.length; y = Lens/sum(Lens); Tips = which(x$edge[,2] <= numNodes);
  goodTips = sample(Tips, size = round(frac * length(Tips))); y[goodTips] = y[goodTips] * scale; y/sum(y)*length(y)})
  someIntW = lapply(Sample, function(x) {Lens = x$edge.length; y = Lens/sum(Lens) * scale; Ints = which(x$edge[,2] > numNodes);
  goodInts = sample(Ints, size = round(frac * length(Ints))); y[goodInts] = y[goodInts] / scale; y/sum(y)*length(y)})
  fullW = list(unscaled = initW, scaledUp = scaledUpW, someTips = someTipW, someInts = someIntW)
  fullResults = c()
  count = 0
  weighted = TRUE
  for (directed in c(FALSE, TRUE)) {
    for (Laplacian in c(FALSE, TRUE)) {
      for (norm in c(FALSE, TRUE)) {
        for (mode in c("in", "out")) {
          if ((!Laplacian && norm) || (!directed && mode == "out")) {
            count = count + 1	
            print(count)
          }
          else {
            cur = driver(Trees = Trees, Weights = fullW, K = K, directed = directed, Laplacian = Laplacian, norm = norm, mode = mode)
            curTitle = paste(c(directed, Laplacian, norm, weighted, mode), collapse = sep)
            print(curTitle)
            fullResults[[curTitle]] = cur
          }
        }
      }
    }
  }
  pValues = performTests(fullResults, firsts = c(1, 1, 1, 2, 2, 3), seconds = c(2, 3, 4, 3, 4, 4), sep = sep)
  output = list(results = fullResults, pValues = pValues)
  output
}

driver = function(Trees, Weights, K, directed, Laplacian, norm, mode) {
  Results = list()
  count = 1
  for (name in names(Weights)) {
    print(paste("Processing dataset number", count))
    curW = Weights[[name]]
    curRes = processItem(Trees, K = K, directed = directed, Laplacian = Laplacian, norm = norm, weights = curW, mode = mode)
    Results[[name]] = curRes
    count = count + 1
  }
  Results
}

### From ArtPaperScript.R
### IMPLEMENTING CAROLINE'S EXAMPLE ###
TreeA = mergeSubtrees(mergeSubtrees(mergeSubtrees(Cherry, Singleton), Singleton), Singleton)
TreeB = mergeSubtrees(mergeSubtrees(Cherry, Cherry), Singleton)
TreeC = mergeSubtrees(mergeSubtrees(Cherry, Singleton), Cherry)
TreeD = TreeC
Tree1 = mergeSubtrees(mergeSubtrees(TreeA, TreeB), mergeSubtrees(TreeC, TreeD))
Tree1$edge.length=rep(1,2*(Tree1$Nnode)+1)
Tree2 = mergeSubtrees(mergeSubtrees(TreeA, TreeC), mergeSubtrees(TreeB, TreeD))
Tree2$edge.length=rep(1,2*(Tree2$Nnode)+1)
H1 = node.depth.edgelength(Tree1)
H2 = node.depth.edgelength(Tree2)
print(all(sort(H1) == sort(H2)))
TreeA1 = mergeSubtrees(TreeA, Singleton)
TreeB1 = mergeSubtrees(TreeB, Singleton)
TreeC1 = mergeSubtrees(TreeC, Singleton)
TreeD1 = mergeSubtrees(TreeD, Singleton)
Tree3 = mergeSubtrees(mergeSubtrees(TreeA1, TreeB1), mergeSubtrees(TreeC1, TreeD1))
Tree3$edge.length=rep(1,2*(Tree3$Nnode)+1)
Tree4 = mergeSubtrees(mergeSubtrees(TreeA1, TreeC1), mergeSubtrees(TreeB1, TreeD1))
Tree4$edge.length=rep(1,2*(Tree4$Nnode)+1)
H3 = node.depth.edgelength(Tree3)
H4 = node.depth.edgelength(Tree4)
print(all(sort(H3) == sort(H4)))
Counter1 = mergeSubtrees(mergeSubtrees(Tree1, mergeSubtrees(Tree3, Singleton)), mergeSubtrees(Tree2, mergeSubtrees(Tree4, Singleton)))
Counter1$edge.length=rep(1,2*(Counter1$Nnode)+1)
Counter2 = mergeSubtrees(mergeSubtrees(Tree2, mergeSubtrees(Tree3, Singleton)), mergeSubtrees(Tree1, mergeSubtrees(Tree4, Singleton)))
Counter2$edge.length=rep(1,2*(Counter2$Nnode)+1)
DMat1 = computeMatrix(Counter1, full = TRUE, lap = FALSE, dist = TRUE)
DMat2 = computeMatrix(Counter2, full = TRUE, lap = FALSE, dist = TRUE)
print(all(sort(DMat1) == sort(DMat2)))
DMat1S = t(apply(DMat1,1,sort))
DMat2S = t(apply(DMat2,1,sort))
order1 = do.call(order,as.data.frame(DMat1S))
order2 = do.call(order,as.data.frame(DMat2S))
print(all(DMat1S[order1,] == DMat2S[order2,]))
Spec1 = eigen(DMat1, symmetric = TRUE, only.values = TRUE)$values
Spec2 = eigen(DMat2, symmetric = TRUE, only.values = TRUE)$values
print(max(abs(Spec1 - Spec2)))
LMat1 = computeMatrix(Counter1, full = TRUE, lap = TRUE, dist = TRUE)
LMat2 = computeMatrix(Counter2, full = TRUE, lap = TRUE, dist = TRUE)
print(all(sort(LMat1) == sort(LMat2)))
LMat1S = t(apply(LMat1,1,sort))
LMat2S = t(apply(LMat2,1,sort))
Order1 = do.call(order,as.data.frame(LMat1S))
Order2 = do.call(order,as.data.frame(LMat2S))
print(all(LMat1S[Order1,] == LMat2S[Order2,]))
Spec3 = eigen(LMat1, symmetric = TRUE, only.values = TRUE)$values
Spec4 = eigen(LMat2, symmetric = TRUE, only.values = TRUE)$values
print(max(abs(Spec3 - Spec4)))

# ####################################################################
# ### WORKING ON THE CO-SPECTRALITY PROBLEM FOR LAPLACIAN DISTANCE ###
# ####################################################################
TrA = createEmptyTree(numNodes = 15)
TrA=createNewNode(TrA, c(2,3), c(1,1), 9);  	TrA=createNewNode(TrA, c(5,6), c(1,1), 10);   	TrA=createNewNode(TrA, c(7,8),  c(1,1), 11)
TrA=createNewNode(TrA, c(1,9), c(1,1), 12); 	TrA=createNewNode(TrA, c(10,11), c(1,1), 13); 	TrA=createNewNode(TrA, c(4,13), c(1,1), 14); 	TrA=createNewNode(TrA, c(12,14), c(1,1), 15);
TreeA=convertPhylo(TrA)
TreeAG=as.igraph(TreeA)
TrB = createEmptyTree(numNodes = 15)
TrB=createNewNode(TrB, c(1,2), c(1,1), 9);  	TrB=createNewNode(TrB, c(3,4), c(1,1), 10);   	TrB=createNewNode(TrB, c(5,6),  c(1,1), 11)
TrB=createNewNode(TrB, c(7,8), c(1,1), 12); 	TrB=createNewNode(TrB, c(9,10), c(1,1), 13); 	TrB=createNewNode(TrB, c(11,12), c(1,1), 14); 	TrB=createNewNode(TrB, c(13,14), c(1,1), 15);
TreeB=convertPhylo(TrB)
TreeBG=as.igraph(TreeB)
AdjA = get.adjacency(TreeAG, sparse = FALSE)
AdjB = get.adjacency(TreeBG, sparse = FALSE)
AdjA = AdjA + t(AdjA)
AdjB = AdjB + t(AdjB)
vec1 = eigen(AdjA)$values
vec2 = eigen(AdjB)$values
all.equal.numeric(vec1, vec2)
LapA = diag(rowSums(AdjA)) - AdjA
LapB = diag(rowSums(AdjB)) - AdjB
Vec1 = eigen(LapA)$values
Vec2 = eigen(LapB)$values
all.equal.numeric(Vec1, Vec2)
allPairsA = distances(TreeAG, mode = "all")
allPairsB = distances(TreeBG, mode = "all")
modLapA = diag(rowSums(allPairsA)) - allPairsA
modLapB = diag(rowSums(allPairsB)) - allPairsB
modVec1 = eigen(modLapA)$values
modVec2 = eigen(modLapB)$values
require("RPANDA")
altVec1 = spectR(TreeA)$eigenvalues
altVec2 = spectR(TreeB)$eigenvalues
all.equal.numeric(modVec1, altVec1)
all.equal.numeric(modVec2, altVec2)
TrA = createEmptyTree(numNodes = 33)
TrA=createNewNode(TrA, c(2,3), c(1,1), 18);		TrA=createNewNode(TrA, c(5,6), c(1,1), 19);		TrA=createNewNode(TrA, c(10,11),  c(1,1), 20);		TrA=createNewNode(TrA, c(14,15), c(1,1), 21); 
TrA=createNewNode(TrA, c(16,17), c(1,1), 22); 	TrA=createNewNode(TrA, c(18,1), c(1,1), 23);	TrA=createNewNode(TrA, c(4,19), c(1,1), 24);   		TrA=createNewNode(TrA, c(9,20), c(1,1), 25);
TrA=createNewNode(TrA, c(21,22), c(1,1), 26); 	TrA=createNewNode(TrA, c(23,24), c(1,1), 27);   TrA=createNewNode(TrA, c(8,25),  c(1,1), 28); 		TrA=createNewNode(TrA, c(13,26), c(1,1), 29); 
TrA=createNewNode(TrA, c(12,29), c(1,1), 30); 	TrA=createNewNode(TrA, c(28,30), c(1,1), 31); 	TrA=createNewNode(TrA, c(7,31), c(1,1), 32); 		TrA=createNewNode(TrA, c(27,32), c(1,1), 33);
TreeA=convertPhylo(TrA)
TreeAG=as.igraph(TreeA)
TrB = createEmptyTree(numNodes = 33)
TrB=createNewNode(TrB, c(1,2), c(1,1), 18);  	TrB=createNewNode(TrB, c(3,4), c(1,1), 19);   	TrB=createNewNode(TrB, c(6,7),  c(1,1), 20); 	TrB=createNewNode(TrB, c(12,13), c(1,1), 21);
TrB=createNewNode(TrB, c(15,16), c(1,1), 22); 	TrB=createNewNode(TrB, c(18,19), c(1,1), 23); 	TrB=createNewNode(TrB, c(8,20), c(1,1), 24); 	TrB=createNewNode(TrB, c(14,21), c(1,1), 25);
TrB=createNewNode(TrB, c(17,22), c(1,1), 26);  	TrB=createNewNode(TrB, c(5,23), c(1,1), 27);   	TrB=createNewNode(TrB, c(9,24),  c(1,1), 28); 	TrB=createNewNode(TrB, c(25,26), c(1,1), 29);
TrB=createNewNode(TrB, c(11,29), c(1,1), 30); 	TrB=createNewNode(TrB, c(10,30), c(1,1), 31); 	TrB=createNewNode(TrB, c(28,31), c(1,1), 32); 	TrB=createNewNode(TrB, c(27,32), c(1,1), 33);
TreeB=convertPhylo(TrB)
TreeBG=as.igraph(TreeB)
AdjA = get.adjacency(TreeAG, sparse = FALSE)
AdjB = get.adjacency(TreeBG, sparse = FALSE)
AdjA = AdjA + t(AdjA)
AdjB = AdjB + t(AdjB)
vec1 = eigen(AdjA)$values
vec2 = eigen(AdjB)$values
all.equal.numeric(vec1, vec2)
LapA = diag(rowSums(AdjA)) - AdjA
LapB = diag(rowSums(AdjB)) - AdjB
Vec1 = eigen(LapA)$values
Vec2 = eigen(LapB)$values
all.equal.numeric(Vec1, Vec2)
allPairsA = distances(TreeAG, mode = "all")
allPairsB = distances(TreeBG, mode = "all")
modLapA = diag(rowSums(allPairsA)) - allPairsA
modLapB = diag(rowSums(allPairsB)) - allPairsB
modVec1 = eigen(modLapA)$values
modVec2 = eigen(modLapB)$values
all.equal.numeric(modVec1, modVec2)
distA = allPairsA[as.character(1:17), as.character(1:17)]
distB = allPairsB[as.character(1:17), as.character(1:17)]
dVec1 = eigen(distA)$values
dVec2 = eigen(distB)$values
all.equal.numeric(dVec1, dVec2)
fVec1 = eigen(allPairsA)$values
fVec2 = eigen(allPairsB)$values
all.equal.numeric(fVec1, fVec2)
V1 = allPairsA["Node1",]
V2 = allPairsB["Node1",]
gVec1 = eigen(rbind(cbind(allPairsA, rep(1, 33)), c(rep(1, 33), 0)))$values
gVec2 = eigen(rbind(cbind(allPairsB, rep(1, 33)), c(rep(1, 33), 0)))$values
all.equal.numeric(gVec1, gVec2)
hVec1 = eigen(rbind(cbind(allPairsA, rep(1, 33)), c(V1, 0)))$values
hVec2 = eigen(rbind(cbind(allPairsB, rep(1, 33)), c(V2, 0)))$values
all.equal.numeric(hVec1, hVec2)
kVec1 = eigen(rbind(cbind(allPairsA, V1), c(V1, 0)))$values
kVec2 = eigen(rbind(cbind(allPairsB, V2), c(V2, 0)))$values
all.equal.numeric(kVec1, kVec2)

### Producing figures 2A and 2B
WantedColumns = rbind(c("between","betweenW","closeness","closenessW","eigen","eigenW"),
                      c("diameter","meanpath","minAdj","maxAdj","minLap","maxLap"),
                      c("cherries","pitchforks","doubcherries","fourprong","num5","num6"),
                      c("colless","sackin","maxwidth","maxheight","stairs","delW"))
MODIFIED = TRUE
Res0 = metaDriverOld("hivdenmeas.Rdata", decreasing = FALSE, K = 0)
save(Res0, file = "ThreeVirusesSummaryEigs.Rdata")
Res1 = metaDriverNew("hivdenmeas.Rdata", decreasing = FALSE, K = 0)
save(Res1, file = "ThreeVirusesSummaryMetrics.Rdata")
Res00 = metaDriverOld("NineViruses.Rdata", decreasing = FALSE, K = 0)
save(Res00, file = "NineVirusesSummaryEigs.Rdata")
Res11 = metaDriverNew("NineViruses.Rdata", decreasing = FALSE, K = 0)
save(Res11, file = "NineVirusesSummaryMetrics.Rdata")
load("treeshapestatsall.Rdata")
allsecond = allsecond[order(allsecond$bugs),] # order by the bug name!
stopifnot(all(unique(allsecond$bugs)==names(Res0[[1]]))) # checking bug order is correct!
Res0Red = lapply(Res0[1:2], function(x) {y = do.call(rbind, x); y = y[,c("min", "max")]})
Res0Red = do.call(cbind, Res0Red)
colnames(Res0Red) = as.vector(outer(c("min","max"), c("Adj","Lap"), paste0))
Res1Red = lapply(Res1[c(1:3, 4, 11, 6, 13, 8, 15)], function(x) {y = do.call(rbind, x); if (ncol(y) > 1) {y = y[,"max"]} else {y}})
Res1Red = do.call(cbind, Res1Red)
colnames(Res1Red) = c(c("diameter","meanpath","assort"), as.vector(outer(c("","W"), c("between","closeness","eigen"), function(x,y) {paste0(y,x)})))
allsecondExt = cbind(allsecond, Res0Red, Res1Red)
fdata = allsecondExt
if (MODIFIED) {
  fdata = fdata[,c("bugs", WantedColumns)]
}
plotlist=lapply(2:length(colnames(fdata)), function(x) ggplot(data=fdata)+aes_string(y=colnames(fdata)[x],x=colnames(fdata)[1],color=colnames(fdata)[1])+geom_boxplot()+theme(text = element_text(size=6) ,axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank())+guides(color=FALSE))
pdf("../Plots/AllMetricsA.pdf"); multiplot(plotlist=plotlist,cols = 7 - MODIFIED); dev.off()
allfirst = allfirst[order(allfirst$bugs),] # order by the bug name!
stopifnot(all(unique(allfirst$bugs)==names(Res00[[1]]))) # checking bug order is correct!
Res00Red = lapply(Res00[1:2], function(x) {y = do.call(rbind, x); y = y[,c("min", "max")]})
Res00Red = do.call(cbind, Res00Red)
colnames(Res00Red) = as.vector(outer(c("min","max"), c("Adj","Lap"), paste0))
Res11Red = lapply(Res11[c(1:3, 4, 11, 6, 13, 8, 15)], function(x) {y = do.call(rbind, x); if (ncol(y) > 1) {y = y[,"max"]} else {y}})
Res11Red = do.call(cbind, Res11Red)
colnames(Res11Red) = c(c("diameter","meanpath","assort"), as.vector(outer(c("","W"), c("between","closeness","eigen"), function(x,y) {paste0(y,x)})))
allfirstExt = cbind(allfirst, Res00Red, Res11Red)
fdata = allfirstExt
if (MODIFIED) {
  fdata = fdata[,c("bugs", WantedColumns)]
}
fluonly = fdata[grep("flu", fdata$bugs), ] # change made by LC for simplicity
plotlistflu=lapply(2:length(colnames(fdata)), function(x) ggplot(data=fluonly)+aes_string(y=colnames(fdata)[x],x=colnames(fdata)[1],color=colnames(fdata)[1])+geom_boxplot()+theme(text = element_text(size=6) ,axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank())+guides(color=FALSE))
pdf("../Plots/AllMetricsB.pdf"); multiplot(plotlist=plotlistflu,cols=7-MODIFIED); dev.off()
save(allfirstExt, allsecondExt, file = "ExtendedDataForCaroline.Rdata")

### This function creates an igraph representation of the tree we use as the running example
constructExample = function() {
  DF = rbind(c("A","a",8), c("A","B",4), c("B","C",2), c("B","D",2), c("C","b",2), c("C","c",2), 
             c("D","E",1), c("D","F",1), c("E","d",1), c("E","e",1), c("F","f",1), c("F","g",1))
  G = graph_from_data_frame(DF)
  E(G)$weight = as.numeric(DF[,3])
  Adir = get.adjacency(G, sparse = FALSE)
  Adir = Adir + t(Adir)
  print(eigen(Adir)$values)
  Adw = get.adjacency(G, attr="weight", sparse = FALSE)
  Adw = Adw + t(Adw)
  print(eigen(Adw)$values)
  ### CONTINUE WITH OTHER SPECTRA AND OTHER METRICS AS NEEDED ###
}