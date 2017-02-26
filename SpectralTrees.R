require(igraph)
require(ape)
require(phytools)

prepareFile = function(fileList, abbrList, outFileName = "PreparedTrees.RData") {
	for (ind in 1:length(fileList)) {
		filename = fileList[ind]
		print(paste("Processing", filename))
		abbr = abbrList[ind]
		allTrees = read.newick(filename)
		assign(abbr, allTrees)
	}
	save(list = abbrList, file = outFileName)
}

metaDriverOld = function(inputDataFile, K = 5, sep = ",", decreasing = TRUE, test = FALSE) {
  fullResults = c()
  count = 0
  for (directed in c(FALSE, TRUE)) {
    for (weighted in c(FALSE, TRUE)) {
      for (Laplacian in c(FALSE, TRUE)) {
        for (norm in c(FALSE, TRUE)) {
          for (mode in c("in", "out")) {
            if ((!Laplacian && norm) || (!directed && mode == "out")) {
              count = count + 1	
              print(count)
            }
            else {
              curResult = driverOld(inputDataFile, K = K, directed = directed, Laplacian = Laplacian, 
                                    norm = norm, weighted = weighted, mode = mode, decreasing = decreasing)
              curTitle = paste(c(directed, Laplacian, norm, weighted, mode), collapse = sep)
              print(curTitle)
              fullResults[[curTitle]] = curResult
            }
          }
        }
      }
    }
  }
  output = fullResults
  if (test) {
    pValues = performTests(fullResults, sep = sep)
    output = list(results = fullResults, pValues = pValues)
  }
  output
}

driverOld = function(inputDataFile, K = 5, directed = FALSE, Laplacian = TRUE, norm = FALSE, weighted = FALSE, mode = "in", decreasing = TRUE) {
  e = new.env()
  load(inputDataFile, e)
  Results = list()
  count = 1
  for (name in ls(e)) {
    print(paste("Processing dataset number", count))
    Item = get(name, e)
    procItem = lapply(Item, function(x) {x$edge})
    W = NULL
    if (weighted) {
      W = lapply(Item, function(x) {x$edge.length/sum(x$edge.length)*length(x$edge.length)})
    }
    curRes = processItem(procItem, K = K, directed = directed, Laplacian = Laplacian, norm = norm, weights = W, mode = mode, decreasing = decreasing)
    Results[[name]] = curRes
    count = count + 1
  }
  Results
}

metaDriverNew = function(inputDataFile, K = 5, sep = ",", decreasing = TRUE) {
	fullResults = c()
	count = 0
	for (directed in c(FALSE, TRUE)) {
		for (weighted in c(FALSE, TRUE)) {
			for (Laplacian in c("diameter", "meanpath", "assortativity", "betweenness", "closeness", "evcent")) {
				for (norm in c(FALSE, TRUE)) {
					if ((Laplacian %in% c("diameter","meanpath","assortativity") && norm) || (Laplacian %in% c("betweenness", "closeness","evcent") && directed) || (Laplacian %in% c("assortativity","meanpath") && weighted)) {
						count = count + 1	
						print(count)
					}
					else {
						curTitle = paste(c(Laplacian, directed, norm, weighted), collapse = sep)
						print(curTitle)
						curResult = driverNew(inputDataFile, K = K, directed = directed, Laplacian = Laplacian, 
												norm = norm, weighted = weighted, decreasing = decreasing)
						fullResults[[curTitle]] = curResult
					}
				}
			}
		}
	}
	fullResults
}

performTests = function(ListOfLists, hyps = 1:5, firsts = c(1, 2), seconds = c(3, 4), tol = 1e-8, file = "Results.pdf", pMax = 0.05, sep = ",", metric = FALSE) {
	pdf(file = file)
	Len = length(ListOfLists)
	allPValues = vector("list", Len)
	names(allPValues) = names(ListOfLists)
	m = length(hyps)
	n = length(firsts)
	for (ind in 1:Len) {
		curList = ListOfLists[[ind]]
		curMatrix = matrix(NA, m, n)
		for (ind2 in 1:n) {
			curFirst = curList[[firsts[ind2]]]
			curSecond = curList[[seconds[ind2]]]
			for (ind1 in 1:m) {
				curHyp = hyps[ind1]
				if (curHyp > 0  && curHyp <= ncol(curFirst)) {
					curFirstVec  = curFirst [, curHyp]
					curSecondVec = curSecond[, curHyp] 
					if ((var(curFirstVec) > tol && var(curSecondVec) > tol)) {
						curMatrix[ind1, ind2] = t.test(curFirstVec, curSecondVec)$p.value
					}
				}
			}
		}
		rownames(curMatrix) = hyps
		colnames(curMatrix) = paste(firsts, seconds, sep = "vs")
		goodHyps = apply(curMatrix < pMax, 1, all)
		if(any(!is.na(goodHyps) & goodHyps)) {
			curName = names(ListOfLists)[ind]
			fullName = convertName(curName, sep = sep, metric = metric)
			for (Hyp in which(goodHyps)) {
				boxplot(as.data.frame(curList)[Hyp + 0:(max(firsts, seconds) - 1) * max(hyps)], main = fullName)
			}
		}
		allPValues[[ind]] = curMatrix
	}
	dev.off()
	allPValues
}

convertName = function(Name, sep, spectrum = FALSE, reduced = FALSE, metric = FALSE, splitBy = " ") {
	splitName = unlist(strsplit(Name, split = sep))
	lastPart = splitName[5]
	if (metric) {
		metName = splitName[1]
	}
	splitName = as.logical(splitName[-c(metric, 5)])
	if (reduced) {
		part0 = ""
		part2 = paste0(rep("Und", !splitName[1]), rep("Dir",  splitName[1]))
		part1 = paste0(rep("Unn", !splitName[3-metric]), rep("Norm", splitName[3-metric]))
		part3 = paste0(rep("Unw", !splitName[4-metric]), rep("Wgt",  splitName[4-metric]))
		part4 = ""
		if (!metric) {
			part4 = paste0(rep("Lap",  splitName[2]), rep("Adj", !splitName[2]))
		}
		else {
			part0 = paste0(toupper(substr(metName,1,1)), substr(metName,2,min(nchar(metName)/2, 4)))
		}
		if (!metric && splitName[3]) {
			part1 = paste0(part1, toupper(substr(lastPart,1,1)), substr(lastPart,2,nchar(lastPart)))
			fullName = paste0(rep(part0, metric), part2, part3, rep(part4, !metric), part1)
		}
		else {
			fullName = paste0(rep(part0, metric), part1, part2, part3, rep(part4, !metric))
		}
	}
	else {
		if (metric) {
			part4 = metName
		}
		else {
			part4 = paste0(rep("Laplacian", splitName[2]), rep("adjacency", !splitName[2]))
		}
		part2 = paste0(rep("un", !splitName[1]), "directed")
		part1 = paste0(rep("un", !splitName[3-metric]), "normalized")
		if (!metric && splitName[3]) {
			part1 = paste0(part1, " by ", lastPart, "degree")
		}
		part3 = paste0(rep("un", !splitName[4-metric]), "weighted")
		if (!metric && splitName[3]) {
			fullName = paste(paste(part2, part3), part4, sep = splitBy)
			if (spectrum) {
				fullName = paste("Spectrum of", fullName, "matrix")
			}
			fullName = paste(fullName, part1)
		}
		else {
			fullName = paste(paste(part1, part2), paste(part3, part4), sep = splitBy)
			if (spectrum) {
				fullName = paste("Spectrum of", fullName, "matrix")
			}
		}
	}
	fullName
}

driverNew = function(inputDataFile, K = 5, directed = TRUE, Laplacian = NULL, norm = TRUE, weighted = TRUE, mode = "in", decreasing = TRUE) {
	e = new.env()
	load(inputDataFile, e)
	Results = list()
	count = 1
	for (name in ls(e)) {
		print(paste("Processing dataset number", count))
		Item = get(name, e)
		procItem = lapply(Item, function(x) {x$edge})
		W = NULL
		if (weighted) {
			W = lapply(Item, function(x) {y=x$edge.length/sum(x$edge.length)*length(x$edge.length);y[y==0]=min(y[y>0]);y})
		}
		curK = K
		if (Laplacian %in% c("diameter", "meanpath", "assortativity")) {
			curK = 1
		}
		curRes = processItem(procItem, K = curK, directed = directed, Laplacian = Laplacian, norm = norm, weights = W, mode = mode, metric = TRUE, decreasing = decreasing)
		Results[[name]] = curRes
		count = count + 1
	}
	Results
}

processItem = function(List, K, directed, Laplacian, norm, weights, mode, metric = FALSE, decreasing = TRUE) {
	L = length(List)
	if (K == 0) {
		output = matrix(NA, L, 3)
		colnames(output) = c("min", "mean", "max")
	}
	else {
		output = matrix(NA, L, K)
		colnames(output) = paste("E", 1:K, sep = "")
	}
	print(paste("There are", L, "graphs to process"))
	for (ind in 1:L) {
		curG = List[[ind]]
		curW = weights[[ind]]
		if (!metric) {
			output[ind,] = findEigs(curG, K = K, directed = directed, Laplacian = Laplacian, norm = norm, weight = curW, mode = mode, decreasing = decreasing)
		}
		else {
			output[ind,] = findMetrics(curG, K = K, directed = directed, norm = norm, weight = curW, metric = Laplacian, decreasing = decreasing)
		}
	}
	output
}

findEigs = function(edgeList, K, directed, Laplacian, norm, weight, mode, decreasing = TRUE, EPS = 1e-15) { # change to FALSE if needed!
	if (!is.null(weight)) {
		df = as.data.frame(cbind(edgeList, weight = weight))
		curGraph = graph.data.frame(df, directed = directed)
		curMat = get.adjacency(curGraph, sparse = FALSE, attr = "weight")
	}
	else {
		curGraph = graph(t(edgeList), directed = directed)
		curMat = get.adjacency(curGraph, sparse = FALSE)
	}
	if (Laplacian) {
		if (mode == "in") {
			degrees = colSums(curMat)
		}
		else {
			degrees = rowSums(curMat)
		}
		curMat = diag(degrees) - curMat
		if (norm) {
			curMat = curMat / sqrt(outer(degrees, degrees))
		}
	}
	curMat[is.infinite(curMat)] = 0
	curMat[is.na(curMat)] = 0
	curEigs = abs(eigen(curMat, symmetric = FALSE, only.values = TRUE)$values)
	if (K == 0) {
		values = c(min(curEigs[curEigs > EPS]), mean(curEigs), max(curEigs))
	}
	else {
		values = sort(curEigs, decreasing = decreasing)[1:K]
	}
	values
}

findMetrics = function(edgeList, K, directed, norm, weight, metric, decreasing = TRUE) { # change to FALSE if needed!
	if (!is.null(weight)) {
		df = as.data.frame(cbind(edgeList, weight = weight))
		curGraph = graph.data.frame(df, directed = directed)
	}
	else {
		curGraph = graph(t(edgeList), directed = directed)
	}
	if (metric == "diameter") {
		value = diameter(curGraph, directed = directed)
	}
	else if (metric == "meanpath") {
		value = average.path.length(curGraph, directed = directed) 
	}
	else if (metric == "assortativity") {
		value = assortativity.degree(curGraph, directed = directed)
	}
	else { # centrality metric - return several values
		if (metric == "betweenness") {
			values = betweenness(curGraph, normalized = norm)
		}
		else if (metric == "closeness") {
			values = closeness(curGraph, normalized = norm)
		}
		else if (metric == "evcent") {
			values = evcent(curGraph, scale = norm)$vector
		}
		if (K == 0) {
			value = c(min(values), mean(values), max(values))
		}
		else {
			value = sort(values, decreasing = decreasing)[1:K]
		}
	}
	value
}