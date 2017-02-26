library(ggplot2)
library(plyr)
library(gridExtra)

# setwd("/Users/leonidc/Documents/HSPHwork/Project with Caroline/newick/") # replace as needed

plotResults = function (sourceFile, varName, metric = FALSE, unitary = FALSE, test = FALSE) {
	dirName = gsub(".Rdata", "", sourceFile)
	if (!(metric && unitary)) {
		dir.create(dirName)
	}
	e = new.env()
	load(sourceFile, e)
	fullResults = get(varName, envir = e)
	if (metric) {
		badInds = c(1:3, 10, 17:20) # unitary metrics
		if (unitary) {
			fullResults = fullResults[badInds]
		}
		else {
			fullResults = fullResults[-badInds]
		}
	}
	else {
		badInds = c(2,5,9,10,15,16) # unnormalized Laplacians
		fullResults = fullResults[-badInds]
	}
	if (test) {
		Q = performTests(fullResults, firsts = 1:3, seconds = 4:6, pMax = 0.01/length(fullResults), metric = TRUE) 
	}
	rearrangedResults = lapply(fullResults, function(item) {do.call(rbind,item)}) # replaces the application of getmatpart
	mynames = names(fullResults[[1]])
	Lengths = sapply(fullResults[[1]], nrow)
	Virus = rep(mynames, Lengths)
	augmentedResults = lapply(rearrangedResults, function(x) {data.frame(x, Virus = Virus)}) # creating augmented data frames
	fullNames = sapply(names(fullResults), function(x) {convertName(x, sep = ",", spectrum = FALSE, reduced = FALSE, metric = metric, splitBy = "\n")})
	names(augmentedResults) = fullNames
	L = length(augmentedResults)
	if (metric && unitary) {
		RANGE = 1
		NCOL = 2
		SEP = "All"
	}
	else {
		RANGE = 1:5
		NCOL = 3
		SEP = "E"
	}
	for (ind in RANGE) {
		p = ggplot()
		plotlist = vector("list", L)
		for (j in 1:L) {
			curData = augmentedResults[[j]]
			curName = names(augmentedResults)[j]
			plotlist[[j]] = makeggsubplot2(p, "Virus", colnames(curData)[ind], Special = TRUE, data = curData, Title = curName)
		}
		if (metric && unitary) {
			plotlist[[7]] = plotlist[[8]]
			plotlist = plotlist[1:7]	
		}
		g <- do.call("arrangeGrob", c(plotlist, ncol = NCOL))
		ggsave(paste0(dirName, "/", "fullResults", SEP, ind, ".pdf"), plot = g)
	}
}

makeggsubplot2 = function(p, xstring, ystring, ylabel=NULL, Special=FALSE, data=NULL, Title=NULL) { # Caroline's script 
   if (Special==FALSE) {
   		plot1 <-p + geom_boxplot(aes_string(x=xstring,y=ystring,color=xstring), data = data)+guides(colour=FALSE)+labs(x=NULL,y=ylabel,title=Title) + theme(plot.title = element_text(size = 10))
   } 
   else {  
   		plot1 <-p + geom_boxplot(aes_string(x=xstring,y=ystring,color=xstring), data = data)+theme(axis.text.x=element_text(angle=45, hjust=1))+guides(colour=FALSE)+labs(x=NULL,y=ylabel,title=Title) +theme(plot.title = element_text(size = 10))
   }
   return(plot1)
}