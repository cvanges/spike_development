#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This runs wgcna with command line arguments passed by 'WGCNA_adjacency_wrapper.py'
#Parameters: subsampled gene file, powers, minModuleSize, merge eigengene, run #
#ex: Rscript wgcna_args.R genes.csv 12 60 1 999


library(WGCNA)

options(stringsAsFactors = FALSE)
options(max.print=500000)

exprData = read.csv(args[1]) #filename is first argument

datExpr0 = as.data.frame(t(exprData[, -(1)]))
names(datExpr0) = exprData$Geneid
rownames(datExpr0) = names(exprData)[-(1)]

enableWGCNAThreads()

bwnet = blockwiseModules(datExpr0, 
						 maxBlockSize = 30000, 
						 power = as.integer(args[2]), 
						 TOMType = "signed", 
						 networkType = "signed", 
						 corType = "bicor", 
						 maxPOutliers = 0.05, 
						 minModuleSize = as.integer(args[3]), 
						 reassignThreshold = 0, 
						 deepSplit = 2, 
						 mergeCutHeight = as.numeric(args[4]), 
						 numericLabels = TRUE, 
						 saveTOMs = FALSE, 
						 verbose = 3)
						 
						 

bwLabels = bwnet$colors
bwModuleColors = labels2colors(bwLabels)

probes = names(datExpr0)
modules = unique(bwLabels)
inModule = is.finite(match(bwLabels, modules))
modProbes = probes[inModule]


write.csv(modProbes, file = paste("1000_output/", args[5], "genes.csv", sep = ""))
write.csv(bwLabels, file = paste("1000_output/", args[5], "modules.csv", sep = ""))
print(paste("COMPLETED", args[1], args[2], args[3], args[4], args[5], sep = " "))

#remove the sample file
file.remove(args[1])

