#script adapted from Shahan, Rachel, et al. "Consensus coexpression network analysis identifies key regulators of flower and fruit development in wild strawberry." Plant physiology 178.1 (2018): 202-216.
#usage: Rscript consensus_cluster.R <min module size> <full path to consensus matrix> <save directory>

library(WGCNA)
library(data.table)

enableWGCNAThreads()

print("loading libraries")

#setup
args<-commandArgs(TRUE)

#gather command line arguments
if (length(args) < 3){
	message("arguments are insufficient")
	quit(save = "no")
} else {
	minModuleSize <- as.numeric(args[1])
	consensusMat_fname <- as.character(args[2])#includes full path
	save_dir <- as.character(args[3])
}

zz <- file(paste0(save_dir,"/log.log"), open="wt")
sink(zz, type = "message")

file.create(paste(save_dir,"/failed", sep =""))

setwd(paste(save_dir,"/..", sep = ""))

adjacency_res <- as.matrix(as.data.frame(fread(consensusMat_fname, header = TRUE))) 
rownames(adjacency_res) <- colnames(adjacency_res)
dist_mat <- 1 - adjacency_res

keepList <- colnames(adjacency_res)

print("building gene Tree")

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dist_mat), method = "average")

# Module identification using dynamic tree cut:
print("cutting modules")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dist_mat, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose= 4)
print("module sizes")
table(dynamicMods)
print("writing module assignments")
write.csv(as.character(dynamicMods), file= paste0(save_dir,"/dynamicMods.csv"))



# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
l <- length(unique(dynamicColors))
table(dynamicColors)

if(l > 2000){		
	file.remove(paste(save_dir, "/failed", sep =""))
	print("too many clusters")
	file.create(paste(save_dir, "/2many", sep =""))
	quit(save = "no")
}

#write out clusters
gene_color <- data.frame(gene = rownames(dist_mat), cluster = dynamicColors)
gene_color <- gene_color[which(gene_color$cluster != "grey"),]
gene_num <- data.frame(gene = rownames(dist_mat), cluster = dynamicMods)
gene_color$cluster <- as.character(gene_color$cluster)

write.csv(gene_color, file = paste0(save_dir, "/gene_cluster_color.csv"), row.names = FALSE)
write.csv(gene_num, file = paste0(save_dir, "/gene_cluster_num.csv"), row.names = FALSE)

###write out connectivity stats
#intramodularConnectiveity with dist_mat based on expression values rather than adjacency mat
exprData = read.csv("Kronos_DEGs_tpms_standard_input_allreps.csv")
datExpr0 = as.data.frame(t(exprData[, -(1)]))
datExpr0 = data.matrix(datExpr0)
names(datExpr0) = exprData$Geneid
rownames(datExpr0) = names(exprData)[-(1)]

alldegrees2 <- intramodularConnectivity.fromExpr(datExpr0, 
                                                dynamicMods,
                                                networkType = "signed"
                                                scaleByMax = TRUE)

sink("scaled_intramodular_connectivity_TPMDEG_FROMEXPR.txt")
alldegrees2
sink()

modMEs <- moduleEigengenes(datExpr0,
							dynamicMods,
							align = 'along average')
sink("moduleEigengenes.csv")
modMEs$eigengenes
sink()
sink("moduleAvgExprNorm.csv")
modMEs$averageExpr
sink()

########
########
########

file.remove(paste(save_dir,"/failed", sep =""))
file.create(paste(save_dir, "/success",sep =""))

for (module in dynamicMods)
{
	x <- names(datExpr0)[dynamicMods$colors==module]
	y <- length(x)
	
	write.table(x, paste("mod_probes/mod_",module,"__",y,"_probes.csv",sep=""), sep="\t")
	}
