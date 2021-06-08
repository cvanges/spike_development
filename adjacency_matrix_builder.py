#This script builds adjacency matrix using output files in '1000_output/'
#A) counts how often each gene-pair was subsampled
#B) counts how many times each gene-pair clusters in the same module
#C) calculates adjacency scores to fill matrix using each gene-pairs co-clustered and co-sampled score 

import csv
from itertools import izip


allGeneList = [] #list containing every gene from original source file
sampleCompendium = {} #gene -> dictionary where gene -> 5 (number of times sampled together
moduleCompendium = {} #gene -> dictionary where gene -> 3 (number of times module'd together)

#Populating allGeneList:
with open('overlap_DEGs_tpms_standard_input_allreps.csv','r') as f:
	reader = csv.reader(f,delimiter=',',quotechar='"')
	reader.next()
	for row in reader:
		allGeneList.append(row[0])

		#initialize compendiums
		sampleCompendium[row[0]] = {}		#22566 items in sampleCompendium dict, each a geneid key with empty dictionary
		moduleCompendium[row[0]] = {}		#22566 items in moduleCompendium dict, each a geneid key with empty dictionary

print 'gene compendium declared' 

#Initialize dictionaries in compendiums
for gene in allGeneList:
	for gene2 in allGeneList:
		sampleCompendium[gene][gene2] = 0 #number of times sampled
		moduleCompendium[gene][gene2] = 0 #number of times module'd
	print gene							  #all possibilities are duplicated


print 'gene compendiums initialized' 

#zip thru the file pair for each wgnca run, add both files into memory to read thru multiple times
for i in range(1000): 						#all 1000 samples
	sampleGenes = set() 					#temporary variable that holds all the genes used in each 'i' wgcna run
	geneModules = [] 						#list containing pairs (gene, module), populated by zipping thru both wgcna outputs

	with open('1000_output/' + str(i) + 'genes.csv', 'r') as geneFile, open('1000_output/' + str(i) + 'modules.csv', 'r') as modFile:
		geneReader = csv.reader(geneFile, delimiter=',', quotechar='"')
		modReader = csv.reader(modFile, delimiter=',', quotechar='"')

		geneReader.next()
		modReader.next()

		for gene, module in izip(geneReader, modReader): 	#load both files all into memory at once, because have to iterate thru multiple times
			geneModules.append((gene[1], module[1]))		#gene and module are read as lists, each with two items. index , gene/module
			sampleGenes.add(gene[1])

	###for each wgcna run 'i', now have set of 18052 geneids + module assignment pairs in geneModules and seperate list of 18052 geneids in sampleGenes
	###first update sampleCompendium by iterating through sampleGenes twice and adding +1 to each gene pair combination in run 'i'
	for gene in sampleGenes:
		print 'sample', gene, '...' 
		for gene2 in sampleGenes:
			sampleCompendium[gene][gene2] += 1			###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WHERETHEMAGICHAPPNS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	#make individual sets for modules, do the same thing. another dictionary incoming
	moduleMap = {} 										#moduleMap['1'] -> set containing all genes in the module 1
	modules = set()										#modules set containing all clusters in run 'i'
	for geneModule in geneModules:						#geneModules are sets with geneModule[0] -> Geneid and geneModule[1] -> x (module)
		if geneModule[1] not in modules:				#checking if module x is a key in moduleMap dict yet
			moduleMap[geneModule[1]] = set()			#initializing module x as an empty set
			modules.add(geneModule[1])					#filling set of all modules in run 'i'
		moduleMap[geneModule[1]].add(geneModule[0])		#add geneid to set of genes in module x of wgcna run 'i'
		
	for module in modules:								#loop through all modules
		print 'module', module, '...'
		for gene in moduleMap[module]:					#loop through all genes in module x
			print 'group gene', gene 	
			for gene2 in moduleMap[module]:				#second nested loop through all genes in module x, so that gene and gene2 now refer to co-clustered gene pairs in run x
				moduleCompendium[gene][gene2] += 1		###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WHERETHEMAGICHAPPNS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#do all the division in compendiums and make a new dictionary with the quotients
geneMatrix = {}						#initialize compendium of 22566 geneids -> dictionary. The nested dictionary will point to the adjacency value for gene-gene2 pair
for gene in allGeneList:
	geneMatrix[gene] = {}
	for gene2 in allGeneList:																								#iterate through all gene-gene2 pairs
		#account for divide by zero errors
		if sampleCompendium[gene][gene2] != 0:																				#check if gene-gene2 pair is sampled at least once
			geneMatrix[gene][gene2] = float(moduleCompendium[gene][gene2]) / float(sampleCompendium[gene][gene2])			#set adjacency value for gene pair in geneMatrix to fraction of coclustered/subsampled
		else:
			geneMatrix[gene][gene2] = 0																						#if never subsampled (shouldnt happen by odds), then adjacency set to 0


#saving the full adjacency matrix
with open('consensus_adjacency_matrix.csv', 'w') as out:
	writer = csv.writer(out, delimiter=',', quotechar='"')

	writer.writerow(['gene'] + allGeneList)	
	
	for gene in allGeneList:
		temprow = [gene]
		for gene2 in allGeneList:
			temprow.append(geneMatrix[gene][gene2])
		writer.writerow(temprow)
