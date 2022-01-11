#This script builds adjacency matrix using output files in '1000_output/'
#A) counts how often each gene-pair was subsampled
#B) counts how many times each gene-pair clusters in the same module
#C) calculates adjacency scores to fill matrix using each gene-pairs co-clustered and co-sampled score 

import csv
from itertools import izip


allGeneList = [] #list containing every gene from original source file
sampleCompendium = {} 
moduleCompendium = {} 

#Populating allGeneList:
with open('overlap_DEGs_tpms_standard_input_allreps.csv','r') as f:
	reader = csv.reader(f,delimiter=',',quotechar='"')
	reader.next()
	for row in reader:
		allGeneList.append(row[0])

		#initialize compendiums
		sampleCompendium[row[0]] = {}		
		moduleCompendium[row[0]] = {}		

print 'gene compendium declared' 

#Initialize dictionaries in compendiums
for gene in allGeneList:
	for gene2 in allGeneList:
		sampleCompendium[gene][gene2] = 0 
		moduleCompendium[gene][gene2] = 0 
	print gene							  


print 'gene compendiums initialized' 

#zip thru the file pair for each wgnca run, add both files into memory to read thru multiple times
for i in range(1000): 						
	sampleGenes = set() 					
	geneModules = [] 						

	with open('1000_output/' + str(i) + 'genes.csv', 'r') as geneFile, open('1000_output/' + str(i) + 'modules.csv', 'r') as modFile:
		geneReader = csv.reader(geneFile, delimiter=',', quotechar='"')
		modReader = csv.reader(modFile, delimiter=',', quotechar='"')

		geneReader.next()
		modReader.next()

		for gene, module in izip(geneReader, modReader): 	
			geneModules.append((gene[1], module[1]))		
			sampleGenes.add(gene[1])

	for gene in sampleGenes:
		print 'sample', gene, '...' 
		for gene2 in sampleGenes:
			sampleCompendium[gene][gene2] += 1	
			
	moduleMap = {} 					
	modules = set()					
	for geneModule in geneModules:				
		if geneModule[1] not in modules:			
			moduleMap[geneModule[1]] = set()		
			modules.add(geneModule[1])			
		moduleMap[geneModule[1]].add(geneModule[0])	
		
	for module in modules:						
		print 'module', module, '...'
		for gene in moduleMap[module]:			
			print 'group gene', gene 	
			for gene2 in moduleMap[module]:				
				moduleCompendium[gene][gene2] += 1	


#do all the division in compendiums and make a new dictionary with the quotients
geneMatrix = {}						
for gene in allGeneList:
	geneMatrix[gene] = {}
	for gene2 in allGeneList:																								#iterate through all gene-gene2 pairs
		#account for divide by zero errors
		if sampleCompendium[gene][gene2] != 0:																				#check if gene-gene2 pair is sampled at least once
			geneMatrix[gene][gene2] = float(moduleCompendium[gene][gene2]) / float(sampleCompendium[gene][gene2])		
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
