#This script is for obtaining TPM counts from a raw read count matrix
#Feature (gene) length is first determined by cumulative exonic length for estimating transcript abundance

import csv
import os
import pandas as pd

countfile = 'kronos_counts.txt'	
outfile = 'Kronos_TPM.txt'		

geneLength = {}

####getting gene length in kilobases,
gffDir = 'Annotations_v1.1/iwgsc_refseqv1.1_genes_2017July06/'
gffFiles = ['IWGSC_v1.1_HC_20170706.gff3', 'IWGSC_v1.1_LC_20170706.gff3']

for filename in gffFiles:
	with open(gffDir+filename, 'r') as f:
		reader = csv.reader(f, delimiter='\t', quotechar='"')

		currGene = ''
		currLen = 0.0
		for row in reader:
			if len(row) > 4:	#skip rows without annotation features
				if row[2] == 'gene':
					if currGene != '':
						geneLength[currGene] = currLen

					currGene = row[8][3:row[8].index(';')]		#getting the name of the gene
					currLen = 0.0

				if row[2] == 'exon':
					currLen += ((float(row[4]) - float(row[3])) / 1000)	#adding to the read-mappable length of the gene only by exon lengths

		geneLength[currGene] = currLen
geneLengthSet = set(geneLength.keys())

####convert the geneLength dictionary to a pd.dataframe
geneLength_df = pd.DataFrame.from_dict(geneLength, orient = 'index')

####read in raw count file from FeatureCounts as pandas.dataframe
readCount_df = pd.read_table(countfile, sep = '\t', index_col = 'Geneid')

####sort gene length dataframe so it aligns with the read count data
geneLength_df = geneLength_df.reindex(index = readCount_df.index)

####produce RPK table: dividing raw counts table by gene length
RPK_df = readCount_df.divide(geneLength_df[0], axis = 0)

####get scaling factor by dividing total RPKs per sample by 1,000,000
scalingFactors = RPK_df.sum(axis = 0)/1000000

####calculate TPM table by dividing RPK table by scaling factor
TPM_df = RPK_df.divide(scalingFactors, axis = 1)

####print TPM_df to outdir
TPM_df.to_csv(outdir, sep = "\t")
