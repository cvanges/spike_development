#This script runs repeated, randomized WGCNA network constructions for constructing adjacency network
#First, this script will establish entire input data set and directories to colluct subsampled genes and cluster assignment for each network construction
#Second, main loop performs two functions 1,000 times 
###A) 80% of variables (gene expression) are subsampled and written to a temporary 'samples/' directory
###B) co-expression network constructed with randomly chosen power, minModSize, and merge passed to secondary script 'wgcna_args.R' which writes cluster assignment to '1000_output/'
#Temporary subsampled expression data is removed after each run, subsampled and co-clustered data used to construct adjacency matrix is stored in "-genes.csv" and "-modules.csv" in '1000_output/'

import subprocess
import os
import random

processes = set()
maxprocesses = 11

lines = []
with open('Kronos_expression.csv', 'r') as f:
	lines = f.readlines()

if not os.path.isdir('samples/'):
	subprocess.call(['mkdir', 'samples'])

if not os.path.isdir('1000_output/'):
	subprocess.call(['mkdir', '1000_output'])

for i in range(0,1000):#number of wgcna runs
	
	sample = random.sample(lines[1:], int((len(lines)-1) * .8)) #random 80% of genes

	with open('samples/sample'+str(i)+'.csv', 'w') as out:
		out.write(lines[0])
		for line in sample:
			out.write(line)

	
	powers = random.choice([1,2,4,8,12,16,20])
	minModuleSize = random.choice([40,60,90,120,150,180,210])
	merge = random.choice([0.15,0.2,0.25,0.3])

	processes.add(subprocess.Popen(['Rscript','wgcna_args.R', 'samples/sample'+str(i)+'.csv', str(powers), str(minModuleSize), str(merge), str(i)]))
	if len(processes) >= maxprocesses:
		os.wait()
		processes.difference_update([p for p in processes if p.poll() is not None])
	
