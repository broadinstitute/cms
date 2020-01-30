## quick script to visualize distributions of results from one or many files (simulated, empirical)
## last updated 11.16.16

oneFile, simFiles, empFiles = True, False, False

filebase = "/idi/sabeti-scratch/jvitti/scores_composite5/"
models = ['default_112115_825am']#, 'nulldefault', 'default_default_101715_12pm','nulldefault_constantsize', 'gradient_101915_treebase_6_best',]
selbins = ["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]
modelpops = [1, 2, 3, 4]
emppops = ['YRI' 'CEU', 'CHB', 'BEB']
chroms = range(1,23)

import numpy as np
import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as plt
import os

def plot_dist(allvals, savefilename= "/web/personal/vitti/test.png", numBins=10000):
	if len(allvals) > 0:
		f, ax = plt.subplots(1)
		ax.hist(allvals, bins=10000)
		plt.savefig(savefilename)
		print('plotted to ' + savefilename)
	return
def readvals_lastcol(filename):
	allvals = []
	openfile = open(filename, 'r')
	for line in openfile:
		entries = line.split()
		value=float(entries[-1])#float(entries[4]) #thisihs, thisihh, thisxpehh, thisfst, thisdelDaf,
		if not np.isnan(value):
			if not np.isnan(np.log(value)):
				allvals.append(np.log(value))
	openfile.close()
	return allvals

def main():
	############################
	## VIEW DIST FROM 1 FILE  ##
	############################
	if oneFile:
		savefilename = "/web/personal/vitti/test2.png"
		infilename = "/idi/sabeti-scratch/jvitti/scores_composite_5/nulldefault/neut/rep1000_4.cms.out"
		#"/idi/sabeti-scratch/jvitti/scores_composite_5/nulldefault/neut/rep1000_4.cms.out"
	#"/idi/sabeti-scratch/jvitti/scores_composite4_b/CHB_gradient_101915_treebase_6_best_neut.chr10.txt.norm"
		#'/idi/sabeti-scratch/jvitti/scores_composite4_b/default_112115_825am/neut/rep500_2.cms.out'
		#'/idi/sabeti-scratch/jvitti/scores_composite4/default_112115_825am/sel1/sel_0.10/rep47_1_vsneut.cms.out'
		if os.path.isfile(infilename):
			values = readvals_lastcol(infilename)
			plot_dist(values, savefilename)


	#########################
	## VIEW DIST FROM SIMS ##
	#########################
	elif simFiles:
		for model in models:		
			for pop in modelpops:
				savefilename = "/web/personal/vitti/" + model +"_" + str(pop) + ".png"
				allvals = []

				#NEUT REPS
				for irep in range(1,501):
					infilename = filebase + model + "/neut/rep" + str(irep) + "_" + str(pop) + ".cms.out"
					if os.path.isfile(infilename):
						values = readvals_lastcol(infilename)
						allvals.extend(values)

				plot_dist(allvals, "/web/personal/vitti/" + model +"_" + str(pop) + "_justneut.png")

				#SEL REPS
				for selbin in selbins:	
					for irep in range(1,501):
						infilename = filebase + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/rep" + str(irep) + "_" + str(pop) + "_vsneut.cms.out"
						if os.path.isfile(infilename):
							values = readvals_lastcol(infilename)
							allvals.extend(values)
			
				plot_dist(allvals, savefilename)

	#############################
	## VIEW DIST FROM EMP DATA ##
	#############################
	elif empFiles:
		for model in models:
			for pop in emppops:
				savefilename = "/web/personal/vitti/" + model +"_" + str(pop) + ".png"
				allvals = []
				for chrom in chroms:
					infilename = filebase + pop + "_" + model + "_" + "neut" + ".chr" + str(chrom) + ".txt"
					if os.path.isfile(infilename):
						values = readvals_lastcol(infilename)
						allvals.extend(values)
				plot_dist(allvals, savefilename)

main()
