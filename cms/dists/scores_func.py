## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 09.22.16 vitti@broadinstitute.org

from math import fabs, sqrt
from random import randint
import numpy as np
import subprocess
import sys
import os

#####################
## CALCS FROM SIMS ## 
#####################
def calc_ihs(inputTped, outputFile, runProgram = "scans.py", numThreads = 7):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_ihs " + inputTped + " " + outputFile + " --threads " + str(numThreads)
	#print(cmdStr)
	subprocess.check_call( cmdStr.split() )
	return
def calc_delihh(readfilename, writefilename):
	"""given a selscan iHS file, parses it and writes an analogous file with delIHH information"""
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		entries = line.split()
		locus, phys, freq_1, ihh_1, ihh_0, ihs, ihh_0_l, ihh_0_r, ihh_1_l, ihh_1_r = entries #[0:10]
					#ancestral - derived
		unstand_delIHH = fabs(float(ihh_1) - float(ihh_0))
		writeline = locus + "\t" + phys + "\t" + freq_1 + "\t" + str(ihs) + "\t" + str(unstand_delIHH) +"\t" + str(unstand_delIHH) +  "\n" #6 columns for selscan norm
		writefile.write(writeline)
	writefile.close()
	readfile.close()
def calc_xpehh(inputTped, inputTped2, outputFile, runProgram = "scans.py", numThreads = 7):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_xpehh " + inputTped + " " + outputFile + " " + inputTped2 + " --threads " + str(numThreads)
	#print(cmdStr)
	subprocess.check_call( cmdStr.split() )
	return	
def calc_fst_deldaf(inputTped, inputTped2, recomFile, outputFile, modelpath):
	if modelpath[-1] != "/":
		modelpath += "/"
	commandstring = modelpath + "calc_fst_deldaf"
	argstring = inputTped + " " + inputTped2 + " " + recomFile + " " + outputFile
	fullcommand = commandstring + " " + argstring
	#print(fullcommand)
	subprocess.check_call( fullcommand.split() )
	return
def read_neut_normfile(neutNormfilename, scoretype ='ihs'):
	"""pulls means and vars so that we can normalize sel scenarios with the same parameters. pulled from func_clean.py"""
	openfile = open(neutNormfilename)
	line = openfile.readline()
	entries = line.split()
	while len(entries) < 1 or entries[0] != "Total":
		line = openfile.readline()
		entries = line.split()
	nsnp = int(entries[-1])
	while "num" not in entries:
		line = openfile.readline()
		entries = line.split()
	if scoretype == "ihs":
		assert entries[0] == "bin"
		bins, nums, means, variances = [], [], [], []
		for line in openfile:
			#if line is not None:
			entries = line.split()
			if (len(entries) >=1):
				if entries[0] != "Normalizing":
					binlimit, numinbin, mean, variance = float(entries[0]), int(entries[1]), float(entries[2]), float(entries[3])
					bins.append(binlimit)
					nums.append(numinbin)
					means.append(mean)
					variances.append(variance)
		openfile.close()
		return bins, nums, means, variances
	elif scoretype == "xp":
		dataline = openfile.readline()
		entries = dataline.split()
		num, mean, var = int(entries[0]), float(entries[1]), float(entries[2])
		#print(str(num) + "\t" +  str(nsnp))
		#assert num == nsnp
		return num, mean, var
	openfile.close()
	return
def norm_neut_ihs(inputScoreFile, outfileName, runProgram = "scans.py"):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_norm_ihs " + inputScoreFile + " > " + outfileName
	print(cmdStr)
	return
def norm_sel_ihs(inputScoreFile, neutNormfilename):
	''' from normalize_Ihs_manually() in func_scores.py''' 
	print("normalizing selection simulates to neutral: \n" + inputScoreFile + "\n" + neutNormfilename)
	bins, nums, means, variances = read_neut_normfile(neutNormfilename, 'ihs')
	#print(str(means))
	#print(str(variances))
	#print(str(bin_bounds))
	nsnps = 0
	normfilename = inputScoreFile + ".norm"	
	openfile = open(inputScoreFile, 'r')
	normfile = open(normfilename, 'w')
	for line in openfile:
		entries = line.split()
		freq_allele1 = float(entries[2]) #1 is ancestral, 0 is derived in tped
		unnormed_ihs_val = float(entries[5]) #locus/phys-pos/1_freq/ihh_1/ihh_0/ihs/derived_ihh_left/derived_ihh_right/ancestral_ihh_left/ancestral_ihh_right
		for ibin in range(len(bins)):
			if freq_allele1 <= bins[ibin]:
				normalizedvalue = (unnormed_ihs_val - means[ibin])/sqrt(variances[ibin])
				assert not(np.isnan(normalizedvalue))
				#if np.isnan(normalizedvalue):
					#print(freq_allele1)
					#print("\t" + str(unnormed_ihs_val) +"\t-\t" + str(means[ibin]) +"\t\\" + str(sqrt(variances[ibin])))
				break
		writeline = line.strip('\n') +"\t" + str(normalizedvalue) + '\n'
		normfile.write(writeline)		
	openfile.close()
	normfile.close()
	print("wrote to: " + normfilename)
	return
def norm_neut_xpehh(inputScoreFile, outfileName, runProgram = "scans.py"):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_norm_xpehh " + inputScoreFile + " > " + outfileName
	print(cmdStr)
	return
def norm_sel_xpehh(inputScoreFile, neutNormfilename):
	''' from normalize_Xp_manually in func_scores.py'''
	num, mean, var = read_neut_normfile(neutNormfilename, scoretype = 'xp')
	normfilename = inputScoreFile + ".norm"
	openfile = open(inputScoreFile, 'r')
	normfile = open(normfilename, 'w')
	header = openfile.readline()
	normfile.write(header)
	#ADJUST EACH VALUE 
	for line in openfile:
		entries = line.split()
		unnormed_xpehh = float(entries[-1])
		normed_xpehh = (unnormed_xpehh - mean)/sqrt(var)
		writeline = line.strip('\n') + "\t" + str(normed_xpehh) + "\n"
		normfile.write(writeline)
	openfile.close()
	normfile.close()
	print("wrote to: " + normfilename)
	return

##################
## HISTOGRAMS ###
###################
def get_hist_bins(score,numBins):
	"""toggle cutoffs of histogram, as well as y-axis limits"""
	if score == "ihs":
		scorerange = [-4., 4.]
		ylims = [0, .25]
	elif score == "delihh":
		scorerange = [-2., 3.]
		ylims = [0, .25]
	elif score == "fst":
		scorerange = [-.05, 1.]#[-2., 4.5]
		ylims = [0, .25]
	elif score == "deldaf":
		scorerange = [-1., 1.]
		ylims = [0, .25]
	elif score == "xp":
		scorerange = [-3., 7.]
		ylims = [0, .25]
	binlen = (scorerange[1] - scorerange[0])/float(numBins-1)
	bins = [scorerange[0] + binlen * i for i in range(numBins)]
	return bins, scorerange, ylims
def get_indices(score, dem_scenario):
	"""CRUCIAL FUNC; tells loadVals which columns to grab and return; changes with filetype. also returns anc_freq_index"""
	if score == "ihs":
		physpos_index, freq_anc_index = 1, 2
		if "sel" in dem_scenario:
			ihs_unnormed_index, ihs_normed_index, expectedlen = 9, 10, 11 #these files have 11 columns
		else: #neutral
			ihs_unnormed_index, ihs_normed_index, expectedlen = 5, 6, 8 #these files have 8 columns; last one is binary variable
		indices = [physpos_index, ihs_normed_index, freq_anc_index]#freq_anc_index, ihs_unnormed_index, ihs_normed_index] 
	elif score == "delihh":
		expectedlen = 8
		physpos_index, freq_anc_index = 1, 2
		delihh_unnormed_index, delihh_normed_index = 5, 6
		indices = [physpos_index, delihh_normed_index, freq_anc_index]#freq_anc_index, delihh_unnormed_index, delihh_normed_index]
	elif score == "fst": #must validate******
		expectedlen = 3#2
		physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "deldaf":
		expectedlen = 3#2
		physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "xp":
		physpos_index, freq_anc_index = 1, 2
		xp_unnormed_index, xp_normed_index, = 7, 8
		if "sel" in dem_scenario:
			expectedlen = 9
		else:
			expectedlen = 10
		indices = [physpos_index, xp_normed_index, freq_anc_index]
	return expectedlen, indices
def load_vals_from_files(filename, numCols, takeindices, stripHeader = False):
	''' from func_clean.py '''
	toreturn, incompleteData = [[] for index in takeindices], 0
	openfile = open(filename, 'r')
	if stripHeader:
		header = openfile.readline()
	for line in openfile:
		entries = line.split()
		if len(entries) != numCols:
			print("ERROR: numCols " + str(numCols) + " " + str(len(entries)) + " " + filename)
			incompleteData +=1
		for iIndex in range(len(takeindices)):
			index = takeindices[iIndex]
			thisValue = float(entries[index])
			toreturn[iIndex].append(thisValue)
	openfile.close()
	return toreturn
def calc_hist_from_scores(causal_scores, linked_scores, neut_scores, xlims, thinToSize = False):
	if thinToSize:
		limiting = min(len(causal_scores), len(linked_scores), len(neut_scores))
		print("Thinning data to " + str(limiting) + " SNPs...")
		short_causal, short_linked, short_neut = [], [], []

		for parentlist in ['causal', 'linked', 'neut']:
			shortlist = eval('short_' + parentlist)
			takenindices = []
			while len(takenindices) < limiting:
				chosenIndex = randint(0, limiting-1)
				if chosenIndex not in takenindices:
					shortlist.append(eval(parentlist + "_scores")[chosenIndex])
					takenindices.append(chosenIndex)

		causal_scores, linked_scores, neut_scores = short_causal, short_linked, short_neut

	#get weights to plot pdf from hist
	weights_causal = np.ones_like(causal_scores)/len(causal_scores)
	weights_linked = np.ones_like(linked_scores)/len(linked_scores)
	weights_neut = np.ones_like(neut_scores)/len(neut_scores)

	causal_scores = np.clip(causal_scores, xlims[0], xlims[1]) #np.clip
	linked_scores = np.clip(linked_scores, xlims[0], xlims[1])
	neut_scores = np.clip(neut_scores, xlims[0], xlims[1])

	n_causal, bins_causal = np.hist(causal_scores, bins=givenBins, weights = weights_causal)
	n_linked, bins_linked = np.hist(linked_scores, bins=givenBins, weights = weights_linked)
	n_neut, bins_neut = np.hist(neut_scores, bins=givenBins, weights = weights_neut)

	return n_causal, n_linked, n_neut, bin_causal, bins_linked, bins_neut
def write_hists_to_files(writePrefix, givenBins, n_causal, n_linked, n_neut):
	for status in ['causal', 'linked', 'neutral']:
		writefilename = writePrefix + "_" + status + ".txt"
		writefile = open(writefilename, 'w')
		n_scores = eval('n_' + status)
		for index in range(len(n_scores)):
			towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(n_scores[index])+ "\n"
			writefile.write(towritestring)
		writefile.close()
	return