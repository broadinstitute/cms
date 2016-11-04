## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 11.4.16 vitti@broadinstitute.org

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
	"""taken from JV func_scores.py. given a selscan iHS file, parses it and writes delihh file"""
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		entries = line.split()
		#handle input with/without ihh decomp
		if len(entries) == 8:
			locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, ihs_normed, lastcol = entries
		elif len(entries) == 10:
			locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, der_ihh_l, der_ihh_r, anc_ihh_l, anc_ihh_r  = entries
		elif len(entries) == 11:
			locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, der_ihh_l, der_ihh_r, anc_ihh_l, anc_ihh_r, manually_normed = entries
					#ancestral - derived
		else:
			locus, phys, freq_1, ihh_1, ihh_0, ihs_unnormed, ihs_normed, lastcol = entries
		unstand_delIHH = fabs(float(ihh_1) - float(ihh_0))
		writeline = locus + "\t" + phys + "\t" + freq_1 + "\t" + str(ihs_unnormed) + "\t" + str(unstand_delIHH) +"\t" + str(unstand_delIHH) +  "\n" #6 columns for selscan norm
		writefile.write(writeline)
	writefile.close()
	#print("wrote to " + writefilename)
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
def norm_neut_ihs(inputScoreFile, outfileName, runProgram = "cms/cms/scans.py"):
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
		normalizedvalue = float('NaN')
		for ibin in range(len(bins)):
			if freq_allele1 <= bins[ibin]:
				normalizedvalue = (unnormed_ihs_val - means[ibin])/sqrt(variances[ibin])
				#assert not(np.isnan(normalizedvalue))
				#if np.isnan(normalizedvalue):
				#	print(freq_allele1)
				#	print("\t" + str(unnormed_ihs_val) +"\t-\t" + str(means[ibin]) +"\t\\" + str(sqrt(variances[ibin])))
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
		if "sel" in dem_scenario:
			expectedlen = 7
		else:
			expectedlen = 8
		physpos_index, freq_anc_index = 1, 2
		delihh_unnormed_index, delihh_normed_index = 5, 6
		indices = [physpos_index, delihh_normed_index, freq_anc_index]#freq_anc_index, delihh_unnormed_index, delihh_normed_index]
	elif score == "nsl":
		physpos_index, freq_anc_index = 1, 2
		normedscoreindex = 6
		if "sel" in dem_scenario:
			expectedlen = 7
		else:
			expectedlen = 8
		indices = [physpos_index, normedscoreindex, freq_anc_index]		
	elif score == "fst": #PULL FROM COMP (.CMS FILE)
		expectedlen = 8
		physpos_index, scoreindex, freq_anc_index = 0, 5, float('nan')
		#expectedlen = 3#2
		#physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "deldaf":#PULL FROM COMP (.CMS FILE)
		expectedlen = 8
		physpos_index, scoreindex, freq_anc_index = 0, 6, float('nan')	
		#expectedlen = 3#2
		#physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "xp": #PULL FROM COMP (.CMS FILE)
		expectedlen = 8
		physpos_index, scoreindex, freq_anc_index = 0, 4, float('nan')	
		#physpos_index, freq_anc_index = 1, 2
		#xp_unnormed_index, xp_normed_index, = 7, 8
		#if "sel" in dem_scenario:
		#	expectedlen = 9
		#else:
		#	expectedlen = 10
		indices = [physpos_index, xp_normed_index, freq_anc_index]

	return expectedlen, indices
def load_vals_from_files(filename, numCols, takeindices, stripHeader = False, printProgress = False, checkCols = False):
	''' if filename is .list, opens and parses multiple files '''

	toreturn, incompleteData = [[] for index in takeindices], 0

	entries = filename.split('.')
	if entries[-1] == "list":
		allfilenames = []
		openfile = open(filename, 'r')
		for line in openfile:
			repfilename = line.strip('\n')
			if os.path.isfile(repfilename):
				allfilenames.append(repfilename)
		openfile.close()

	else:
		allfilenames = [filename]
	print('loading data from ' + str(len(allfilenames)) + ' files...')

	for ifilename in range(len(allfilenames)):
		filename = allfilenames[ifilename]
		if printProgress:
			if ifilename % 100 == 0: 
				print("now file " + str(ifilename) + " out of " + str(len(allfilenames)))
		openfile = open(filename, 'r')
		if stripHeader:
			header = openfile.readline()
		for line in openfile:
			entries = line.split()
			if entries[0] != "chrom": #quick patch for calc_fst_deldaf printing chrom-wide average
				if checkCols and (len(entries) != numCols):
					print("ERROR: numCols " + str(numCols) + " " + str(len(entries)) + " " + filename)
					incompleteData +=1
					break
				if takeindices[-1] ==float('nan'):
					fullrange = len(takeindices) - 1
				else:
					fullrange = len(takeindices)
				for iIndex in range(fullrange):
					index = takeindices[iIndex]
					thisValue = float(entries[index])
					toreturn[iIndex].append(thisValue)
		openfile.close()
	return toreturn
def choose_vals_from_files(filename, numCols, takeindices, comp, stripHeader = False, checkCols = False, method = "max"):
	''' expects a .list file with multiple records (i.e., same replicate, different poppairs) on the same line as input '''

	entries = filename.split('.')
	if entries[-1] == "list":
		allfilenames = []
		openfile = open(filename, 'r')
		for line in openfile:
			thisrepfiles = []
			repfilenames = line.split()
			for repfilename in repfilenames:
				if os.path.isfile(repfilename):
					thisrepfiles.append(repfilename)
			allfilenames.append(thisrepfiles)
		openfile.close()

	else:
		print('check input file.')
		sys.exit(0)
		#allfilenames = [filename]
	print('loading data from ' + str(len(allfilenames)) + ' replicates...')
	
	#assert len(takeindices) == 3
	#print(takeindices)
	#sys.exit()

	posIndex, takeIndex, ancfreqIndex = takeindices
	#if comp == "fst":
	#	posIndex, takeIndex = 0, 1
	#elif comp == "deldaf":
	#	posIndex, takeIndex = 0, 2
	#else:
	#	posIndex, takeIndex = takeindices[0], takeindices[1]
	
	allpos, allscores = [], []

	for ireplicate in range(len(allfilenames)):
		filelist = allfilenames[ireplicate]
		repscores = []
		reppositions = []
		repanc = []
		for filename in filelist:
			#print(filename)
			#print(str(takeIndex))
			vals = []
			positions =[]
			ancs = []
			openfile = open(filename, 'r')
			if stripHeader:
				header = openfile.readline()
			for line in openfile:
				entries = line.split()
				if entries[0] != "chrom": #quick patch for calc_fst_deldaf printing chrom-wide average
					if checkCols and (len(entries) != numCols):
						print("ERROR: numCols " + str(numCols) + " " + str(len(entries)) + " " + filename)
						incompleteData +=1
						break
					val = float(entries[takeIndex])
					pos = int(entries[posIndex])
					anc = float(entries[ancfreqIndex])
					vals.append(val)
					positions.append(pos)
					ancs.append(anc)
			openfile.close()
			repscores.append(vals)
			reppositions.append(positions)
			repanc.append(ancs)
		#choose here
		chosen, chosenpos = choose_from_reps(repscores, reppositions,  repanc, mode = method)
		allscores.extend(chosen)
		allpos.extend(chosenpos)
	alltoreturn = [allpos, allscores]
	return alltoreturn

def choose_from_reps(repscores, reppositions, repanc, mode="max"):
	'''flexible function to choose for likes '''
	ncomp = len(repscores)
	allpositions = []
	for reppositionlist in reppositions:
		allpositions.extend(reppositionlist)
	allpositions = set(allpositions)
	toreturn = []
	for position in allpositions:
		scores = []
		freqs = []
		for icomp in range(ncomp):
			if position in reppositions[icomp]:
				thisindex = reppositions[icomp].index(position)
				scores.append(repscores[icomp][thisindex])		
				freqs.append(repanc[icomp][thisindex])
		#this is where choosing happens
	
		if mode == "max":
			itemtoreturn = max(scores)
		elif mode in ['mean', 'ave', 'average']:
			itemtoreturn = np.mean(scores)
		elif mode == "min":
			itemtoreturn = min(scores)
		elif mode == "daf" or mode == "deldaf": #
			itemtoreturn = np.mean(scores)
			#print('testing...')
			#thispop_anc = float(freqs[0])
			#print(thispop_anc)
			#thispop_der = 1. - thispop_anc
			#otherpops_anc = freqs[1:]
			#otherpops_der = [1. - float(item) for item in otherpops_anc]
			#otherpops_ave = np.mean(otherpops_der)
			#daf_val = thispop_der - otherpops_ave
			#itemtoreturn = daf_val
			#print(freqs)
			#print(daf_val)

		toreturn.append(itemtoreturn)
		#print(scores)
		#print(itemtoreturn)
	return toreturn, allpositions


def calc_hist_from_scores(causal_scores, linked_scores, neut_scores, xlims, givenBins, thinToSize = False):
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

	#print(causal_scores)

	#check for nans
	causal_scores = np.array(causal_scores)
	linked_scores = np.array(linked_scores)
	neut_scores = np.array(neut_scores)
	causal_scores = causal_scores[~np.isnan(causal_scores)]
	linked_scores = linked_scores[~np.isnan(linked_scores)]
	neut_scores = neut_scores[~np.isnan(neut_scores)]

	#get weights to plot pdf from hist
	weights_causal = np.ones_like(causal_scores)/len(causal_scores)
	weights_linked = np.ones_like(linked_scores)/len(linked_scores)
	weights_neut = np.ones_like(neut_scores)/len(neut_scores)

	causal_scores = np.clip(causal_scores, xlims[0], xlims[1]) #np.clip
	linked_scores = np.clip(linked_scores, xlims[0], xlims[1])
	neut_scores = np.clip(neut_scores, xlims[0], xlims[1])

	n_causal, bins_causal = np.histogram(causal_scores, range=xlims, bins=givenBins, weights = weights_causal)
	n_linked, bins_linked = np.histogram(linked_scores, range=xlims,  bins=givenBins, weights = weights_linked)
	n_neut, bins_neut = np.histogram(neut_scores,range=xlims, bins=givenBins, weights = weights_neut)


	#debug_array = [n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut]
	#print(debug_array)

	totalNsnps = len(causal_scores) + len(linked_scores) + len(neut_scores)

	#add pseudocount (pseudoprob) for empty bins
	for ibin in range(len(n_causal)):
		if n_causal[ibin] == 0:
			n_causal[ibin] = (1./totalNsnps)
		if n_linked[ibin] == 0:
			n_linked[ibin] = (1./totalNsnps)
		if n_neut[ibin] == 0:
			n_neut[ibin] = (1./totalNsnps)

	return n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut


def write_hists_to_files(writePrefix, givenBins, n_causal, n_linked, n_neut):
	assert len(givenBins) == (len(n_causal) + 1)
	for status in ['causal', 'linked', 'neut']:
		writefilename = writePrefix + "_" + status + ".txt"
		writefile = open(writefilename, 'w')
		n_scores = eval('n_' + status)
		for index in range(len(n_scores)):
			towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(n_scores[index])+ "\n"
			writefile.write(towritestring)
		writefile.close()
	return