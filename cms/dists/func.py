## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.03.16 vitti@broadinstitute.org

import subprocess
import sys
import os

##################
### BOOKKEEPING ##
##################
def give_command(commandstr):
	commands = commandstr.split(' ')
	subprocess.check_output(commands)
def check_make_dir(dirpath):
	"""function to handle creation of (sub)folders, avoiding duplication.
	"""
	if os.path.exists(dirpath):
		contents = os.listdir(dirpath)
		if len(contents) == 0:
			return dirpath
		else:
			alt, iAlt = False, 1
			while alt == False:			
				altpath = dirpath.strip('/') + "_alt" + str(iAlt) #best way to avoid recursive issues?
				if os.path.exists(altpath):
					iAlt +=1
				else:
					alt = True
			mkdircommand = "mkdir " + altpath
			give_command(mkdircommand)
			return altpath
	else:
		mkdircommand = "mkdir " + dirpath
		give_command(mkdircommand)
	return dirpath
def check_file_len(filename):
	'''counts number of lines in file'''
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		openfile = open(filename, 'r')
		iline = 0
		for line in openfile:
			iline +=1
		openfile.close()
		return iline
	else:
		return 0

###################
### SELFREQ BINS ##
###################
def get_bin_strings(bin_medians):
	'''given input bin medians, returns minimal len strs needed to create unique bin labels'''
	stringlen = 3 #minimum is '0.x' -- not necessarily the most descriptive way to round these, but it's functional.
	unique = False
	while unique == False:
		labels = [str(x)[:stringlen] for x in bin_medians]
		if len(labels) == len(set(labels)):
			unique = True
		else:
			stringlen +=1
	return labels
def get_bins(freqRange=".05-.95", numbins=9):
	fullrange = [float (x) for x in freqRange.split('-')]
	binlen = (fullrange[1] - fullrange[0]) / float(numbins)
	bin_starts = [fullrange[0] + binlen*i for i in range(numbins)]
	bin_ends = [fullrange[0] + binlen*i for i in range(1,numbins+1)]
	bin_medians = [(float(bin_starts[i]) + float(bin_ends[i]))/2. for i in range(len(bin_starts))]
	bin_medians_str = get_bin_strings(bin_medians)
	return fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str	
def check_bin_filled(directory, numsims):
	'''counts all non-empty files in directory, returns True if >= numsims'''
	filenames = os.listdir(directory)
	nonempty = [x for x in filenames if check_file_len(x) !=0]
	if len(filenames) < numsims:
		return False
	else:
		return True

#####################
## CALCS FROM SIMS ## ## gratuitous. should probably nix (?)
#####################
def calc_ihs(inputTped, outputFile, runProgram = "scans.py", numThreads = 7):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_ihs " + inputTped + " " + outputFile + " --threads " + str(numThreads)
	print cmdStr
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
def calc_xp(inputTped, inputTped2, outputFile, runProgram = "scans.py", numThreads = 7):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_xpehh " + inputTped + " " + outputFile + " " + inputTped2 + " --threads " + str(numThreads)
	print cmdStr
	return	
def calc_fst_deldaf(inputTped, inputTped2, outputFile):
	print "assumes C program compiled "
	commandstring = "./calc_fst_deldaf_tped"
	argstring = inputTped1 + " " + inputTped2 + " " + outputFile
	print commandstring + " " + argstring
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
			entries = line.split()
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
		#print str(num)
		#print str(nsnp)
		#assert num == nsnp
		return num, mean, var
	openfile.close()
	return
def norm_neut_ihs(inputScoreFile, outfileName, runProgram = "scans.py"):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_norm_ihs " + inputScoreFile + " > " + outfileName
	print cmdStr
	return
def norm_sel_ihs(inputScoreFile, neutNormfilename, bin_bounds):
	''' from normalize_Ihs_manually() in func_scores.py''' 
	print "normalizing selection simulates to neutral: " 
	print inputScoreFile
	print neutNormfilename
	bins, nums, means, variances = read_neut_normfile(neutNormfilename, 'ihs')
	nsnps = 0
	normfilename = inputScoreFile + ".norm"	
	openfile = open(writefilename, 'r')
	normfile = open(normfilename, 'w')
	for line in openfile:
		entries = line.split()
		freq_allele1 = float(entries[2]) #1 is ancestral, 0 is derived in tped
 		unnormed_ihs_val = float(entries[5]) #locus           phys-pos        1_freq          ihh_1           ihh_0           ihs             derived_ihh_left     derived_ihh_right    ancestral_ihh_left      ancestral_ihh_right
		for ibin in range(len(bin_bounds)):
			if freq_allele1 <= bin_bounds[ibin]:
				normalizedvalue = (unnormed_ihs_val - means_bin[ibin])/sqrt(vars_bin[ibin])
				break
		writeline = line.strip('\n') +"\t" + str(normalizedvalue) + '\n'
		normfile.write(writeline)		
	openfile.close()
	normfile.close()
	print "wrote to: " + normfilename
	retur
def norm_neut_xpehh(inputScoreFile, outfileName, runProgram = "scans.py"):
	'''from func_clean.py'''
	cmdStr = "python " + runProgram + " selscan_norm_xpehh " + inputScoreFile + " > " + outfileName
	print cmdStr
	return
def norm_sel_xpehh(inputScoreFile, neutNormfilename, bin_bounds):
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
	print "wrote to: " + normfilename
	retur

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
	"""CRUCIAL FUNC; tells loadVals which columns to grab and return; changes with filetype"""
	if score == "ihs":
		physpos_index, freq_anc_index = 1, 2
		if "sel" in dem_scenario:
			ihs_unnormed_index, ihs_normed_index, expectedlen = 9, 10, 11 #these files have 11 columns
		else: #neutral
			ihs_unnormed_index, ihs_normed_index, expectedlen = 5, 6, 8 #these files have 8 columns; last one is binary variable
		indices = [physpos_index, ihs_normed_index]#freq_anc_index, ihs_unnormed_index, ihs_normed_index] 
	elif score == "delihh":
		expectedlen = 8
		physpos_index, freq_anc_index = 1, 2
		delihh_unnormed_index, delihh_normed_index = 5, 6
		indices = [physpos_index, delihh_normed_index]#freq_anc_index, delihh_unnormed_index, delihh_normed_index]
	elif score == "fst":
		expectedlen = 2
		physpos_index, scoreindex = 0, 1
		indices = [physpos_index, scoreindex]
	elif score == "deldaf":
		expectedlen = 2
		physpos_index, scoreindex = 0, 1
		indices = [physpos_index, scoreindex]
	elif score == "xp":
		physpos_index, freq_anc_index = 1, 2
		xp_unnormed_index, xp_normed_index, = 7, 8
		if "sel" in dem_scenario:
			expectedlen = 9
		else:
			expectedlen = 10
		indices = [physpos_index, xp_normed_index]
		#[physpos_index, freq_anc_index, xp_unnormed_index, xp_normed_index]
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
			print "ERROR: numCols " + str(numCols) + " " + str(len(entries)) + " " + filename
			incompleteData +=1
		for iIndex in range(len(takeindices)):
			index = takeindices[iIndex]
			thisValue = float(entries[index])
			toreturn[iIndex].append(thisValue)
	openfile.close()
	return toreturn