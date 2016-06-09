## {{DEM MODEL --> COMPONENT SCORE PROBABILITY DISTS}} -- ancillary functions.
## last updated 05.31.16 	vitti@broadinstitute.org

import os
import subprocess
import numpy as np
import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as plt
import pylab
from math import sqrt
from random import randint


#################
## BOOKKEEPING ##
#################
def givecommand(commandstr):
	commands = commandstr.split(' ')
	subprocess.check_output(commands)
def checkMakeDir(dirpath, overwrite = False, createAlt = False, verbose = False):
	"""flexible function to handle creation of (sub)folders, avoiding duplication as specified.
	by default: will create DIRPATH if it doesn't exist; else do nothing
	"""
	if os.path.exists(dirpath):
		contents = os.listdir(dirpath)
		if len(contents) == 0:
			if verbose:
				print dirpath + " already exists and is empty"
			return
		else:
			if verbose:
				print dirpath + " already exists and has " + str(len(contents)) +  " files"
			if overwrite == True:
				if verbose:
					print "overwriting..."
				rmdircommand = "rm -R " + dirpath
				givecommand(rmdircommand)
				mkdircommand = "mkdir " + dirpath
				givecommand(mkdircommand)
			elif createAlt == True:
				altpath = dirpath + "_alt" + str(randint(0, 99999)) #best way to avoid recursive issues?
				if verbose:
					print "writing to alternate folder: " + altpath
				mkdircommand = "mkdir " + altpath
				givecommand(mkdircommand)
				return
	else:
		if verbose:
			print dirpath + " does not exist; creating..."
		mkdircommand = "mkdir " + dirpath
		givecommand(mkdircommand)
	return
def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-q long', '-P sabeti_lab'], dependencies = ['.python-2.7.8-sqlite3-rtrees', '.zlib-1.2.6', '.gsl-1.14', '.numpy-1.8.0-python-2.7.1-sqlite3-rtrees', '.samtools-1.2', '.tabix-0.2.6', '.biopython-1.64-python-2.7.1-sqlite3-rtrees', '.gcc-5.2.0', 'GCC-5.1', '.glib-2.42.2']):
	"""flexible function for generating scripts and argfiles to dispatch parallel jobs on Univa Grid Engine.
	arglist: a list of strings. each string represents one argstring, for one task on the array."""
	#write file with arguments for all tasks
	t = 0
	taskargfilename = filewd + fileprefix + "_taskargs.txt"
	taskargfile = open(taskargfilename, "w")
	for argstring in arguments:
		taskargfile.write(argstring + "\n")
		t += 1
	taskargfile.close() 

	#write script for one job that points to task-argument-file
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("#$ -t 1-" + str(t) + "\n")
	scriptfile.write("\n")
	#scriptfile.write('eval `/broad/software/dotkit/init -b`\n')	
	scriptfile.write('source /broad/software/scripts/useuse\n')
	for dependency in dependencies:
		scriptfile.write('reuse -q ' + dependency + '\n')
	scriptfile.write('source /idi/sabeti-scratch/jvitti/cms_venv/bin/activate\n')
	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	subcommand = "qsub " + scriptname
	print subcommand
def checkFileLen(filename, minNumLines = 1, warn = True):
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		openfile = open(filename, 'r')
		iline = 0
		for line in openfile:
			iline +=1
			if iline > minNumLines:
				return filename
		openfile.close()
		if warn == True:
			print "Check possible file truncation: " + filename + "\t" + str(iline) 
		return filename + "\t" + str(iline) 
	else:
		#print "Could not find: " + filename
		return "MISSING"

#################
## SIMULATIONS ##
#################
def writeSweepToParamFile(neutralparamfilename, sweepline, sweeplabel):
	"""input: a parameter file of a neutral model, and variables for a selective sweep in coalescent_msweep format
	output: a new parameter file with the sweep added in."""
	writefilename = neutralparamfilename[:-4] + "_" + sweeplabel + ".par"
	writefile = open(writefilename, 'w')
	openfile = open(neutralparamfilename, 'r')
	for line in openfile:
		writefile.write(line)
	openfile.close()
	writefile.write(sweepline + "\n") 
	writefile.close()
	return writefilename
def getDictFromParamFile(paramfilename):
	#initialize empty lists so we can append from paramfile as we parse
	paramDict = {'numPops':4,'labels':{1:'YRI', 2:'CEU', 3:'CHB', 4:'BEB'}, 'singrate':.25, 'presentSizes':[], 'num_indivs_per_sample':200}
	for pop in range(1,5):
		for parameter in ['Ne', 'Ne_times', 'Ne_change', 'bn', 'bn_times']:
			key = (parameter, pop)
			paramDict[key] = []
		for jpop in range(1,5):
			if pop != jpop:
				migstring = str(pop) + "->" + str(jpop)
				key1, key2 = ('mig', migstring), ('mig_time', migstring)
				paramDict[key1] = []
				paramDict[key2] = []
	openfile = open(paramfilename, 'r')
	for line in openfile:
		if line[0] != "#":
			entries = line.split()
			if "length" in entries:
				chromlength = int(entries[-1])
				paramDict['chromlength'] = chromlength
			elif "mutation_rate" in entries:
				mutation_rate = float(entries[-1])
				paramDict['mutation_rate'] = mutation_rate
			elif 'recomb_file' in entries:
				recomb_file = entries[-1]
				paramDict['recomfilename'] = recomb_file
			elif "gene_conversion_relative_rate" in entries:
				gene_conversion_relative_rate = float(entries[-1])
				paramDict['gene_conv_rel_rate'] = gene_conversion_relative_rate
			elif "pop_size" in entries:
				currentsize = int(entries[2])
				paramDict['presentSizes'].append(currentsize)
				
			elif "migration_rate" in entries:
				sourcepop, targetpop, time, rate = int(entries[-4]), int(entries[-3]), float(entries[-2]), float(entries[-1])
				migstring = str(sourcepop) + "->" + str(targetpop)
				key1 = ('mig', migstring)
				key2 = ('mig_time', migstring)
				paramDict[key1].append(rate)
				paramDict[key2].append(time)
			elif "bottleneck" in entries:
				pop, time, rate = int(entries[-3]), float(entries[-2]), float(entries[-1])
				key1 = ('bn', pop)
				key2 = ('bn_times', pop)
				paramDict[key1].append(rate)
				paramDict[key2].append(time)
			elif "exp_change_size" in entries or "exp_change_size2" in entries:
				pop = int(entries[-5])
				finaltime, starttime = float(entries[-4]), float(entries[-3]) #forward-looking concept of 'final/start'
				finalsize, startsize = float(entries[-2]), float(entries[-1])
				key1 = ('Ne', pop)
				key2 = ('Ne_times', pop)
				key3 = ('Ne_change', pop)
				paramDict[key1].append(finalsize)
				paramDict[key2].append(finaltime)
				paramDict[key3].append(1)
			elif "change_size" in entries:
				pop, time, size = int(entries[-3]), float(entries[-2]), float(entries[-1])
				key1 = ('Ne', pop)
				key2 = ('Ne_times', pop)
				key3 = ('Ne_change', pop)
				paramDict[key1].append(size)
				paramDict[key2].append(time)
				paramDict[key3].append(0)
			elif "split" in entries:
				parentpop, childpop, time = int(entries[-3]), int(entries[-2]), float(entries[-1])
				assert parentpop == childpop - 1
				key = ('split', childpop)
				paramDict[key] = [time]
			elif "admix" in entries:
				targetpop, sourcepop, time, rate = int(entries[-4]), int(entries[-3]), float(entries[-2]), float(entries[-1])
				paramDict['admix_rate'] = [rate]
				paramDict['admix_time'] = [time]
	openfile.close()
	#print paramDict
	return paramDict
def rewriteDistributionString(diststring):
	"""gets rid of spaces so that we can use these variables in our filenaming"""
	toReturn = ""
	entries = diststring.split()
	for item in entries:
		for character in item:
			if character.isalnum():
				toReturn += character 
			elif character == ".":
				toReturn += "-"
		toReturn += "_"
	return toReturn[:-1]
def getFreqTrajs(trajfilename):
	openfile = open(trajfilename, 'r')
	openfile.readline() #header
	all_gens, all_f1s, all_f2s, all_f3s, all_f4s = [], [], [], [], [] #store by replicate
	gens, f1s, f2s, f3s, f4s = [], [], [], [], [] #for first replicate
	lastrep = 1
	for line in openfile:
		rep, gen, f1, f2, f3, f4 = line.split()
		rep = int(rep)
		if rep != lastrep: #past final generation of replicate 
			all_gens.append(gens)
			all_f1s.append(f1s)
			all_f2s.append(f2s)
			all_f3s.append(f3s)
			all_f4s.append(f4s)
			gens, f1s, f2s, f3s, f4s = [], [], [], [], []
		gens.append(float(gen))
		f1s.append(float(f1))
		f2s.append(float(f2))
		f3s.append(float(f3))	
		f4s.append(float(f4))					
		lastrep = rep
	#final replicate
	all_gens.append(gens)
	all_f1s.append(f1s)
	all_f2s.append(f2s)
	all_f3s.append(f3s)
	all_f4s.append(f4s)
	openfile.close()
	return all_gens, all_f1s, all_f2s, all_f3s, all_f4s
def getFinalSelFreq(trajfilename, selpop = 1):
	"""Assumes a singular trajectory in the file!"""
	all_gens, all_f1s, all_f2s, all_f3s, all_f4s = getFreqTrajs(trajfilename)
	return eval('all_f' + str(selpop)+'s')[0][-1]

##################
## CALCULATIONS ##
##################
def launch_calc_ihs_batch(scenario, inputs, outputs, numThreads = 35, runProgram = "/idi/sabeti-scratch/jvitti/cms_venv/cms/cms/scans.py"):
	taskargstrings = []
	for i in range(len(inputs)):
		inputTped = inputs[i]
		outputFile = outputs[i]
		argstring = "selscan_ihs " + inputTped + " " + outputFile + " --threads " + str(numThreads)
		taskargstrings.append(argstring)
	commandstring = "python2.7 " + runProgram
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_ihs_" + scenario)	
def normalizeiHS(scenario, nRep, ihsloc, pops = [1, 2, 3, 4], runProgram = "/idi/sabeti-scratch/jvitti/cms_venv/cms/cms/scans.py"):
	taskargstrings = []
	for pop in pops:
		nsnps = 0
		writefilename = ihsloc  + "pop" + str(pop) + ".ihs.concat"
		writefile = open(writefilename, 'w')
		for iRep in range(nRep):
			ihsfile = ihsloc + "rep" + str(iRep) + "_pop" + str(pop) + ".ihs.out"
			if os.path.isfile(ihsfile):
				readfile = open(ihsfile, 'r')
				for line in readfile:
					writefile.write(line)
					nsnps +=1
				readfile.close()
		writefile.close()
		print "for pop " + str(pop) + " we logged " + str(nsnps) + " scores to " + writefilename
		taskargstring = "selscan_norm_ihs " + writefilename 
		taskargstrings.append(taskargstring)
	commandstring = "python2.7 " + runProgram
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "norm_ihs_" + scenario)		
	return	
def readNeutNormFile(neutnormfilename, scoretype ='xp'):
	"""pulls means and vars so that we can normalize sel scenarios with the same parameters"""
	openfile = open(neutnormfilename)
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
def normalize_Ihs_manually(ihsdir, nRep, pop, means_bin, vars_bin, bin_bounds):
	"""allows us to pass the bin means and vars from neutral simulations to normalize selection sims accordingly
	means, vars: array of values
	"""
	nsnps = 0
	writefilename = ihsdir + "n" + str(nRep) + ".ihs.concat"
	writefile = open(writefilename, 'w')
	for iRep in range(nRep):
		ihsfile = ihsdir + "/rep" + str(iRep) + "_pop" + str(pop) + ".ihs.out"
		if os.path.isfile(ihsfile):
			readfile = open(ihsfile, 'r')
			for line in readfile:
				writefile.write(line)
				nsnps +=1
			readfile.close()
	writefile.close()
	print "wrote to: " + writefilename
	normfilename = writefilename + ".norm"	
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

def launch_calc_xp_batch(scenario, inputs_1, inputs_2, outputs, numThreads = 35, runProgram = "/idi/sabeti-scratch/jvitti/cms_venv/cms/cms/scans.py"):
	taskargstrings = []
	for i in range(len(inputs_1)):
		inputTped1 = inputs_1[i]
		inputTped2 = inputs_2[i]
		outputFile = outputs[i]
		argstring = "selscan_xpehh " + inputTped1 + " " + outputFile +  " " + inputTped2 + " --threads " + str(numThreads)
		taskargstrings.append(argstring)
	commandstring = "python2.7 " + runProgram
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_xp_" + scenario)	

#################
## LIKELIHOODS ##
#################
def getNeutralFilename(score, dem_scenario, selpop, altpop = ''):
	neutral_value_filename = ""	
	popCompString = "pop" + str(selpop) + "_pop" + str(altpop)
	if score == 'ihs':
		neutral_value_filename = "scores/cosi_ihs/" + dem_scenario + "/neut/pop" + str(selpop) + ".ihs.concat.100bins.norm" 
	elif score == "delihh":
		neutral_value_filename = "scores/cosi_delihh/" + dem_scenario + "/neut/pop" + str(selpop) + ".delihh.concat.100bins.norm" 
	elif score == "fst":
		neutral_value_filename = "scores/cosi_fst_deldaf/" + dem_scenario + "/neut/" + popCompString + "_fst.txt"
	elif score == "deldaf":
		neutral_value_filename = "scores/cosi_fst_deldaf/" + dem_scenario + "/neut/" + popCompString + "_deldaf.txt"
	elif score == "xp":
		neutral_value_filename = "scores/cosi_xp/" + dem_scenario + "/neut/" + popCompString +".xpehh.concat.norm"
	else:
		print "ERROR: scoretype " + str(score)
	if not os.path.isfile(neutral_value_filename):
		print "MISSING FILE:" + neutral_value_filename
	return neutral_value_filename
def getSelFilename(score, dem_scenario, selpop, selfreq, altpop = ''):
	#sel_value_filename = ""	
	popCompString = "pop" + str(selpop) + "_pop" + str(altpop)
	if score == 'ihs':
		sel_value_filename = "scores/cosi_ihs/" + dem_scenario + "/sel" + str(selpop) + "/sel" + str(selfreq) + "/n300.ihs.concat.norm"
	elif score == "delihh":
		sel_value_filename = "scores/cosi_delihh/" + dem_scenario + "/sel" + str(selpop) + "/sel" + str(selfreq) + "/n300.delihh.concat.100bins.norm"
	elif score == "fst":
		sel_value_filename = "scores/cosi_fst_deldaf/" + dem_scenario + "/sel" + str(selpop) + "/sel" + str(selfreq) + "/" + popCompString + "_fst.txt"
	elif score == "deldaf":
		sel_value_filename = "scores/cosi_fst_deldaf/" + dem_scenario + "/sel" + str(selpop) + "/sel" + str(selfreq) + "/" + popCompString + "_deldaf.txt"
	elif score == "xp":
		sel_value_filename = "scores/cosi_xp/" + dem_scenario+ "/sel" + str(selpop) + "/sel" + str(selfreq) + "/pop" + str(selpop) + "_pop" + str(altpop) +".xpehh.concat.norm"
	else:
		print "ERROR: scoretype " + str(score)
	if not os.path.isfile(sel_value_filename):
		print "MISSING FILE:" + sel_value_filename
	return sel_value_filename
def getIndices(score, dem_scenario):
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
def getBins(score, numBins):
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
def loadValsFromFile(filename, numCols, takeindices, stripHeader = False, verbose = False):
	toreturn, incompleteData = [[] for index in takeindices], 0
	if verbose:
		print "~loading vals from file~"
		print filename
	openfile = open(filename, 'r')
	if stripHeader:
		header = openfile.readline()
		#if verbose:
		#	print "strip header: "
		#	print header
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
	if verbose:
		print "missing data for " + str(incompleteData) + " snps..."
		print "loaded " + str(len(toreturn[0])) + " datapoints..."
	return toreturn
def plotFig(ax, givenBins, causal_scores, linked_scores, neut_scores, xlabelstring, ylabelstring, xlims, ylims, transparency, thinToSize = True, show = False, save = True, savename = "test.png", titlestring = ''):
	#if desired; downsample linked/neut variants according to how many selected variants we have data for
	if thinToSize == True: 
		limiting = min(len(causal_scores), len(linked_scores), len(neut_scores))
		print "Thinning data to " + str(limiting) + " SNPs..."
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

	causal_scores = np.clip(causal_scores, xlims[0], xlims[1])
	linked_scores = np.clip(linked_scores, xlims[0], xlims[1])
	neut_scores = np.clip(neut_scores, xlims[0], xlims[1])
	
	#plot 
	n_causal, bins_causal, patches_causal = ax.hist(causal_scores, bins=givenBins, weights = weights_causal, color = "red", alpha=transparency, label="causal; n = " + str(len(causal_scores)))
	n_linked, bins_linked, patches_linked = ax.hist(linked_scores,  bins=givenBins, weights=weights_linked, color = "green",alpha=transparency,  label="linked; n = " + str(len(linked_scores)))
	n_neut, bins_neut, patches_neut = ax.hist(neut_scores, bins=givenBins, weights=weights_neut, color = "blue",alpha=transparency,  label="neutral; n = " + str(len(neut_scores)))

	#print bins_causal
	#print givenBins

	writefilename = 'scores/likes/' + titlestring + "_causal.txt"
	writefile = open(writefilename, 'w') 
	#writefile.write(str(givenBins) +"\n")
	for index in range(len(n_causal)):
		towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(n_causal[index])+ "\n"
		writefile.write(towritestring)
	writefile.close()
	print "wrote to file: " + writefilename


	writefilename = 'scores/likes/' + titlestring + "_linked.txt"
	writefile = open(writefilename, 'w') 
	#writefile.write(str(givenBins) +"\n")
	for index in range(len(n_causal)):
		towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(n_linked[index])+ "\n"
		writefile.write(towritestring)
	writefile.close()
	print "wrote to file: " + writefilename

	writefilename = 'scores/likes/' + titlestring + "_neut.txt"
	writefile = open(writefilename, 'w') 
	#writefile.write(str(givenBins) +"\n")
	for index in range(len(n_causal)):
		towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(n_neut[index])+ "\n"
		writefile.write(towritestring)
	writefile.close()
	print "wrote to file: " + writefilename


	#annotate and display/save
	#plt.title(titlestring)
	ax.set_xlabel(xlabelstring, size=8)
	ax.set_ylabel(ylabelstring, size=8)
	ax.set_ylim(ylims)
	#plt.legend()
	if show:
		plt.show()
	if save:
		plt.savefig(savename)
		print "saved figure to: " + savename
	return