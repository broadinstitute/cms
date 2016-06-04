## {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
## lifted from /idi/sabeti-scratch/jvitti/synth/ e.g. func_files.py 01.14.16 ; func_aws_calcs.py 04.04.16
## last updated 06.04.16 	vitti@broadinstitute.org

from func_clean import *
import matplotlib.pyplot as plt
from math import sqrt, fabs
import subprocess
import pylab
import os

######################
## MANIPULATE FILES ##
######################
def givecommand(commandstr):
	commands = commandstr.split(' ')
	subprocess.check_output(commands)
def checkFileExists(filename):
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		return filename
	else:
		if ".gz" in filename:
			filename = filename.strip('.gz')
			if os.path.isfile(filename) and os.path.getsize(filename) > 0:
				return filename
		else:
			filename = filename + ".gz"
			if os.path.isfile(filename) and os.path.getsize(filename) > 0:
				return filename
		#print "missing: " + filename
		return "MISSING"\
def getFileNames(pop1, pop2, chrom, vers="v5a", mask = 'strictMask', runtime = '062215_2pm_nomultialleles', suffix = "_merged_121915"):
	"""returns CMS filenames for passing to combine_cms_scores.c to generate pairwise comparisons.
	tracks pop order to facilitate comparisons with Booleans rev_xp and rev_fstdeldaf
	"""

	ihs1filename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/" + vers + "_autosomes_" + mask + "_" + pop1 + "_" + runtime + suffix + ".ihs.out.100bins.norm.chr" + str(chrom) + ".txt"
	ihs2filename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/" + vers + "_autosomes_" + mask + "_" + pop2 + "_" + runtime + suffix + ".ihs.out.100bins.norm.chr" + str(chrom) + ".txt"

	delihh1filename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/merged/" + vers + "_autosomes_" + mask + "_" + pop1 + "_" + runtime + suffix + ".delihh.out.100bins.norm.chr" + str(chrom) + ".txt"
	delihh2filename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/merged/" + vers + "_autosomes_" + mask + "_" + pop2 + "_" + runtime + suffix + ".delihh.out.100bins.norm.chr" + str(chrom) + ".txt"

	xpfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask +"_" + runtime + suffix + "_" + pop1 +"_" + pop2 +".xpehh.norm.chr" + str(chrom) + ".txt"
	xpehh_loc = checkFileExists(xpfilename)
	rev_xp = False
	if xpehh_loc == "MISSING":
		xpfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask +"_" + runtime + suffix + "_" + pop2 +"_" + pop1 +".xpehh.norm.chr" + str(chrom) + ".txt"
		rev_xp = True

	rev_fstdeldaf = False
	fstdeldaffilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop1 + "_" + pop2 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
	fst_loc = checkFileExists(fstdeldaffilename)
	if fst_loc == "MISSING":
		fstdeldaffilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop2 + "_" + pop1 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
		rev_fstdeldaf = True

	for filename in [ihs1filename, ihs2filename, delihh1filename, delihh2filename, xpfilename, fstdeldaffilename]:
		if checkFileExists(filename) == "MISSING":
			print "HE--? DOKO NI IT-TA NO--?: " + filename +"\n"
	return ihs1filename, ihs2filename, delihh1filename, delihh2filename, xpfilename, fstdeldaffilename, rev_xp, rev_fstdeldaf
def plotStatsHist(stats, title, xlabel, ylabel, savename, nbins = 60,show=True, save=False):
	#weights_stats = np.ones_like(stats)/len(stats)
	plt.figure()
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	n_stats, bin_stats, patches = pylab.hist(stats,  bins=nbins, histtype = 'step') #weights = weights_stats,
	if show == True:
		plt.show()
	if save == True:
		plt.savefig(savename)
	plt.close()
	return n_stats, bin_stats
def gatherStats(infiles, takeIndex = -2, stripHeader = False, verbose = False):
	toTake = []
	for infile in infiles:
		openfile = open(infile, 'r')
		if stripHeader == True:
			openfile.readline()
		for line in openfile:
			entries = line.split() 
			toTake.append(entries[takeIndex])
		openfile.close()
	if verbose:
		print "\tloaded " + str(len(toTake)) + " values from " + str(len(infiles)) + " files, using index = " + str(takeIndex)
		print infiles[0]
	return toTake

######################
## CALCULATE SCORES ##
######################

def launch_calc_ihs(scenario, nRep, numThreads = 35, ihsloc = '/SIMS2/cosi_ihs/', tpedloc ="/SIMS2/cosi_tpeds/", pops = [1, 2, 3, 4]):
	mkdircommand = "mkdir " + ihsloc + scenario
	givecommand(mkdircommand)
	taskargstrings = []
	for pop in pops:
		for iRep in range(1,nRep+1):
			fileprefix = tpedloc + scenario + "/rep" + str(iRep) 
			inputTped = fileprefix + "_pop" + str(pop) + ".tped"
			outputFile = ihsloc + scenario + "/rep" +  str(iRep) + "_pop" + str(pop) + "_ihs"
			argstring = "selscan_ihs " + inputTped + " " + outputFile + " --threads " + str(numThreads)
			taskargstrings.append(argstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_ihs_" + scenario)		
def launch_calc_ihs_batch(scenario, inputs, outputs, numThreads = 35):
	taskargstrings = []
	for i in range(len(inputs)):
		inputTped = inputs[i]
		outputFile = outputs[i]
		argstring = "selscan_ihs " + inputTped + " " + outputFile + " --threads " + str(numThreads)
		taskargstrings.append(argstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_ihs_" + scenario)	
def write_delIHH_file(readfilename, writefilename):
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
def launch_calcdelihh(scenario, ihsloc, delihhloc):
	"""n.b. must have normalized ihs values first to get concatenated file"""
	for pop in pops:
		inputFilename = ihsloc + "/pop" + str(pop) + ".ihs.concat"		
		assert os.isfile(inputFilename)		
		outputFile = delihhloc + "/pop" + str(pop) + ".delihh.concat"
		write_delIHH_file(inputFilename, outputFile)
		print "wrote file to: " + outputFile
	return
def launch_calc_xp(scenario, nRep,numThreads = 35,  Xploc = '/SIMS2/cosi_xpehh/',tpedloc ="/SIMS2/cosi_tpeds/", popPairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]):
	mkdircommand = "mkdir " + Xploc + scenario
	givecommand(mkdircommand)
	taskargstrings = []
	for popPair in popPairs:
		pop1, pop2 = popPair[0], popPair[1]
		for iRep in range(1, nRep+1):
			fileprefix = tpedloc + scenario + "/rep" + str(iRep) 
			inputTped1 = fileprefix + "_pop" + str(pop1) + ".tped"
			inputTped2 = fileprefix + "_pop" + str(pop2) + ".tped"
			outputFile = Xploc + scenario + "/rep" +  str(iRep) + "_pop" + str(pop1) + "_pop" + str(pop2) + "_xpehh"
			argstring = "selscan_xpehh " + inputTped1 + " " + outputFile +  " " + inputTped2 + " --threads " + str(numThreads)
			taskargstrings.append(argstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_xp_" + scenario)
def launch_calc_xp_batch(scenario, inputs_1, inputs_2, outputs, numThreads = 35):
	taskargstrings = []
	for i in range(len(inputs_1)):
		inputTped1 = inputs_1[i]
		inputTped2 = inputs_2[i]
		outputFile = outputs[i]
		argstring = "selscan_xpehh " + inputTped1 + " " + outputFile +  " " + inputTped2 + " --threads " + str(numThreads)
		taskargstrings.append(argstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_xp_" + scenario)	
def launch_fst_deldaf(scenario, nRep, tpedloc, fst_loc, popPairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]):
	taskargstrings = []
	for popPair in popPairs:
		pop1, pop2 = popPair[0], popPair[1]
		for iRep in range(1, nRep+1):
			fileprefix = tpedloc + "rep" + str(iRep) 
			inputTped1 = fileprefix + "_pop" + str(pop1) + ".tped"
			inputTped2 = fileprefix + "_pop" + str(pop2) + ".tped"
			outputFile = fst_loc  + "rep" +  str(iRep) + "_pop" + str(pop1) + "_pop" + str(pop2) + "_fst_deldaf.txt"
			argstring = inputTped1 + " " + inputTped2 + " " + outputFile
			taskargstrings.append(argstring)
	commandstring = "./calc_Fst_delDAF_cosi_tped"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "calc_fst_deldaf_" + scenario)
	return


######################
## NORMALIZE SCORES ##
######################

def readNeutNormFile(neutnormfilename, scoretype ='xp'):
def normalizeiHS(scenario, nRep, ihsloc, pops = [1, 2, 3, 4]):
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
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "norm_ihs_" + scenario)		
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

	#print "now normalizing with the following neutral parameters: " 
	#print str(means_bin)
	#print str(vars_bin)

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
	return
def normalize_delIhh(scenario, delihhloc, pops = [1, 2, 3, 4]):
	taskargstrings = []
	for pop in pops:
		delihhfilename = delihhloc  + "pop" + str(pop) + ".delihh.concat"
		assert os.path.isfile(delihhfilename)
		taskargstring = delihhfilename 
		taskargstrings.append(taskargstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py selscan_norm_ihs"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "norm_delihh_" + scenario)		
	return
def normalize_delihh_manually(unnormedfilename, means_bin, vars_bin, bin_bounds):
	normfilename = unnormedfilename + ".norm"
	openfile = open(unnormedfilename, 'r')
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
	return
def normalizeXpehh(scenario, nRep, Xploc, popPairs = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]):
	taskargstrings = []
	for popPair in popPairs:
		pop1, pop2 = popPair[0], popPair[1]
		nsnps = 0
		writefilename = Xploc  + "pop" + str(pop1) + "_pop" + str(pop2) + ".xpehh.concat"
		writefile = open(writefilename, 'w')
		for iRep in range(nRep):
			xpehhfile = Xploc  + "rep" + str(iRep) + "_pop" + str(pop1) + "_pop" + str(pop2) + ".xpehh.out"
			#print xpehhfile
			if os.path.isfile(xpehhfile):
				readfile = open(xpehhfile, 'r')
				header=readfile.readline() #strip header
				if iRep == 0:
					writefile.write(header)
				for line in readfile:
					writefile.write(line)
					nsnps +=1
				readfile.close()
		writefile.close()
		print "for pop pair " + str(popPair) + " we logged " + str(nsnps) + " scores to " + writefilename
		taskargstring = "selscan_norm_xpehh " + writefilename 
		taskargstrings.append(taskargstring)
	commandstring = "python2.7 /SIMS/cms_venv/cms/cms/scans.py"
	writeJobsToTaskArray_fromArglist(commandstring, taskargstrings, "norm_xp_" + scenario)	
	return	
def	normalize_Xp_manually(xpdir, nRep, mean, var, pop1, pop2):
	"""allows us to pass the mean and var from neutral simulations to normalize selection sims accordingly"""

	nsnps = 0
	writefilename = xpdir  + "pop" + str(pop1) + "_pop" + str(pop2) + ".xpehh.concat"
	writefile = open(writefilename, 'w')
	for iRep in range(nRep):
		xpehhfile = xpdir + "rep" + str(iRep) + "_pop" + str(pop1) + "_pop" + str(pop2) + ".xpehh.out"
		#print xpehhfile
		if os.path.isfile(xpehhfile):
			readfile = open(xpehhfile, 'r')
			header=readfile.readline() #strip header
			if iRep == 0:
				writefile.write(header)
			for line in readfile:
				writefile.write(line)
				nsnps +=1
			readfile.close()
	writefile.close()
	print "for pops " + str(pop1) + "-" +str(pop2) + " we logged " + str(nsnps) + " scores to " + writefilename

	print "now normalizing with the following neutral parameters: " 
	print str(mean)
	print str(var)

	normfilename = writefilename + ".norm"
	openfile = open(writefilename, 'r')
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