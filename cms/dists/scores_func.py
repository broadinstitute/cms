## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 03.24.17 vitti@broadinstitute.org

from scipy.stats.kde import gaussian_kde
from math import fabs, sqrt
from random import randint
import numpy as np
import subprocess
import sys
import os

#####################
## CALCS FROM SIMS ## are these necessary? maybe replace with 
##################### all-in-one function
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
		unstand_delIHH = fabs(float(ihh_1) - float(ihh_0)) #WAIT WHY FABS?????????
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
	print(fullcommand)
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
def get_sim_compscore_files(pop, repNum, basedir):
	''' pulls all replicates with complete normalized component score files'''
	pops = [1, 2, 3, 4]
	altpops = pops[:]
	altpops.remove(int(pop))
	check_files = []
	ihs_normedfile = basedir + "ihs/rep" + str(repNum) + "_" + str(pop) + ".ihs.out.norm"
	delihh_normedfile =  basedir + "delihh/rep" + str(repNum) + "_" + str(pop) + ".txt.norm"
	nsl_normedfileprefix = basedir + "nsl/rep" + str(repNum) + "_" + str(pop) + ".nsl.out.norm"
	check_files.extend([ihs_normedfile, delihh_normedfile, nsl_normedfileprefix])
	for altpop in altpops:
		xpehh_normedfile = basedir + "xpehh/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out.norm"
		fstdeldaf_outfilename = basedir + "fst_deldaf/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop)
		check_files.extend([xpehh_normedfile, fstdeldaf_outfilename])
	return check_files
def get_scores_from_files(all_completed_neut, all_completed_sel, scoreindex, sel_bin_index):
	''' '''
	neut_files1 = [all_completed_neut[0][irep][scoreindex] for irep in range(len(all_completed_neut[0]))]
	neut_files2 = [all_completed_neut[1][irep][scoreindex] for irep in range(len(all_completed_neut[1]))]
	neut_files3 = [all_completed_neut[2][irep][scoreindex] for irep in range(len(all_completed_neut[2]))]
	neut_files4 = [all_completed_neut[3][irep][scoreindex] for irep in range(len(all_completed_neut[3]))]
	neut_values1 = load_from_files(neut_files1)
	neut_values2 = load_from_files(neut_files2)
	neut_values3 = load_from_files(neut_files3)
	neut_values4 = load_from_files(neut_files4)
	print("loaded " + str(len(neut_values1)) + " neutral values for pop 1...")
	print("loaded " + str(len(neut_values2)) + " neutral values for pop 2...")	
	print("loaded " + str(len(neut_values3)) + " neutral values for pop 3...")			
	print("loaded " + str(len(neut_values4)) + " neutral values for pop 4...")

	sel_files1 = [all_completed_sel[0][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
	sel_files2 = [all_completed_sel[1][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
	sel_files3 = [all_completed_sel[2][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
	sel_files4 = [all_completed_sel[3][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[3][sel_bin_index]))]		
	causal_values1, linked_values1 = load_from_files_discriminate_causal(sel_files1)
	causal_values2, linked_values2 = load_from_files_discriminate_causal(sel_files2)
	causal_values3, linked_values3 = load_from_files_discriminate_causal(sel_files3)
	causal_values4, linked_values4 = load_from_files_discriminate_causal(sel_files4)
	print("loaded " + str(len(causal_values1)) + " causal SNP and " + str(len(linked_values1)) + " linked SNP iHS values for pop 1...")
	print("loaded " + str(len(causal_values2)) + " causal SNP and " + str(len(linked_values2)) + " linked SNP iHS values for pop 2...")
	print("loaded " + str(len(causal_values3)) + " causal SNP and " + str(len(linked_values3)) + " linked SNP iHS values for pop 3...")
	print("loaded " + str(len(causal_values4)) + " causal SNP and " + str(len(linked_values4)) + " linked SNP iHS values for pop 4...")

	all_score_values = [[neut_values1, causal_values1, linked_values1],
						[neut_values2, causal_values2, linked_values2],
						[neut_values3, causal_values3, linked_values3],
						[neut_values4, causal_values4, linked_values4],]
	return all_score_values
def get_compscores_from_files(all_completed_neut, all_completed_sel, scorestring, sel_bin_index):
	''' '''
	if scorestring in ['fst', 'deldaf']:
		physIndex = 0
		neut_files1a = [all_completed_neut[0][irep][4] for irep in range(len(all_completed_neut[0]))]
		neut_files1b = [all_completed_neut[0][irep][6] for irep in range(len(all_completed_neut[0]))]
		neut_files1c = [all_completed_neut[0][irep][8] for irep in range(len(all_completed_neut[0]))]
		neut_files2a = [all_completed_neut[1][irep][4] for irep in range(len(all_completed_neut[1]))]
		neut_files2b = [all_completed_neut[1][irep][6] for irep in range(len(all_completed_neut[1]))]
		neut_files2c = [all_completed_neut[1][irep][8] for irep in range(len(all_completed_neut[1]))]
		neut_files3a = [all_completed_neut[2][irep][4] for irep in range(len(all_completed_neut[2]))]
		neut_files3b = [all_completed_neut[2][irep][6] for irep in range(len(all_completed_neut[2]))]
		neut_files3c = [all_completed_neut[2][irep][8] for irep in range(len(all_completed_neut[2]))]
		neut_files4a = [all_completed_neut[3][irep][4] for irep in range(len(all_completed_neut[3]))]
		neut_files4b = [all_completed_neut[3][irep][6] for irep in range(len(all_completed_neut[3]))]
		neut_files4c = [all_completed_neut[3][irep][8] for irep in range(len(all_completed_neut[3]))]
		sel_files1a = [all_completed_sel[0][sel_bin_index][irep][4] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files1b = [all_completed_sel[0][sel_bin_index][irep][6] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files1c = [all_completed_sel[0][sel_bin_index][irep][8] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files2a = [all_completed_sel[1][sel_bin_index][irep][4] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files2b = [all_completed_sel[1][sel_bin_index][irep][6] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files2c = [all_completed_sel[1][sel_bin_index][irep][8] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files3a = [all_completed_sel[2][sel_bin_index][irep][4] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files3b = [all_completed_sel[2][sel_bin_index][irep][6] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files3c = [all_completed_sel[2][sel_bin_index][irep][8] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files4a = [all_completed_sel[3][sel_bin_index][irep][4] for irep in range(len(all_completed_sel[3][sel_bin_index]))]
		sel_files4b = [all_completed_sel[3][sel_bin_index][irep][6] for irep in range(len(all_completed_sel[3][sel_bin_index]))]
		sel_files4c = [all_completed_sel[3][sel_bin_index][irep][8] for irep in range(len(all_completed_sel[3][sel_bin_index]))] 
	else:
		assert(scorestring == "xpehh") #this could be made efficient
		physIndex = 1
		neut_files1a = [all_completed_neut[0][irep][3] for irep in range(len(all_completed_neut[0]))]
		neut_files1b = [all_completed_neut[0][irep][5] for irep in range(len(all_completed_neut[0]))]
		neut_files1c = [all_completed_neut[0][irep][7] for irep in range(len(all_completed_neut[0]))]
		neut_files2a = [all_completed_neut[1][irep][3] for irep in range(len(all_completed_neut[1]))]
		neut_files2b = [all_completed_neut[1][irep][5] for irep in range(len(all_completed_neut[1]))]
		neut_files2c = [all_completed_neut[1][irep][7] for irep in range(len(all_completed_neut[1]))]
		neut_files3a = [all_completed_neut[2][irep][3] for irep in range(len(all_completed_neut[2]))]
		neut_files3b = [all_completed_neut[2][irep][5] for irep in range(len(all_completed_neut[2]))]
		neut_files3c = [all_completed_neut[2][irep][7] for irep in range(len(all_completed_neut[2]))]
		neut_files4a = [all_completed_neut[3][irep][3] for irep in range(len(all_completed_neut[3]))]
		neut_files4b = [all_completed_neut[3][irep][5] for irep in range(len(all_completed_neut[3]))]
		neut_files4c = [all_completed_neut[3][irep][7] for irep in range(len(all_completed_neut[3]))]
		sel_files1a = [all_completed_sel[0][sel_bin_index][irep][3] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files1b = [all_completed_sel[0][sel_bin_index][irep][5] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files1c = [all_completed_sel[0][sel_bin_index][irep][7] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
		sel_files2a = [all_completed_sel[1][sel_bin_index][irep][3] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files2b = [all_completed_sel[1][sel_bin_index][irep][5] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files2c = [all_completed_sel[1][sel_bin_index][irep][7] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
		sel_files3a = [all_completed_sel[2][sel_bin_index][irep][3] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files3b = [all_completed_sel[2][sel_bin_index][irep][5] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files3c = [all_completed_sel[2][sel_bin_index][irep][7] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
		sel_files4a = [all_completed_sel[3][sel_bin_index][irep][3] for irep in range(len(all_completed_sel[3][sel_bin_index]))]
		sel_files4b = [all_completed_sel[3][sel_bin_index][irep][5] for irep in range(len(all_completed_sel[3][sel_bin_index]))]
		sel_files4c = [all_completed_sel[3][sel_bin_index][irep][7] for irep in range(len(all_completed_sel[3][sel_bin_index]))]

	neut_values1a = load_from_files(neut_files1a, stripHeader=True)
	neut_values1b = load_from_files(neut_files1b, stripHeader=True)
	neut_values1c = load_from_files(neut_files1c, stripHeader=True)
	neut_values2a = load_from_files(neut_files2a, stripHeader=True)
	neut_values2b = load_from_files(neut_files2b, stripHeader=True)
	neut_values2c = load_from_files(neut_files2c, stripHeader=True)
	neut_values3a = load_from_files(neut_files3a, stripHeader=True)
	neut_values3b = load_from_files(neut_files3b, stripHeader=True)
	neut_values3c = load_from_files(neut_files3c, stripHeader=True)
	neut_values4a = load_from_files(neut_files4a, stripHeader=True)
	neut_values4b = load_from_files(neut_files4b, stripHeader=True)
	neut_values4c = load_from_files(neut_files4c, stripHeader=True)
	causal_values1a, linked_values1a = load_from_files_discriminate_causal(sel_files1a, stripHeader=True, physIndex=physIndex)
	causal_values1b, linked_values1b = load_from_files_discriminate_causal(sel_files1b, stripHeader=True, physIndex=physIndex)
	causal_values1c, linked_values1c = load_from_files_discriminate_causal(sel_files1c, stripHeader=True, physIndex=physIndex)
	causal_values2a, linked_values2a = load_from_files_discriminate_causal(sel_files2a, stripHeader=True, physIndex=physIndex)
	causal_values2b, linked_values2b = load_from_files_discriminate_causal(sel_files2b, stripHeader=True, physIndex=physIndex)
	causal_values2c, linked_values2c = load_from_files_discriminate_causal(sel_files2c, stripHeader=True, physIndex=physIndex)
	causal_values3a, linked_values3a = load_from_files_discriminate_causal(sel_files3a, stripHeader=True, physIndex=physIndex)
	causal_values3b, linked_values3b = load_from_files_discriminate_causal(sel_files3b, stripHeader=True, physIndex=physIndex)
	causal_values3c, linked_values3c = load_from_files_discriminate_causal(sel_files3c, stripHeader=True, physIndex=physIndex)
	causal_values4a, linked_values4a = load_from_files_discriminate_causal(sel_files4a, stripHeader=True, physIndex=physIndex)
	causal_values4b, linked_values4b = load_from_files_discriminate_causal(sel_files4b, stripHeader=True, physIndex=physIndex)
	causal_values4c, linked_values4c = load_from_files_discriminate_causal(sel_files4c, stripHeader=True, physIndex=physIndex)
	print("loaded " + str(sum([len(neut_values1a), len(neut_values1b), len(neut_values1c)])) + " neutral values for pop 1 vs three outgroups...")
	print("loaded " + str(sum([len(neut_values2a), len(neut_values2b), len(neut_values2c)])) + " neutral values for pop 2 vs three outgroups...")
	print("loaded " + str(sum([len(neut_values3a), len(neut_values3b), len(neut_values3c)])) + " neutral values for pop 3 vs three outgroups...")
	print("loaded " + str(sum([len(neut_values4a), len(neut_values4b), len(neut_values4c)])) + " neutral values for pop 4 vs three outgroups...")
	
	print("loaded " + str(sum([len(linked_values1a), len(linked_values1b), len(linked_values1c)])) + " linked values for pop 1 vs three outgroups...")
	print("loaded " + str(sum([len(linked_values2a), len(linked_values2b), len(linked_values2c)])) + " linked values for pop 2 vs three outgroups...")
	print("loaded " + str(sum([len(linked_values3a), len(linked_values3b), len(linked_values3c)])) + " linked values for pop 3 vs three outgroups...")
	print("loaded " + str(sum([len(linked_values4a), len(linked_values4b), len(linked_values4c)])) + " linked values for pop 4 vs three outgroups...")

	print("loaded " + str(sum([len(causal_values1a), len(causal_values1b), len(causal_values1c)])) + " causal values for pop 1 vs three outgroups...")
	print("loaded " + str(sum([len(causal_values2a), len(causal_values2b), len(causal_values2c)])) + " causal values for pop 2 vs three outgroups...")
	print("loaded " + str(sum([len(causal_values3a), len(causal_values3b), len(causal_values3c)])) + " causal values for pop 3 vs three outgroups...")
	print("loaded " + str(sum([len(causal_values4a), len(causal_values4b), len(causal_values4c)])) + " causal values for pop 4 vs three outgroups...")

	all_score_values = [[[neut_values1a, causal_values1a, linked_values1a],[neut_values1b, causal_values1b, linked_values1b],
						[neut_values1c, causal_values1c, linked_values1c]],
						[[neut_values2a, causal_values2a, linked_values2a],[neut_values2b, causal_values2b, linked_values2b],
						[neut_values2c, causal_values2c, linked_values2c]],
						[[neut_values3a, causal_values3a, linked_values3a],[neut_values3b, causal_values3b, linked_values3b],
						[neut_values3c, causal_values3c, linked_values3c]],
						[[neut_values4a, causal_values4a, linked_values4a],[neut_values4b, causal_values4b, linked_values4b],
						[neut_values4c, causal_values4c, linked_values4c]],
						] #[pop][outgroup][type]
	return all_score_values
def load_from_files(files, takeIndex = -1, stripHeader = False):
	#print('loading from ' + str(len(files)) + " files...")
	values = []
	for file in files:
		openfile = open(file, 'r')
		if stripHeader:
			openfile.readline()
		for line in openfile:
			entries = line.split()
			value = float(entries[takeIndex])
			values.append(value)
		openfile.close()
	return values
def load_from_files_discriminate_causal(files, takeIndex = -1, stripHeader = False, causalLoc = 750000, physIndex = 1):
	#print('loading from ' + str(len(files)) + " files...")
	causal_values, linked_values = [], []
	for file in files:
		openfile = open(file, 'r')
		if stripHeader:
			openfile.readline()
		for line in openfile:
			entries = line.split()
			value = float(entries[takeIndex])
			if int(entries[physIndex]) == causalLoc:
				causal_values.append(value)
			else:
				linked_values.append(value)
		openfile.close()
	return causal_values, linked_values	
def plot_pdf_comparison_from_scores(ax, neutvals, causalvals, linkedvals, minVal, maxVal, nProbBins, annotate = False):
	neutvals = np.array(neutvals)
	neutvals = neutvals[~np.isnan(neutvals)]
	causalvals = np.array(causalvals)
	causalvals = causalvals[~np.isnan(causalvals)]
	linkedvals = np.array(linkedvals)
	linkedvals = linkedvals[~np.isnan(linkedvals)]
	if annotate:
		mean = np.mean(neutvalues)
		var = np.var(neutvalues)
		sd = np.sqrt(var)
		annotationstring = "NEUT max: " + str(max(values)) + "\nmin: " + str(min(values)) + "\nmean: " + str(np.mean(values)) + "\nvar: " + str(np.var(values))
		xlims = ax.get_xlim()
		ylims = ax.get_ylim()
		fullxrange = xlims[1] - xlims[0]
		fullyrange = ylims[1] - ylims[0]
		anno_x = xlims[1] - (fullxrange/3)
		anno_y = ylims[1] - (fullyrange/2)
		ax.annotate(annotationstring, xy=(anno_x, anno_y), fontsize=6)
	kde_neut = gaussian_kde( neutvals )
	kde_causal = gaussian_kde( causalvals )
	kde_linked = gaussian_kde( linkedvals )
	dist_space = np.linspace(minVal, maxVal, nProbBins)
	ax.plot( dist_space, kde_neut(dist_space) , color='blue')
	ax.plot( dist_space, kde_causal(dist_space) , color='red')
	ax.plot( dist_space, kde_linked(dist_space) , color='green')
	return
def get_plot_pdf_params(score):
	#for now a placeholder -- merge with old hist func? 
	return 3, 3, 60, False

"""
def get_indices(score, dem_scenario):
	#CRUCIAL FUNC; tells loadVals which columns to grab and return; changes with filetype. also returns anc_freq_index
	if score == "ihs":
		physpos_index, freq_anc_index = 1, 2
		#if "sel" in dem_scenario:
		if True: #PULL FROM NORMED FILE / flexible to old format
			ihs_unnormed_index, ihs_normed_index, expectedlen = -2, -1, 11 #these files have 11 columns
		#else: #neutral
		#	ihs_unnormed_index, ihs_normed_index, expectedlen = 5, 6, 10 #these files have 8 columns; last one is binary variable
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
	elif score == "fst": #revert#PULL FROM COMP (.CMS FILE)
		#expectedlen = 8
		#physpos_index, scoreindex, freq_anc_index = 0, 5, float('nan')
		expectedlen = 3#2
		physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "deldaf":#revert#PULL FROM COMP (.CMS FILE)
		#expectedlen = 8
		#physpos_index, scoreindex, freq_anc_index = 0, 6, float('nan')	
		expectedlen = 3#2
		physpos_index, scoreindex, freq_anc_index = 0, 1, 2
		indices = [physpos_index, scoreindex, freq_anc_index]
	elif score == "xp": #revert#PULL FROM COMP (.CMS FILE)
		#expectedlen = 8
		#physpos_index, scoreindex, freq_anc_index = 0, 4, float('nan')	
		physpos_index, freq_anc_index = 1, 2
		xp_unnormed_index, xp_normed_index, = 7, 8
		expectedlen=9
		scoreindex = xp_normed_index
		#if "sel" in dem_scenario:
		#	expectedlen = 9
		#else:
		#	expectedlen = 10
		indices = [physpos_index, scoreindex, freq_anc_index]

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
		if os.path.getsize(filename) > 0:
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
					if np.isnan(takeindices[-1]):
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
	#weights_causal = np.ones_like(causal_scores)/len(causal_scores)
	#weights_linked = np.ones_like(linked_scores)/len(linked_scores)
	#weights_neut = np.ones_like(neut_scores)/len(neut_scores)

	causal_scores = np.clip(causal_scores, xlims[0], xlims[1]) #np.clip
	linked_scores = np.clip(linked_scores, xlims[0], xlims[1])
	neut_scores = np.clip(neut_scores, xlims[0], xlims[1])

	n_causal, bins_causal = np.histogram(causal_scores, range=xlims, bins=givenBins)#, weights = weights_causal)
	n_linked, bins_linked = np.histogram(linked_scores, range=xlims,  bins=givenBins)#, weights = weights_linked)
	n_neut, bins_neut = np.histogram(neut_scores,range=xlims, bins=givenBins)#, weights = weights_neut)

	#debug_array = [n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut]
	#print(debug_array)

	#totalNsnps = len(causal_scores) + len(linked_scores) + len(neut_scores)

	#for ibin in range(len(n_causal)):
	#	if n_causal[ibin] <= 1:
	#		n_causal[ibin] = 1e-10
	#	if n_linked[ibin] <= 1:
	#		n_linked[ibin] = 1e-10
	#	if n_neut[ibin] <=1:
	#		n_neut[ibin] = 1e-10

	return n_causal, n_linked, n_neut, bins_causal, bins_linked, bins_neut
def write_hists_to_files(writePrefix, givenBins, n_causal, n_linked, n_neut):
	assert len(givenBins) == (len(n_causal) + 1)
	for status in ['causal', 'linked', 'neut']:
		writefilename = writePrefix + "_" + status + ".txt"
		writefile = open(writefilename, 'w')
		n_scores = eval('n_' + status)
		for index in range(len(n_scores)):
			num_in_bin = n_scores[index]
			if (num_in_bin <= 1):
				writeprob = 1e-10
			else:
				if (index < (len(n_scores) - 1) and index > 0): #check neighbors
					if (n_scores[index+1] == n_scores[index-1]) and (n_scores[index+1] <=1):
						writeprob = 1e-10
					else:
						writeprob = float(n_scores[index])/(sum(n_scores)) #
				else:
					writeprob = float(n_scores[index])/(sum(n_scores)) #
			towritestring =  str(givenBins[index]) + "\t" + str(givenBins[index+1]) + "\t" + str(writeprob)+ "\n"
			writefile.write(towritestring)
		writefile.close()
	return
"""