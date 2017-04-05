## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 04.05.2017 vitti@broadinstitute.org

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
def get_pbs_from_dafs(seldaf, outdaf_1, outdaf_2):
	in_fst_1 = quick_fst_estimator(seldaf, outdaf_1) #could build in a check here to ensure
	in_fst_2 = quick_fst_estimator(seldaf, outdaf_2) #that these align with previous calculations
	out_fst = quick_fst_estimator(outdaf_1, outdaf_2)
	in_t_1 = -1 * np.log(1.-in_fst_1)
	in_t_2 = -1 * np.log(1.-in_fst_2)
	out_t = -1 * np.log(1.-out_fst)
	pbs = (in_t_1 + in_t_2 - out_t) / 2.
	return pbs
def quick_fst_estimator(daf1, daf2, nperPop=172):
	''' helper method to get_pbs_from_dafs '''
	daf_mean = (daf1 + daf2) / 2.
	msp = nperPop * (daf1 - daf_mean) * (daf1 - daf_mean) + nperPop * (daf2 - daf_mean) * (daf2 - daf_mean)
	msg = (nperPop * daf1 * (1. - daf1) + nperPop * daf2 * (1. - daf2)) / (nperPop - 1 + nperPop - 1)
	num = msp - msg
	denom = msp + ((nperPop) - 1) * msg
	if denom != 0:
		fst_hat = num / denom
	else:
		fst_hat = 0
	return fst_hat

#############################
## MANIPULATE SCORE FILES ###
#############################
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
def get_scores_from_files(all_completed_neut, all_completed_sel, scoreindex, sel_bin_index, startbound, endbound, foldDists = False):
	''' '''
	neut_files1 = [all_completed_neut[0][irep][scoreindex] for irep in range(len(all_completed_neut[0]))]
	neut_files2 = [all_completed_neut[1][irep][scoreindex] for irep in range(len(all_completed_neut[1]))]
	neut_files3 = [all_completed_neut[2][irep][scoreindex] for irep in range(len(all_completed_neut[2]))]
	neut_files4 = [all_completed_neut[3][irep][scoreindex] for irep in range(len(all_completed_neut[3]))]
	neut_values1 = load_from_files(neut_files1,  startbound, endbound, absVal = foldDists)
	neut_values2 = load_from_files(neut_files2,  startbound, endbound, absVal = foldDists)
	neut_values3 = load_from_files(neut_files3,  startbound, endbound, absVal = foldDists)
	neut_values4 = load_from_files(neut_files4,  startbound, endbound, absVal = foldDists)
	print("loaded " + str(len(neut_values1)) + " neutral values for pop 1...")
	print("loaded " + str(len(neut_values2)) + " neutral values for pop 2...")	
	print("loaded " + str(len(neut_values3)) + " neutral values for pop 3...")			
	print("loaded " + str(len(neut_values4)) + " neutral values for pop 4...")

	sel_files1 = [all_completed_sel[0][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[0][sel_bin_index]))]
	sel_files2 = [all_completed_sel[1][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[1][sel_bin_index]))]
	sel_files3 = [all_completed_sel[2][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[2][sel_bin_index]))]
	sel_files4 = [all_completed_sel[3][sel_bin_index][irep][scoreindex] for irep in range(len(all_completed_sel[3][sel_bin_index]))]		
	causal_values1, linked_values1 = load_from_files_discriminate_causal(sel_files1,  startbound, endbound, absVal = foldDists)
	causal_values2, linked_values2 = load_from_files_discriminate_causal(sel_files2,  startbound, endbound, absVal = foldDists)
	causal_values3, linked_values3 = load_from_files_discriminate_causal(sel_files3,  startbound, endbound, absVal = foldDists)
	causal_values4, linked_values4 = load_from_files_discriminate_causal(sel_files4,  startbound, endbound, absVal = foldDists)
	print("loaded " + str(len(causal_values1)) + " causal SNP and " + str(len(linked_values1)) + " linked SNP iHS values for pop 1...")
	print("loaded " + str(len(causal_values2)) + " causal SNP and " + str(len(linked_values2)) + " linked SNP iHS values for pop 2...")
	print("loaded " + str(len(causal_values3)) + " causal SNP and " + str(len(linked_values3)) + " linked SNP iHS values for pop 3...")
	print("loaded " + str(len(causal_values4)) + " causal SNP and " + str(len(linked_values4)) + " linked SNP iHS values for pop 4...")

	all_score_values = [[neut_values1, causal_values1, linked_values1],
						[neut_values2, causal_values2, linked_values2],
						[neut_values3, causal_values3, linked_values3],
						[neut_values4, causal_values4, linked_values4],]
	return all_score_values
def get_compscores_from_files_flatten(all_completed_neut, all_completed_sel, scorestring, sel_bin_index, startbound, endbound, foldDists = False):
	""" for each snp in each replicate, consider all population comparisons. return a single likelihood distribution for the putative selpop.
	n.b., this should correspond to whatever manner of choice is implemented in a given run of combine_scores(_alt)&c """
	if scorestring in ['fst', 'deldaf']:
		physIndex = 0
		selDafIndex = 3
		if scorestring == "fst":
			takeIndex = 2 #counter-intuitive; Fst value is actually in index 1. HOWEVER pass DELDAF to facilitate getting outgroup-outgroup Fst-->PBS.
		elif scorestring == "deldaf":
			takeIndex = 2
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
		selDafIndex = 2
		takeIndex = -1
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

	#THis is not handling pops right????
	neut_values1 = load_from_files_flatten(neut_files1a, neut_files1b, neut_files1c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(neut_values1)) + " neutral values for pop 1 ... (chosen from among three pop comps)")
	neut_values2 = load_from_files_flatten(neut_files2a, neut_files2b, neut_files2c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(neut_values2)) + " neutral values for pop 2 ... (chosen from among three pop comps)")
	neut_values3 = load_from_files_flatten(neut_files3a, neut_files3b, neut_files3c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(neut_values3)) + " neutral values for pop 3 ... (chosen from among three pop comps)")
	neut_values4 = load_from_files_flatten(neut_files4a, neut_files4b, neut_files4c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(neut_values4)) + " neutral values for pop 4 ... (chosen from among three pop comps)")


	causal_values1, linked_values1 = load_from_files_discriminate_causal_flatten(sel_files1a, sel_files1b, sel_files1c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(causal_values1)) + " causal values for pop 1 ... (chosen from among three pop comps)")
	print("loaded " + str(len(linked_values1)) + " linked values for pop 1 ... (chosen from among three pop comps)")
	causal_values2, linked_values2 = load_from_files_discriminate_causal_flatten(sel_files2a, sel_files2b, sel_files2c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(causal_values2)) + " causal values for pop 2 ... (chosen from among three pop comps)")
	print("loaded " + str(len(linked_values2)) + " linked values for pop 2 ... (chosen from among three pop comps)")
	causal_values3, linked_values3 = load_from_files_discriminate_causal_flatten(sel_files3a, sel_files3b, sel_files3c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)
	print("loaded " + str(len(causal_values3)) + " causal values for pop 3 ... (chosen from among three pop comps)")
	print("loaded " + str(len(linked_values3)) + " linked values for pop 3 ... (chosen from among three pop comps)")
	causal_values4, linked_values4 = load_from_files_discriminate_causal_flatten(sel_files4a, sel_files4b, sel_files4c, startbound, endbound, scorestring, stripHeader=True, physIndex=physIndex, absVal = foldDists, takeIndex =takeIndex, selDafIndex=selDafIndex)	
	print("loaded " + str(len(causal_values4)) + " causal values for pop 4 ... (chosen from among three pop comps)")
	print("loaded " + str(len(linked_values4)) + " linked values for pop 4 ... (chosen from among three pop comps)")
	
	all_score_values = [[neut_values1, causal_values1, linked_values1],
						[neut_values2, causal_values2, linked_values2],
						[neut_values3, causal_values3, linked_values3],
						[neut_values4, causal_values4, linked_values4]
						] 
	return all_score_values
def load_from_files(files, startbound, endbound, stripHeader = False, takeIndex = -1, physIndex = 1, absVal = False):
	#print('loading from ' + str(len(files)) + " files...")
	#print("startbound: " + str(startbound) + ", endbound: " + str(endbound))
	values = []
	for file in files:
		openfile = open(file, 'r')
		if stripHeader:
			openfile.readline()
		for line in openfile:
			entries = line.split()
			value = float(entries[takeIndex])
			if absVal:
				value = fabs(value)
			thisPhysPos = int(entries[physIndex])
			if thisPhysPos >= startbound and thisPhysPos <= endbound:
				values.append(value)
		openfile.close()
	return values
def load_posvals_from_files(files, startbound, endbound, stripHeader = False, takeIndex = -1, physIndex = 1, selDafIndex=2, absVal = False):
	#as above but include parallel-indexed physical pos
	pos, values, seldafs = [], [], []
	for file in files:
		openfile = open(file, 'r')
		if stripHeader:
			openfile.readline()
		for line in openfile:
			entries = line.split()
			value = float(entries[takeIndex])
			if absVal:
				value = fabs(value)
			thisPhysPos = int(entries[physIndex])
			thisSelDaf = float(entries[selDafIndex])
			if thisPhysPos >= startbound and thisPhysPos <= endbound:
				values.append(value)
				pos.append(thisPhysPos)
				seldafs.append(thisSelDaf)
		openfile.close()
	return pos, values, seldafs
def load_from_files_discriminate_causal(files,  startbound, endbound, causalLoc = 750000, stripHeader = False, takeIndex = -1, physIndex = 1, absVal = False):
	#print('loading from ' + str(len(files)) + " files...")
	#print("startbound: " + str(startbound) + ", endbound: " + str(endbound))
	causal_values, linked_values = [], []
	for file in files:
		openfile = open(file, 'r')
		if stripHeader:
			openfile.readline()
		for line in openfile:
			entries = line.split()
			value = float(entries[takeIndex])
			if absVal:
				value = fabs(value)
			thisPhysPos = int(entries[physIndex])
			if thisPhysPos >= startbound and thisPhysPos <= endbound:
				if thisPhysPos == causalLoc:
					causal_values.append(value)
				else:
					linked_values.append(value)
		openfile.close()
	return causal_values, linked_values	

#these must go replicate by replicate
def load_from_files_flatten(filesa, filesb, filesc, startbound, endbound, score, stripHeader = False,  takeIndex = -1, physIndex = 1, selDafIndex=2,  absVal = False):
	""" as above, but implements a method to select xpop scores for the same snp(/replicate) from comparisons with multiple populations.
	this must match the method implemented in combine_scores"""
	values = [] 
	nfiles = len(filesa)
	assert(len(filesb) == nfiles)
	assert(len(filesc) == nfiles)
	for ifile in range(nfiles):
		print("neut:" + str(ifile))
		filea = filesa[ifile]
		fileb = filesb[ifile]
		filec = filesc[ifile] 
		pos_a, values_a, seldafs_a = load_posvals_from_files([filea], startbound, endbound, stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		pos_b, values_b, seldafs_b = load_posvals_from_files([fileb], startbound, endbound,  stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		pos_c, values_c, seldafs_c = load_posvals_from_files([filec], startbound, endbound,  stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		all_snps_this_rep = pos_a[:]
		all_snps_this_rep.extend(pos_b)
		all_snps_this_rep.extend(pos_c)
		all_snps_this_rep = set(all_snps_this_rep)
		all_snps_this_rep = list(all_snps_this_rep)
		for snp in all_snps_this_rep:
			thissnp_all_seldafs = []
			availValues = []
			if snp in pos_a:
				index_a = pos_a.index(snp)
				value_a = values_a[index_a]
				seldaf_a = seldafs_a[index_a]
				thissnp_all_seldafs.append(seldaf_a)	
				availValues.append(value_a)			
			if snp in pos_b:
				index_b = pos_b.index(snp)
				value_b = values_b[index_b]
				seldaf_b = seldafs_b[index_b]		
				thissnp_all_seldafs.append(seldaf_b)		
				availValues.append(value_b)			

			if snp in pos_c:
				index_c = pos_c.index(snp)
				value_c = values_c[index_c]
				seldaf_c = seldafs_c[index_c]
				thissnp_all_seldafs.append(seldaf_c)
				availValues.append(value_c)			

			#all_seldafs = [seldaf_a, seldaf_b, seldaf_c]
			all_seldafs = set(thissnp_all_seldafs)
			all_seldafs = list(all_seldafs)
			#assert(len(all_seldafs) == 1)
			if (len(all_seldafs) != 1): # I need to figure out what's wrong with my process here.
				print(thissnp_all_seldafs)		
			value = choose_from_options(availValues, score, all_seldafs[0])
			values.append(value)
	return values
def load_from_files_discriminate_causal_flatten(filesa, filesb, filesc, startbound, endbound, score, causalLoc = 750000, stripHeader = False, takeIndex = -1,  physIndex = 1, selDafIndex=2, absVal = False):
	""" as above, but implements a method to select xpop scores for the same snp(/replicate) from comparisons with multiple populations.
	this must match the method implemented in combine_scores"""
	causal_values, linked_values = [], []
	nfiles = len(filesa)
	assert(len(filesb) == nfiles)
	assert(len(filesc) == nfiles)
	for ifile in range(nfiles):
		print("sel:" + str(ifile))
		filea = filesa[ifile]
		fileb = filesb[ifile]
		filec = filesc[ifile] 
		pos_a, values_a, seldafs_a = load_posvals_from_files([filea], startbound, endbound, stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		pos_b, values_b, seldafs_b = load_posvals_from_files([fileb], startbound, endbound, stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		pos_c, values_c, seldafs_c = load_posvals_from_files([filec], startbound, endbound, stripHeader = stripHeader, physIndex = physIndex,  takeIndex = takeIndex, selDafIndex =selDafIndex, absVal = absVal)
		all_snps_this_rep = pos_a[:]
		all_snps_this_rep.extend(pos_b)
		all_snps_this_rep.extend(pos_c)
		all_snps_this_rep = set(all_snps_this_rep)
		all_snps_this_rep = list(all_snps_this_rep)
		
		for snp in all_snps_this_rep:
			thissnp_all_seldafs = []
			availValues = []
			if snp in pos_a:
				index_a = pos_a.index(snp)
				value_a = values_a[index_a]
				seldaf_a = seldafs_a[index_a]
				thissnp_all_seldafs.append(seldaf_a)
				availValues.append(value_a)
			if snp in pos_b:
				index_b = pos_b.index(snp)
				value_b = values_b[index_b]
				seldaf_b = seldafs_b[index_b]
				thissnp_all_seldafs.append(seldaf_b)
				availValues.append(value_b)
			if snp in pos_c:
				index_c = pos_c.index(snp)
				value_c = values_c[index_c]
				seldaf_c = seldafs_c[index_c]
				thissnp_all_seldafs.append(seldaf_c)
				availValues.append(value_c)

			#all_seldafs = [seldaf_a, seldaf_b, seldaf_c]
			all_seldafs = set(thissnp_all_seldafs)
			all_seldafs = list(all_seldafs)
			#assert(len(all_seldafs) == 1)
			if (len(all_seldafs) != 1):
				print(thissnp_all_seldafs)
			value = choose_from_options(availValues, score, all_seldafs[0]) #Only compare if the snp has a score for this pop-pair!
			if int(snp) == causalLoc: #causal SNP
				causal_values.append(value)
			else: #noncausal
				linked_values.append(value)
	return causal_values, linked_values	
def choose_from_options(values, score, seldaf):
	""" function essential for 'flattening' population comparisons. you have a putative selPop and a bunch of outgroups, so
	how do you use statistics that are defined as a property of population pairs (e.g. XP-EHH, Fst, delDAF)?. this function defines this  
	procedure in PYTHON for the purpose of defining likelihood tables (per-pop) from simulated data. The methods used here for defining
	a given set of likelihood tables should correspond to the actual comparisons being made when calculating empirical CMS values. (These 
	functions are described in C in cms_data.c, but can be toggled from the command line using the python scripts we provide.) """ ##JV DOUBLE CHECK AND MAKE SURE THIS IS FULLY DOCUMENTED

	#### XP-EHH: cf method compareXp()
	#### in cms_data.c
	if score == "xpehh" or score == "xp": #take the maximum, the most positive
		value = max(values)

	#### delDAF: cf method comparedelDaf_outgroup_ave() in cms_data.c
	#### or, take simple ave with comparedelDaf()
	elif score == "deldaf": #compare selpop to worldpop, i.e. set of all outgroups
		if len(values) == 3:
			value_a, value_b, value_c = values
			daf_a = seldaf - value_a #CONFIRM CORRECT 
			daf_b = seldaf - value_b #INCL SIGN
			daf_c = seldaf - value_c
			ave_daf_outgroups = np.mean([daf_a, daf_b, daf_c])
			value = seldaf - ave_daf_outgroups
		elif len(values) == 2:
			value_a, value_b  = values
			daf_a = seldaf - value_a #CONFIRM CORRECT 
			daf_b = seldaf - value_b #INCL SIGN
			ave_daf_outgroups = np.mean([daf_a, daf_b])
			value = seldaf - ave_daf_outgroups
		else:
			value = values[0]


	#### Fst: cf method compareFst_PBS() in cms_data.c
	#### or, take simple ave with compareFst()
	elif score == "fst": #//get PBS for each pair of outgroups with selpop and take the maximum
		## Population-Branch Statistic; an population-specific generalization of Fst for three populations
		## Yi et al., Science 2013
		if len(values) == 3:
			value_a, value_b, value_c = values
			daf_a = seldaf - value_a #BUT THIS NEEDS TO BE DAF NOT FST
			daf_b = seldaf - value_b
			daf_c = seldaf - value_c

			pbs_1 = get_pbs_from_dafs(seldaf, daf_a, daf_b)
			pbs_2 = get_pbs_from_dafs(seldaf, daf_b, daf_c)
			pbs_3 = get_pbs_from_dafs(seldaf, daf_a, daf_c)

			value = max(pbs_1, pbs_2, pbs_3)

		elif len(values) == 2:
			value_a, value_b = values
			daf_a = seldaf - value_a #BUT THIS NEEDS TO BE DAF NOT FST
			daf_b = seldaf - value_b
			pbs = get_pbs_from_dafs(seldaf, daf_a, daf_b)
			value = pbs
		else:
			value = values[0] #just Fst and not PBS, n.b.
			#print('hmm, rethink this.')
			#pass
		#print("VALUE: " + str(value))
	return value


"""
def get_compscores_from_files(all_completed_neut, all_completed_sel, scorestring, sel_bin_index, startbound, endbound, foldDists = False):
	#original version: generates separate likelihood distributions for separate population comparisons 
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

	neut_values1a = load_from_files(neut_files1a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values1b = load_from_files(neut_files1b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values1c = load_from_files(neut_files1c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values2a = load_from_files(neut_files2a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values2b = load_from_files(neut_files2b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values2c = load_from_files(neut_files2c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values3a = load_from_files(neut_files3a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values3b = load_from_files(neut_files3b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values3c = load_from_files(neut_files3c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values4a = load_from_files(neut_files4a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values4b = load_from_files(neut_files4b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	neut_values4c = load_from_files(neut_files4c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values1a, linked_values1a = load_from_files_discriminate_causal(sel_files1a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values1b, linked_values1b = load_from_files_discriminate_causal(sel_files1b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values1c, linked_values1c = load_from_files_discriminate_causal(sel_files1c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values2a, linked_values2a = load_from_files_discriminate_causal(sel_files2a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values2b, linked_values2b = load_from_files_discriminate_causal(sel_files2b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values2c, linked_values2c = load_from_files_discriminate_causal(sel_files2c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values3a, linked_values3a = load_from_files_discriminate_causal(sel_files3a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values3b, linked_values3b = load_from_files_discriminate_causal(sel_files3b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values3c, linked_values3c = load_from_files_discriminate_causal(sel_files3c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values4a, linked_values4a = load_from_files_discriminate_causal(sel_files4a, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values4b, linked_values4b = load_from_files_discriminate_causal(sel_files4b, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
	causal_values4c, linked_values4c = load_from_files_discriminate_causal(sel_files4c, startbound, endbound, stripHeader=True, physIndex=physIndex, absVal = foldDists)
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
"""
