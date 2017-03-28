## helper functions for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 03.28.2017 vitti@broadinstitute.org

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
def get_compscores_from_files(all_completed_neut, all_completed_sel, scorestring, sel_bin_index, startbound, endbound, foldDists = False):
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
def load_from_files(files, startbound, endbound, takeIndex = -1, stripHeader = False, physIndex = 1, absVal = False):
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
def load_from_files_discriminate_causal(files,  startbound, endbound, takeIndex = -1, stripHeader = False, causalLoc = 750000, physIndex = 1, absVal = False):
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
