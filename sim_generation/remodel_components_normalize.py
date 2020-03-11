##	for each model/regime, MUST GET /THESE/ NEUTRAL NORMALIZATION PARAMETERS.
##	02.25.2019 		# quick hack for def15, use defdef15 norm params
#basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run2/"
basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run4/"

import os
import sys
#import importlib.machinery
#scores_func = importlib.machinery.SourceFileLoader('scores_func','/idi/sabeti-scratch/jvitti/cms/cms/dists/scores_func.py').load_module()
#from scores_func import norm_sel_ihs, norm_sel_xpehh
import subprocess
from math import fabs, sqrt
from random import randint

def read_neut_normfile(neutNormfilename, scoretype ='ihs'):
	"""pulls means and vars so that we can normalize sel scenarios with the same parameters. """
	openfile = open(neutNormfilename)
	line = openfile.readline()
	entries = line.split()
	#while len(entries) < 1 or entries[0] != "Total":
	#	line = openfile.readline()
	#	entries = line.split()
	#nsnp = int(entries[-1])
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

def norm_sel_ihs(inputScoreFile, neutNormfilename):
	''' normalize iHS component score for selection scenarios to neutral parameters ''' 
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
def norm_neut_xpehh(inputScoreFile, outfileName, runProgram = "scans.py"):##JV: I found the way to properly write to file; should implement here
	'''wraps call to scans.py'''
	cmdStr = "python " + runProgram + " selscan_norm_xpehh " + inputScoreFile + " > " + outfileName #return to this
	print(cmdStr)
	return
def norm_sel_xpehh(inputScoreFile, neutNormfilename):
	''' normalize XPEHH component score for selection scenarios to neutral parameters ''' 
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


def get_binfile(score, pop, altpop, model):
	if model == "def15":
		model = "defdef15" #QUICK HACK
	if score != "xpehh":
		binfilename = "/idi/sabeti-scratch/jvitti/remodel/run/normalize_params/" + model + "_" + score + "_OUTGROUP" + str(pop) + ".bins"
	else:
		binfilename = "/idi/sabeti-scratch/jvitti/remodel/run/normalize_params/" + model + "_" + score + "_OUTGROUP" + str(pop) + "_vs" + str(altpop) + ".bins"
	if not os.path.isfile(binfilename):
		zip_binfilename = binfilename +  ".gz"
		if os.path.isfile(zip_binfilename):
			gunzip_cmd = "gunzip " + zip_binfilename
			subprocess.check_output( gunzip_cmd.split() )

	return binfilename

def main():

	model = sys.argv[1]
	regime = sys.argv[2]
	pop = int(sys.argv[3]) 
	irep = int(sys.argv[4])
	#basedir = basewritedir + model + "_" + regime + "_sel" + str(pop) + "/"
	basedir = basewritedir + model + "_" + regime + "/"
	replicate_idstring = "rep" + str(irep) + "_pop"+ str(pop)
	
	#########
	## IHS ##
	#########
	altpop = ""
	ihs_repfilename = basedir + "ihs/" + replicate_idstring + ".ihs.out"
	if not os.path.isfile(ihs_repfilename):
		ihs_repfilename += ".gz"
	if (os.path.isfile(ihs_repfilename) and not os.path.isfile(ihs_repfilename + ".norm")):
		binfilename = get_binfile("ihs", pop,  altpop, model)
	
		norm_sel_ihs(ihs_repfilename, binfilename)
	else:
		print(ihs_repfilename)
	
	############
	## DELIHH ##
	############
	delihh_repfilename =  basedir + "delihh/" + replicate_idstring
	if not os.path.isfile(delihh_repfilename):
		delihh_repfilename += ".gz"

	if (os.path.isfile(delihh_repfilename) and not os.path.isfile(delihh_repfilename + ".norm")):
		binfilename = get_binfile("delihh", pop, altpop, model)
		norm_sel_ihs(delihh_repfilename, binfilename)
	
	#########
	## NSL ##
	#########
	nsl_repfilename = basedir + "nsl/" + replicate_idstring + ".nsl.out"
	if not os.path.isfile(nsl_repfilename):
		nsl_repfilename += ".gz"

	if (os.path.isfile(nsl_repfilename) and not os.path.isfile(nsl_repfilename + ".norm")):
		binfilename = get_binfile("nsl", pop, altpop, model)
		norm_sel_ihs(nsl_repfilename, binfilename)
	
	###########
	## XPEHH ##
	###########
	
	altpops = [1, 2, 3, 4]
	if "gravel" in basedir:
		altpops.remove(4)
	altpops.remove(pop)
	for altpop in altpops:
		xpehh_repfilename = basedir + "xpehh/" + replicate_idstring + "_vs" + str(altpop) + ".xpehh.out"
		if not os.path.isfile(xpehh_repfilename):
			xpehh_repfilename += ".gz"
		if (os.path.isfile(xpehh_repfilename) and not os.path.isfile(xpehh_repfilename + ".norm")):
			binfilename = get_binfile("xpehh", pop, altpop, model)
			print(binfilename)
			norm_sel_xpehh(xpehh_repfilename, binfilename,)
					
					
main()
