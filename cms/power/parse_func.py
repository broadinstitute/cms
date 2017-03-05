##	functions for manipulating empirical/simulated CMS output
##	last updated 02.13.2017	vitti@broadinstitute.org

import matplotlib as mp 
mp.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os

###################
## COLLATE FILES ##
###################
def write_pair_sourcefile(writefilename, ihsfilename, delihhfilename, nslfilename, xpehhfilename, freqsfilename):
	openfile = open(writefilename, 'w')
	openfile.write(ihsfilename+ "\n")
	openfile.write(delihhfilename+ "\n")
	openfile.write(nslfilename+ "\n")
	openfile.write(xpehhfilename+ "\n")
	openfile.write(freqsfilename+ "\n")
	openfile.close()
	return writefilename
def write_run_paramfile(writefilename, ihs_master_likesfile, nsl_master_likesfile, delihh_master_likesfile, xpehh_master_likesfile,
	fst_master_likesfile, deldaf_master_likesfile, cutoffline, includeline):
	if not os.path.isfile(writefilename):
		openfile = open(writefilename, 'w')
		openfile.write(ihs_master_likesfile + "\n")
		openfile.write(nsl_master_likesfile + "\n") #CHANGE ORDER
		openfile.write(delihh_master_likesfile + "\n")
		openfile.write(xpehh_master_likesfile + "\n")
		openfile.write(fst_master_likesfile + "\n")
		openfile.write(deldaf_master_likesfile + "\n")	
		openfile.write(cutoffline + "\n")
		openfile.write(includeline + "\n")		
		openfile.close()
	return writefilename

##################
## LOCATE FILES ##
##################
def get_emp_component_score_files(filebase = "/n/regal/sabeti_lab/jvitti/clear-synth/1kg_scores/", suffix = "clear-synth-"):
	''' ? version control issue? ''' #NO I just need to patch this in; I swapped to dispatch_1kg_composite.py on RC; being lazy with dispatch_cms_local.py for now but come back to this
	in_ihs_file = basedir + "ihs/chr" + str(chrom) + "_strictMask_" + str(pop) + ".ihs.out"
	in_delihh_file =  basedir + "delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
	in_nsl_file = basedir + "nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
	in_xp_file = basedir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
	in_fst_deldaf_file = basedir + "fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
	return
def get_sim_component_score_files(model, irep, pop, altpop, selbin = "neut", filebase = "/idi/sabeti-scratch/jvitti/clean/scores/", normed = False):
	""" locates component score files for simulated data """
	if selbin == 'neut':
		basedir =  filebase + model + "/neut/"
	else:
		basedir = filebase + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/"
	in_ihs_file = basedir + "ihs/rep" + str(irep) + "_" + str(pop) + ".ihs.out"
	in_delihh_file =  basedir + "delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
	in_nsl_file = basedir + "nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
	in_xp_file = basedir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
	in_fst_deldaf_file = basedir + "fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
	if normed:
		in_ihs_file += ".norm" 
		in_nsl_file += ".norm"
		in_delihh_file += ".norm" 
		in_xp_file += ".norm"
	for returnfile in [in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file]:
		if not os.path.isfile(returnfile):
			print("Missing: " + returnfile)
	return in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file 
def get_concat_files(model, pop, score, altpop = '', basedir = "/idi/sabeti-scratch/jvitti/clean/scores/"):
	""" locates concatenated component score files to facilitate normalization to neutral replicates """
	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_" + str(altpop) + "_"
	else:
		concatfilebase = ""
	concatfilename = concatfilebase + score + ".txt"
	binfilename = concatfilebase + score + ".bins"
	return concatfilename, binfilename
def get_neut_repfile_name(model, irep, pop, normed = False, basedir = "/idi/sabeti-scratch/jvitti/clean/scores/", suffix = ""):
	""" locates neutral simulated composite score file """
	repfilename = basedir +"composite/"+ model + "/neut/rep" + str(irep) + "_" + str(pop) + ".cms" + suffix
	if normed:
		repfilename += ".norm"
	return repfilename
def get_sel_repfile_name(model, irep, pop, freqbin, normed = False, basedir = "/idi/sabeti-scratch/jvitti/scores/", suffix =""):
	""" locate simulated composite score file in scenario with selection """
	repfilename = basedir + "composite/"+ model + "/sel" + str(pop) + "/sel_" + str(freqbin) + "/rep" + str(irep) + "_" + str(pop) + ".cms" + suffix
	if normed:
		repfilename += ".norm"
	return repfilename
def get_likesfiles(model, selpop, likesdir, allfreqs = True, likessuffix= "neut"):
	""" locates master likelihood files for CMS component scores """
	if not allfreqs:
		hi_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_hi_vs_" + likessuffix + ".txt"
		mid_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_mid_vs_" + likessuffix + ".txt"
		low_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_low_vs_" + likessuffix + ".txt"
		return hi_likesfile, mid_likesfile, low_likesfile 
	else:
		likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_allfreq_vs_" + likessuffix + ".txt_2"
		return likesfile
	return
def get_pr_filesnames(key, basedir):
	""" locates files with true/false positive rate data for CMS 2.0 """
	regionlen, percentage, cutoff, model, pop = key
	fprfile = basedir + "/fpr/" + model + "/sel" + str(pop) + "/fpr_" + str(regionlen) + "_" + str(percentage) + "_" + str(cutoff)
	tprfile = basedir + "/tpr/" + model + "/sel" + str(pop) + "/tpr_" + str(regionlen) + "_" + str(percentage) + "_" + str(cutoff)
	for filename in [tprfile, fprfile]:
		if not os.path.isfile(filename):
			print("missing: " + filename)
		else:
			pass
	return fprfile, tprfile 
def get_emp_cms_file(selpop, model, chrom, normed = False, basedir = "/n/regal/sabeti_lab/jvitti/clear-synth/1kg_composite/", suffix = ""):
	""" locates CMS files for empirical data """
	filename = basedir + "chr" + str(chrom) + "_" + str(selpop) + "_strictMask_" + model + ".cms" + suffix
	if normed:
		filename += ".norm"
	if not os.path.isfile(filename):
		print("MISSING empirical file : " + filename)
	return filename

##################
## FILE PARSING ##
##################
def read_cms_repfile(infilename):
	""" loads into memory all data from .cms file """
	physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = [], [], [], [], [], [], [], [], [], [], []
	if not os.path.isfile(infilename):
		return physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed
	openfile = open(infilename, 'r')
	for line in openfile:
		entries = line.split()
		physpos.append(int(entries[0]))
		genpos.append(float(entries[1]))
		seldaf.append(float(entries[2]))
		ihs_normed.append(float(entries[3]))
		delihh_normed.append(float(entries[4]))
		nsl_normed.append(float(entries[5]))
		xpehh_normed.append(float(entries[6]))
		fst.append(float(entries[7]))
		deldaf.append(float(entries[8]))
		cms_unnormed.append(float(entries[9]))
		if ".norm" in infilename:
			cms_normed.append(float(entries[10]))
	openfile.close()
	return physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed
def read_pr(infilename):
	""" reads true/false positive rate from simple file """
	openfile = open(infilename, 'r')
	firstline = openfile.readline()
	pr = float(firstline)
	openfile.close()
	return pr
def read_vals_lastcol(filename):
	""" quick/generic function to pull from last column of a datafile """
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
def load_simscores(model, pop, numRep = 500, normed = False, takeIndex = -1):
	""" for all simulated replicates, loads a value according to takeIndex """
	all_simscores = []
	for irep in range(1, numRep + 1):
		neutfilename = get_neut_repfile_name(model, irep, pop, normed=normed)
		if os.path.isfile(neutfilename):
			#print(neutfilename)
			openfile = open(neutfilename, 'r')
			for line in openfile:
				entries = line.split()
				val = float(entries[takeIndex])
				all_simscores.append(val)
			openfile.close()
	infostring = 'loaded a total of ' + str(len(all_simscores)) 
	if normed:
		infostring += " normed"
	infostring += " scores."
	print(infostring)
	return all_simscores
def load_empscores(model, selpop, normed = False, suffix = '', takeIndex = -1):
	""" for genome-wide empirical CMS data, loads a value according to takeIndex """	
	chroms = range(1,23)
	scores = []
	for chrom in chroms:
		scorefile = get_emp_cms_file(selpop, model, chrom, normed=normed, suffix=suffix)
		print('loading from ' + scorefile)
		assert os.path.isfile(scorefile)
		openfile = open(scorefile, 'r')
		for line in openfile:
			entries=line.split()
			scores.append(float(entries[takeIndex]))
		openfile.close()
	return scores
def load_regions(regionfile):
	""" loads in data from a file with intervals, e.g. denoting regions above significance thresshold """
	openfile = open(regionfile, 'r')
	allchroms, allstarts, allends = [], [], []
	for line in openfile:
		entries = line.split()
		if len(entries) == 3:
			chrom, startpos, endpos = entries[0], int(entries[1]), int(entries[2])
			allchroms.append(chrom)
			allstarts.append(startpos)
			allends.append(endpos)
	openfile.close()
	return allchroms, allstarts, allends

#################
## BOOKKEEPING ##
#################
def execute(commandstring):
	subprocess.check_output(commandstring.split())
	return
def check_create_file(filename, checkOverwrite):
	""" useful to prevent overwriting, but easy to toggle in case one wants to mass-replace files """
	if checkOverwrite == False:
		return True #make it anyway
	else:
		if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
			return True
		else:
			return False
	return False
def check_create_dir(directory):
	""" ensure that the directory exists; create it if it doesn't """
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', '-p', directory])
