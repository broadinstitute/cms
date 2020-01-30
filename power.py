## script to manipulate and analyze empirical/simulated CMS output
## last updated 11.20.16

from power_parser import full_parser_power
from power_func import normalize, merge_windows, get_window, check_outliers, check_rep_windows, calc_pr, get_pval, plotManhattan, \
						plotManhattan_extended, loadregions, quick_plot, get_causal_rank, get_cdf_from_causal_ranks
from parse_func import get_component_score_files, get_neut_repfile_name, get_sel_repfile_name, get_emp_cms_file, read_cms_repfile, load_simscores, \
						load_empscores, writeJobsToTaskArray_fromArglist, check_create_dir, execute, read_pr, write_master_likesfile, get_likesfiles, \
						get_neutinscorefiles, check_create_file, get_concat_files, plot_dist, readvals_lastcol

import matplotlib as mp 
mp.use('agg')
import matplotlib.pyplot as plt
from tempfile import TemporaryFile
#from xlwt import Workbook, easyxf 
#from pybedtools import BedTool
import numpy as np
import argparse
import sys
import os

#################################
## DEFINE EXECUTIVE FUNCTIONS ###
#################################

########	Manipulate simulated data 
########	and provide composite scores.
def execute_run_neut_sims(args):
	'''from run_additional_sims.py'''
	cmd = "coalescent"
	model = args.model
	nrep = args.nrep
	nStart = args.startrep
	nEnd = nStart + nrep
	writefolder = args.writedir
	paramfolder = args.paramfolder
	paramfilename = paramfolder + model + ".par"
	
	for irep in range(nStart, nEnd + 1):
		outbase = writefolder + "rep" + str(irep)
		argstring = "-p " + paramfilename + " --genmapRandomRegions --drop-singletons .25 --tped " + outbase
		cosicreatefilename = outbase + "_0_1.tped"

		proceed = check_create_file(cosicreatefilename, args.checkOverwrite)
		if proceed:
			fullCmd = cmd + " " + argstring
			#print(fullCmd)
			execute(fullCmd)
			for ipop in [1, 2, 3, 4]:
				torenamefile = outbase + "_0_" + str(ipop) + ".tped"
				#print(torenamefile)
				assert os.path.isfile(torenamefile)
				renamed = outbase +"_" + str(ipop) + ".tped"
				renamecmd = "mv "  + torenamefile + " " + renamed
				execute(renamecmd)
	print("wrote simulates to e.g. " + renamed)
	return
def execute_run_neut_repscores(args):
	''' adapted from rerun_scores_onerep.py'''
	model = args.model
	pop = args.simpop
	repNum = args.irep
	basedir = args.writedir
	cmsdir = args.cmsdir
	tpeddir = args.tpedfolder
	simRecomFile = args.simRecomFile

	tped = tpeddir + "rep" + str(repNum) + "_" + str(pop) + ".tped"
	assert os.path.isfile(tped)
	#in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file  = get_component_score_files(model, repNum, pop, altpop, scenario = "neut", filebase = basedir)
		
	####### Calculate per-population
	####### scores: iHS, delIHH, nSL
	ihs_commandstring = "python " + cmsdir + "scans.py selscan_ihs"
	ihs_outfileprefix = basedir + "ihs/rep" + str(repNum) + "_" + str(pop) 
	ihs_unnormedfile = ihs_outfileprefix + ".ihs.out"
	ihs_argstring = tped + " " + ihs_outfileprefix + " --threads 7 --truncOk"
	ihs_fullcmd = ihs_commandstring + " " + ihs_argstring
	proceed = check_create_file(ihs_unnormedfile, args.checkOverwrite)
	if proceed:
		print(ihs_fullcmd)
		execute(ihs_fullcmd)
	delihh_commandstring = "python " + cmsdir + "likes_from_model.py scores_from_sims --delIhh"
	delihh_unnormedfile =  basedir + "delihh/rep" + str(repNum) + "_" + str(pop) + ".txt"
	delihh_argstring = ihs_unnormedfile + " "+ delihh_unnormedfile
	delihh_fullcmd = delihh_commandstring + " " + delihh_argstring
	proceed = check_create_file(delihh_unnormedfile, args.checkOverwrite)
	if proceed:
		print(delihh_fullcmd)
		execute(delihh_fullcmd)		
	nsl_commandstring = "python " + cmsdir + "scans.py selscan_nsl" 
	nsl_unnormedfileprefix = basedir + "nsl/rep" + str(repNum) + "_" + str(pop)
	nsl_argstring = tped + " " + nsl_unnormedfileprefix
	nsl_fullcmd = nsl_commandstring + " " + nsl_argstring
	nsl_unnormedfilename = nsl_unnormedfileprefix + ".nsl.out"
	proceed = check_create_file(nsl_unnormedfilename, args.checkOverwrite)
	if proceed:
		print(nsl_fullcmd)
		execute(nsl_fullcmd)	

	####### Calculate per-population-pair
	####### scores: XP-EHH, Fst, delDAF
	pops = [1, 2, 3, 4]
	altpops = pops[:]
	altpops.remove(int(pop))
	for altpop in altpops:
		xpehh_commandstring = "python " + cmsdir + "scans.py selscan_xpehh --threads 7"
		tped2 = tpeddir + "rep" + str(repNum) + "_" + str(altpop) + ".tped"
		xpehh_outfileprefix = basedir + "xpehh/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop)
		xpehh_unnormedfile = basedir + "xpehh/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
		xpehh_argumentstring = tped + " " + xpehh_outfileprefix + " " + tped2
		xpehh_fullcmd = xpehh_commandstring + " " + xpehh_argumentstring
		proceed = check_create_file(xpehh_unnormedfile, args.checkOverwrite)
		if proceed:
			print(xpehh_fullcmd)
			execute(xpehh_fullcmd)	
		fstdeldaf_commandstring = "python " + cmsdir + "likes_from_model.py scores_from_sims"
		fstdeldaf_outfilename = basedir + "fst_deldaf/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop)
		fstdeldaf_argumentstring = tped + " --fst_deldaf " + tped2 + " " + fstdeldaf_outfilename + " --recomfile " + simRecomFile
		fstdeldaf_fullcmd = fstdeldaf_commandstring + " " + fstdeldaf_argumentstring
		proceed = check_create_file(fstdeldaf_outfilename, args.checkOverwrite)
		if proceed:
			print(fstdeldaf_fullcmd)
			execute(fstdeldaf_fullcmd)
	return
def execute_run_norm_neut_repscores(args):
	''' creates a concatenated file and saves bin output if run for the first time;
		if these already exist, run per-rep'''
	model = args.model
	pop = args.simpop
	numReps = args.nrep
	basedir = args.writedir
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop

	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		altpop = args.altpop
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_" + str(altpop) + "_"

	concatfilename = concatfilebase + score + ".txt"
	binfilename = concatfilebase + score + ".bins"

	if not os.path.isfile(binfilename):
		if not os.path.isfile(concatfilename):
			repfiles = []
			for irep in range(1, numReps+1):
				if score == 'ihs':
					unnormedfile = basedir + model + "/neut/ihs/rep" + str(irep) + "_" + str(pop) + ".ihs.out"
				elif score == "delihh":
					unnormedfile = basedir + model + "/neut/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
				elif score == "nsl":
					unnormedfile = basedir + model + "/neut/nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
				elif score == "xpehh":
					unnormedfile = basedir + model + "/neut/xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
				elif score == "fst":
					unnormedfile = basedir + model + "/neut/fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
				if os.path.isfile(unnormedfile):
					repfiles.append(unnormedfile)
				else:
					print('missing: ' + unnormedfile)
			concatfile = open(concatfilename, 'w')
			for irepfile in range(len(repfiles)):
				repfile = repfiles[irepfile]
				readfile = open(repfile, 'r')
				firstline = readfile.readline()
				if score in ['xpehh', 'fst']: #header
					if irepfile == 0:
						concatfile.write(firstline)
					else:
						pass
				else:
					concatfile.write(firstline)
				for line in readfile:
					concatfile.write(line)
				readfile.close()
			concatfile.close()
			print('wrote to: ' + concatfilename)

		#already have concatfile
		infilename = concatfilename
		outfilename = concatfilename + ".norm"
		argstring = infilename #+ " " + outfilename

		if score in ['ihs', 'delihh']:
			cmd = "python " + args.cmsdir + "scans.py selscan_norm_ihs"
		elif score in ['nsl']:
			cmd = "python " + args.cmsdir + "scans.py selscan_norm_nsl"
		elif score in ['xpehh']:
			cmd = "python " + args.cmsdir + "scans.py selscan_norm_xpehh"
		else:
			cmd = ""

		fullcmd = cmd + " " + argstring
		execute(fullcmd)
		print("use rename_binfiles.py")
		#print(fullcmd)
	return
def execute_norm_from_binfile(args):
	model = args.model
	pop = args.simpop
	numReps = args.nrep
	basedir = args.writedir
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop
	concatfilename, binfilename = get_concat_files(model, pop, score, altpop)
	
	for irep in range(1, numReps+1):
		if score in ['ihs']:
			unnormedfile = basedir  + model + "/neut/ihs/rep" + str(irep) + "_" + str(pop)  + ".ihs.out"
			normedfilename = unnormedfile + ".norm"
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['delihh']:
			unnormedfile = basedir + model + "/neut/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
			normedfilename = unnormedfile + ".norm"	
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['nsl']:
			unnormedfile = basedir + model + "/neut/nsl/rep" + str(irep) + "_" + str(pop)+ ".nsl.out"
			normedfilename = unnormedfile + ".norm"			
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['xpehh']:
			unnormedfile = basedir + model + "/neut/xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
			normedFilename = unnormedfile + ".norm"
			argstring = "--normalizeXpehh --neutXpehhNormParams " + binfilename + " " + unnormedfile + " " + normedFilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['fst']:
			print('currently handling this manually: rewrite_fst_bins.py')
			pass
		fullcmd = commandstring + " " + argstring
		#print(fullcmd)
		execute(fullcmd)
	return
def execute_run_poppair(args):
	''' from run_additioanal_poppair.py'''
	model = args.model
	selpop = args.simpop
	altpop = args.altpop
	repNum = args.nrep
	writedir = args.writedir #just for poppairs
	cmsdir = args.cmsdir
	cmd = "python " + cmsdir + "composite.py poppair"

	modeldir = writedir + model + "/"
	check_create_dir(modeldir)
	neutdir = modeldir + "neut/"
	check_create_dir(neutdir)

	for irep in range(1, repNum+1):
		in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_component_score_files(model, irep, selpop, altpop, "neut")
		outfile = neutdir + "pairs/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		
		alreadyExists = False
		if args.checkOverwrite:
			if not os.path.isfile(outfile): #check for overwrite
				alreadyExists = False
			else:
				alreadyExists = True				
		if alreadyExists == False:
			argstring = in_ihs_file + " " + in_nsl_file + " " + in_delihh_file + " " + in_xp_file + " " + in_fst_deldaf_file + " " + outfile
			fullcmd = cmd + " " +  argstring
			#print(fullcmd)
			execute(fullcmd)
	return
def execute_composite_sims(args):
	model = args.model
	selpop = args.simpop
	likessuffix = args.likessuffix
	likesdir = args.likes_basedir
	cmsdir = args.cmsdir
	writedir = args.writedir
	numPerBin = args.nrep

	cmd = "python " + cmsdir + "composite.py outgroups"

	hi_likesfile = get_likesfiles(model, selpop, likessuffix, likesdir, allfreqs=True)
	mid_likesfile, low_likesfile = hi_likesfile, hi_likesfile #not using likesfreqs for now


	#ALL SEL SIMS
	"""
	sel_freq_bins = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
	for sel_freq_bin in sel_freq_bins:
		writeselfreqbin = writeseldir + "/sel_" + str(sel_freq_bin)
		check_create_dir(writeselfreqbin)
		for irep in range(1, numPerBin +1):
			inscorefilelist = get_inscorefiles(model, selpop, sel_freq_bin, irep, pairbasedir) 
			if len(inscorefilelist) > 0:
				scorefilelist = ""
				for filename in inscorefilelist:
					scorefilelist += filename + ","
				scorefilelist = scorefilelist[:-1]
				outfile = basedir + model + "/sel" + str(selpop) + "/sel_" + str(sel_freq_bin) + "/rep" + str(irep) + "_" + str(selpop) +"_vsneut.cms.out"
				alreadyExists = False
				if args.checkOverwrite:
					if not os.path.isfile(outfile): #check for overwrite
						alreadyExists = False
					else:
						alreadyExists = True				
				if alreadyExists == False:
					argstring = scorefilelist + " " + hi_likesfile + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + str(selpop) + " " + outfile 
					fullcmd = cmd + argstring
					print(fullcmd)
					#execute(fullcmd)
		"""
	#ALL NEUT
	for irep in range(1, numPerBin +1):	
		neutinscorefilelist = get_neutinscorefiles(model, selpop, irep)#, pairbasedir) 
		if len(neutinscorefilelist) > 0:
			scorefilelist = ""
			for filename in neutinscorefilelist:
				scorefilelist += filename + ","
			scorefilelist = scorefilelist[:-1]
			outfile = writedir  + "rep" + str(irep) + "_" + str(selpop) +".cms.out"
			
			alreadyExists = False
			if args.checkOverwrite:
				if not os.path.isfile(outfile): #check for overwrite
					alreadyExists = False
				else:
					alreadyExists = True				
			if alreadyExists == False:
				argstring = scorefilelist + " " + hi_likesfile + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + str(selpop) + " " + outfile
				fullcmd = cmd + " " + argstring
				#print(fullcmd)
				execute(fullcmd)
	return
def execute_normsims(args):
	model = args.model
	selpop = args.simpop
	numPerBin = args.nrep
	if args.likessuffix == "linked":
		vsNeut = False
	else:
		vsNeut = True

	sel_freq_bins = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']

	values = []
	##############################
	## LOAD STATS FROM NEUT SIMS #
	##############################
	for irep in range(1, numPerBin +1):	
		outfile  = get_neut_repfile_name(model, irep, selpop, vsNeut)
		#print(outfile)
		if os.path.isfile(outfile):
			openfile = open(outfile, 'r')
			for line in openfile:
				entries = line.split()
				rawscore = np.log(float(entries[-1]))
				values.append(rawscore)
			openfile.close()

	print('loaded ' + str(len(values)) + ' values from neutral sims...')

	#check for nans
	values = np.array(values)
	values = values[~np.isnan(values)]
	mean = np.mean(values)
	var = np.var(values)
	sd = np.sqrt(var)

	print("max: " + str(max(values)))
	print("min: " + str(min(values)))
	print("mean: " + str(np.mean(values)))
	print("var: " + str(np.var(values)))

	###############
	## NORMALIZE ##
	###############
	#ALL NEUT
	
	for irep in range(1, numPerBin +1):	

		outfile  = get_neut_repfile_name(model, irep, selpop, vsNeut,)

		if os.path.isfile(outfile):
			normedfile = outfile + ".norm"
			if not os.path.isfile(normedfile):
				openfile = open(outfile, 'r')
				writefile = open(normedfile, 'w')
				for line in openfile:
					entries = line.split()
					rawscore = np.log(float(entries[-1]))
					normalized = normalize(rawscore, mean, sd)
					writeline = line.strip('\n') + "\t" + str(normalized)+ "\n"
					writefile.write(writeline)
				openfile.close()
				writefile.close()
	print("wrote to eg: " + normedfile)	

	#ALL SEL SIMS
	#for sel_freq_bin in sel_freq_bins:
	#	for irep in range(1, numPerBin +1):
	#		rawfile = get_sel_repfile_name(model, irep, selpop, sel_freq_bin, vsNeut = True, normed = False)
	#		if os.path.isfile(rawfile):
	#			normedfile = rawfile + ".norm"
	#			if not os.path.isfile(normedfile):
	#				openfile = open(rawfile, 'r')
	#				writefile = open(normedfile, 'w')
	#				for line in openfile:
	#					entries = line.split()
	#					rawscore = np.log(float(entries[-1]))
	#					normalized = normalize(rawscore, mean, sd)
	#					writeline = line.strip('\n') + "\t" + str(normalized) + "\n"
	#					writefile.write(writeline)
	#				openfile.close()
	#				writefile.close()
	#print("wrote to eg: " + normedfile)	

	return
def execute_write_master_likes(args):
	''' from write_master_likes.py'''
	basedir = args.likes_basedir
	models = ['nulldefault_constantsize', 'default_112115_825am', 'gradient_101915_treebase_6_best', 'nulldefault', 'default_default_101715_12pm']
	selpops = [1, 2, 3, 4]
	freqs = ['allfreq']#'hi', 'low', 'mid']
	misses = ["neut"]#, "linked"]
	for model in models:
		for selpop in selpops:
			for freq in freqs:
				for miss in misses:					
					writefiledir = basedir + model + "/master/"
					check_create_dir(writefiledir)
					writefilename = writefiledir + "likes_" + str(selpop) + "_" + str(freq) + "_vs_" + str(miss) + ".txt"
					write_master_likesfile(writefilename, model, selpop, freq, basedir, miss)
	return

########	Manipulate empirical data
########	and provide composite scores
def execute_composite_emp(args):
	''' from empirical_composite.py '''
	cmsdir = args.cmsdir
	basedir = args.basedir
	likesdir = args.likes_basedir
	likessuffix = args.likessuffix
	model = args.model
	selPop = args.emppop
	modelPop = args.simpop
	cmd = "python " + cmsdir + "composite.py outgroups"
	chroms = range(1,23)

	hi_likesfile = get_likesfiles(model, modelPop, likessuffix, likesdir, allfreqs=True)
	mid_likesfile, low_likesfile = hi_likesfile, hi_likesfile
	for chrom in chroms:
		inscorefilelist = get_pairfiles(chrom, selpop) 
		if len(inscorefilelist) > 0:
			scorefilelist = ""
			for filename in inscorefilelist:
				scorefilelist += filename + ","
			scorefilelist = scorefilelist[:-1]
			outfile = basedir + selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt"
			alreadyExists = False
			if args.checkOverwrite:
				if not os.path.isfile(outfile): #check for overwrite
					alreadyExists = False
				else:
					alreadyExists = True				
			if alreadyExists == False:
				argstring = scorefilelist + " " + hi_likesfile + " " + str(modelPop) + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + outfile
				fullcmd = cmd + argstring
				print(fullcmd)
				#execute(fullcmd)	
def execute_normemp(args):
	selpop = args.emppop
	model = args.model
	likessuffix = args.likessuffix
	
	scores2 = load_empscores(model, selpop, likessuffix, normed=False)
	scores = [np.log(item) for item in scores2]

	#check for nans
	scores = np.array(scores)
	scores = scores[~np.isnan(scores)]

	print('loaded ' + str(len(scores)) + " scores")
	print("max: " + str(max(scores)))
	print("min: " + str(min(scores)))
	print("mean: " + str(np.mean(scores)))
	print("var: " + str(np.var(scores)))

	mean = np.mean(scores)
	var = np.var(scores)
	sd = np.sqrt(var)

	##############
	## NORMALIZE #
	##############
	chroms = range(1,23)
	for chrom in chroms:
		unnormedfile = get_emp_cms_file(selpop, model, likessuffix, chrom, normed=False)
		assert os.path.isfile(unnormedfile)
		normedfile = unnormedfile + ".norm"

		readfile = open(unnormedfile, 'r')
		writefile = open(normedfile, 'w')
		for line in readfile:
			line = line.strip('\n')
			entries=line.split()
			rawscore = np.log(float(entries[-1]))
			normedscore = normalize(rawscore, mean, sd)
			writeline = line + "\t" + str(normedscore) + '\n'
			writefile.write(writeline)
		readfile.close()
		writefile.close
		print('wrote to '  + normedfile)

	return

########	Visualize composite score 
########	output for a given CMS run.
def execute_repviz(args):
	'''visualize CMS scores for simulated data as individual replicates (w/ component and combined scores vs. pos)'''
	causal_index = -1
	model, pop = args.model, args.simpop
	reps = args.nrep
	savefilename = args.savefilename
	scenars = ['0.50']#'neut', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']

	for scenar in scenars:
		for irep in range(1, reps+1):
			if scenar == "neut":
				cmsfilename = get_neut_repfile_name(model, irep, pop, normed = True)
			else:
				cmsfilename = get_sel_repfile_name(model, irep, pop, scenar, normed =True, vsNeut = True)
				#print(cmsfilename)
			if os.path.isfile(cmsfilename):
				print('loading from... ' + cmsfilename)
				physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(cmsfilename)

				f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, sharex = True)
				quick_plot(ax1, physpos, ihs_normed, "ihs_normed", causal_index)
				quick_plot(ax2, physpos, delihh_normed, "delihh_normed", causal_index)
				quick_plot(ax3, physpos, xpehh_normed, "xpehh_normed", causal_index)
				quick_plot(ax4, physpos, fst, "fst", causal_index)
				quick_plot(ax5, physpos, deldaf, "deldaf", causal_index)
				quick_plot(ax6, physpos, cms_unnormed, "cms_unnormed", causal_index)
				quick_plot(ax7, physpos, cms_normed, "cms_normed", causal_index)				
				savefilename = args.savefilename + "_"+ model + "_" + str(pop) + "_" + str(scenar) + "_" + str(irep) + ".png"

				plt.savefig(savefilename)
				if irep == 1:
					print("plotted to " + savefilename)
				plt.close()
	return
def execute_distviz(args):
	"""from explore_p3_dists.py"""
	print("MUST UPDATE")
	############################
	## VIEW DIST FROM 1 FILE  ##
	############################
	oneFile = False
	simFiles = True
	models = [ 'gradient_101915_treebase_6_best', 'default_default_101715_12pm', 'default_112115_825am','nulldefault_constantsize', 'nulldefault', ]
	modelpops = [1, 2, 3, 4]

	if oneFile:
		savefilename = "/web/personal/vitti/test2.png"
		infilename = "/idi/sabeti-scratch/jvitti/scores_composite_5/nulldefault/neut/rep1000_4.cms.out"
		#"/idi/sabeti-scratch/jvitti/scores_composite_5/nulldefault/neut/rep1000_4.cms.out"
	#"/idi/sabeti-scratch/jvitti/scores_composite4_b/CHB_gradient_101915_treebase_6_best_neut.chr10.txt.norm"
		#'/idi/sabeti-scratch/jvitti/scores_composite4_b/default_112115_825am/neut/rep500_2.cms.out'
		#'/idi/sabeti-scratch/jvitti/scores_composite4/default_112115_825am/sel1/sel_0.10/rep47_1_vsneut.cms.out'
		if os.path.isfile(infilename):
			values = readvals_lastcol(infilename)
			plot_dist(values, savefilename)


	#########################
	## VIEW DIST FROM SIMS ##
	#########################
	elif simFiles:
		for model in models:		
			for pop in modelpops:
				savefilename = "/web/personal/vitti/" + model +"_" + str(pop) + ".png"
				allvals = []

				#NEUT REPS
				for irep in range(1,1001):
					filebase = "/idi/sabeti-scratch/jvitti/clean/scores/"
					infilename = filebase + model + "/neut/composite/rep" + str(irep) + "_" + str(pop) + ".cms.out"
					#print(infilename)
					if os.path.isfile(infilename):
						values = readvals_lastcol(infilename)
						allvals.extend(values)
				#print(allvals)
				plot_dist(allvals, "/web/personal/vitti/" + model +"_" + str(pop) + "_justneut.png")

				#SEL REPS
				#for selbin in selbins:	
				#	for irep in range(1,501):
				#		infilename = filebase + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/rep" + str(irep) + "_" + str(pop) + "_vsneut.cms.out"
				#		if os.path.isfile(infilename):
				#			values = readvals_lastcol(infilename)
				#			allvals.extend(values)
			
				#plot_dist(allvals, savefilename)

	#############################
	## VIEW DIST FROM EMP DATA ##
	#############################
	elif empFiles:
		for model in models:
			for pop in emppops:
				savefilename = "/web/personal/vitti/" + model +"_" + str(pop) + ".png"
				allvals = []
				for chrom in chroms:
					infilename = filebase + pop + "_" + model + "_" + "neut" + ".chr" + str(chrom) + ".txt"
					if os.path.isfile(infilename):
						values = readvals_lastcol(infilename)
						allvals.extend(values)
				plot_dist(allvals, savefilename)



	"""
	''' view distribution of (currently) unnormalized cms_gw scores '''
	model = args.model
	savename = args.savefilename
	if args.normed_cms:
		normed = True
	else:
		normed = False
	if args.sim:
		pop = args.simpop
		values = load_simscores(model, pop, normed=normed)
	elif args.emp:
		pop = args.emppop
		likessuffix = args.likessuffix
		values = load_empscores(model, pop, likessuffix, normed=normed)

	if normed==False:
		print('taking log of raw scores...')
		#print(values[0])
		logvalues = [np.log(item) for item in values if not np.isnan(item)]
		#print(logvalues[0])
		values = logvalues

	f, ax = plt.subplots(1)
	ax.hist(values, bins=1000)
	plt.savefig(savename)
	print('plotted to ' + savename)
	"""
	return

########	Quantify and visualize power
########	across significance cutoffs.
def execute_cdf(args):
	reps = args.nrep
	savefilename = args.savefilename
	scenars = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
	model = args.model
	causal_ranks = []
	for pop in [1, 2, 3, 4]:
		for scenar in scenars:
			for irep in range(1, reps+1):
				cmsfilename = get_sel_repfile_name(model, irep, pop, scenar, normed =True, vsNeut = True)
				if os.path.isfile(cmsfilename):
					physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(cmsfilename)
					
					causalPos = 500000
					if causalPos in physpos:
						causal_index = physpos.index(causalPos)
						causal_unnormed = cms_unnormed[causal_index]
						causal_rank = get_causal_rank(cms_unnormed, causal_unnormed)
						#print('causal rank: ' + str(causal_rank)) 
						causal_ranks.append(causal_rank)

	cdf_fig, cdf_ax = plt.subplots()
	cdf_bins, cdf = get_cdf_from_causal_ranks(causal_ranks)
	cdf_ax.plot(cdf_bins[1:], cdf)
	cdf_ax.set_xlim([0, 50])
	plt.title(model + ", " + str(len(causal_ranks)) + " selection replicates")
	plt.ylabel('probability that the causal variant is captured')
	plt.xlabel('significance thresshold (i.e., examining the top x variants)')
	plt.savefig(savefilename)
	plt.close()
	print('plotted to ' + savefilename)
	return
def execute_fpr(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	numReps = args.nrep
	pop = args.simpop

	all_scores = []
	all_percentages = []
	
	if True:
		for irep in range(1, numReps + 1):
			repfilename = get_neut_repfile_name(model, irep, pop, normed=True)
			physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
			if len(cms_normed) > 0:
				all_scores.append(cms_normed)
				rep_percentages = check_rep_windows(physpos, cms_normed, regionlen, cutoff = cutoff)
				all_percentages.append(rep_percentages)		

	print('loaded ' + str(len(all_scores)) + " replicates populations for model " + model + "...")
	fpr = calc_pr(all_percentages, thresshold)
	print('false positive rate: ' + str(fpr) + "\n")

	if args.saveLog	is not None:
		writefilename = args.saveLog
		writefile = open(writefilename, 'w')
		writefile.write(str(fpr)+'\n')

		writefile.write(model + "\t" + str(regionlen) + "\t" + str(thresshold) + '\t' + str(cutoff) + '\n')
		writefile.close()
		print('wrote to :  ' + str(writefilename))
	return
def execute_tpr(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	numReps = args.nrep
	pop = args.simpop

	all_scores = []
	all_percentages = []
	
	if args.saveLog	is not None:
		writefilename = args.saveLog
		if os.path.isfile(writefilename):
			print(writefilename + " already exists; aborting.")
			sys.exit(0)

	allrepfilenames = []
	for selbin in ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']:
		for irep in range(1, numReps + 1):
			repfilename = get_sel_repfile_name(model, irep, pop, selbin, normed=True, vsNeut=True)
			if os.path.isfile(repfilename):
				allrepfilenames.append(repfilename)
	chosen = np.random.choice(allrepfilenames, 1000, replace=False) #take random sample	
	for repfilename in chosen:
		physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
		if len(cms_normed) > 0:
			all_scores.append(cms_normed)
			rep_percentages = check_rep_windows(physpos, cms_normed, regionlen, cutoff = cutoff)
			all_percentages.append(rep_percentages)		

	print('loaded ' + str(len(all_scores)) + " replicates populations for model " + model + "...")
	tpr = calc_pr(all_percentages, thresshold)
	print('true positive rate: ' + str(tpr) + "\n")

	if args.saveLog	is not None:
		writefilename = args.saveLog
		writefile = open(writefilename, 'w')
		writefile.write(str(tpr)+'\n')

		writefile.write(model + "\t" + str(regionlen) + "\t" + str(thresshold) + '\t' + str(cutoff) + '\n')
		writefile.close()
		print('wrote to :  ' + str(writefilename))
	return	
def execute_roc(args):
	''' from quick_roc.py '''
	plot_roc = args.plot_curve
	find_opt = args.find_opt
	maxFPR = args.maxFPR
	fprloc = args.fprloc
	tprloc = args.tprloc

	regionlens = [10000, 25000, 50000, 100000]
	thressholds = [25, 30, 35, 40, 45, 50]
	cutoffs = [2, 2.5, 3, 3.5, 4, 4.5, 5]
	models = ['default_112115_825am', 'nulldefault', 'nulldefault_constantsize', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm']#
	pops = [1, 2, 3, 4]

	###############
	## LOAD DATA ##
	###############
	allfpr = {}
	alltpr = {}
	for regionlen in regionlens:
		for percentage in thressholds:
			for cutoff in cutoffs:
				for model in models:
					allpops_fpr, allpops_tpr = [], []
					for pop in pops:
						this_key = (regionlen, percentage, cutoff, model, pop)
						fprfile = fprloc + "fpr_111716_" + model + "_" + str(pop) + "_" + "neut" + "_" + str(regionlen) +"_" + str(percentage) + "_"+ str(cutoff) +".txt"
						tprfile = tprloc + "tpr_111716_" + model + "_" + str(pop) + "_" + "neut" + "_" + str(regionlen) +"_" + str(percentage) + "_"+ str(cutoff) +".txt"
						#print(tprfile)
						if os.path.isfile(fprfile) and os.path.isfile(tprfile) and os.path.getsize(fprfile) > 0 and os.path.getsize(tprfile) > 0:
							fpr = read_pr(fprfile)
							tpr = read_pr(tprfile)
							allfpr[this_key] = fpr
							alltpr[this_key] = tpr
							allpops_fpr.append(fpr)
							allpops_tpr.append(tpr)
						else:
							print("missing " + tprfile + "\t" + fprfile)
					#assert len(allpops_fpr) == 4
					ave_fpr = np.average(allpops_fpr)
					#assert len(allpops_tpr) == 4
					ave_tpr = np.average(allpops_tpr)

					popave_key = (regionlen, percentage, cutoff, model, "ave")
					allfpr[popave_key] = ave_fpr
					alltpr[popave_key] = ave_tpr

	###############
	## PLOT DATA ##
	###############
	if plot_roc:
		fig, ax = plt.subplots(1)
		for model in models:
			#need to plotfpr from dict
			plotfpr, plottpr = [], []
			for regionlen in regionlens:
				for percentage in thressholds:
					for cutoff in cutoffs:
						this_key = (regionlen, percentage, cutoff, model, "ave")
						plotfpr.append(allfpr[this_key])
						plottpr.append(alltpr[this_key])

			if (len(allfpr)) > 0:
				#ax.scatter(plotfpr, plottpr, label=model)
				plotfpr, plottpr = zip(*sorted(zip(plotfpr, plottpr)))
				ax.plot(plotfpr, plottpr, label=model)

		ax.set_xlabel('FPR')
		ax.set_ylabel('TPR')
		ax.set_xlim([0,.1])
		ax.set_ylim([0,1])
		savefilename = "/web/personal/vitti/roc/" +"compare.png"#+ model + ".png"
		plt.legend()
		plt.savefig(savefilename)
		print("plotted to " + savefilename)
		plt.close()
			
	#####################
	## COMPARE CUTOFFS ##
	#####################
	if find_opt:
		for model in models:
			print(model)
			for regionlen in regionlens:
				for percentage in thressholds:
					for cutoff in cutoffs:
						avekey = (regionlen, percentage, cutoff, model, "ave")
						#FPR passes thresshold. What's the best TPR we can get?
						if allfpr[avekey] <= maxFPR:
							print(avekey)
							print(str(alltpr[avekey]))

	

	return

########	Apply significance cutoff 
########	to empirical results.
def execute_gw_regions(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	pop = args.emppop
	windowlen = args.regionlen
	vs = args.likessuffix

	chroms = range(1,23)
	signif_windows = []
	####################
	## LOOP OVER CHRS ##
	####################
	for chrom in chroms:
		chrom_signif = []
		normedempfilename = get_emp_cms_file(pop, model, vs, chrom, normed=True)
		if not os.path.isfile(normedempfilename):
			print("missing: " + normedempfilename)
		else:
			physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(normedempfilename)
			for iPos in range(len(physpos)):
				##################
				## CHECK REGION ##
				##################
				window_scores = get_window(iPos, physpos, cms_normed, windowlen)
				percentage = check_outliers(window_scores, cutoff)
				if percentage > thresshold:
					chrom_signif.append(physpos[iPos])
		signif_windows.append(chrom_signif)

	##############################
	## MERGE CONTIGUOUS WINDOWS ##
	##############################
	final_starts = []
	final_ends = []
	print('merging regions')
	for chrom_signif in signif_windows:
		starts, ends =  merge_windows(chrom_signif, windowlen)
		final_starts.append(starts)
		final_ends.append(ends)
	#print(final_starts)
	#print(final_ends)

	###################
	## WRITE TO FILE ##
	###################
	if args.saveLog is not None:
		writefilename = args.saveLog 
		writefile = open(writefilename, 'w')
		for ichrom in range(len(final_starts)):
			chromnum = ichrom + 1

			starts = final_starts[ichrom]
			ends = final_ends[ichrom]

			for iregion in range(len(starts)-1):
				writeline = "chr" + str(chromnum) + "\t" + str(starts[iregion]) + "\t" + str(ends[iregion]) + '\n'
				writefile.write(writeline)
		writefile.close()
		print('wrote to ' + writefilename)
	
	return
def execute_regionlog(args):
	''' writes an excel file '''
	regionfiles = []
	pops = ['YRI', 'CEU', 'CHB', 'BEB']
	takepops = []
	for pop in pops:
		#regionfilename = "/idi/sabeti-scratch/jvitti/cms2_power/regions/regions_103116_" + str(pop) +"_" + (args.model) + "_" + str(args.likessuffix) +"_" + str(int(args.regionlen)) + "_" + str(int(args.thresshold)) + "_" + str(args.cutoff) + ".txt"
		print(regionfilename)
		if os.path.isfile(regionfilename):
			regionfiles.append(regionfilename)
			takepops.append(pop)
	if len(regionfiles) == 0:
		return
	else:
		totalselregions = 0
		print('loaded regions from ' + str(len(regionfiles)) + " files...")
		boldstyle = easyxf('font: bold 1;')
		wrapstyle = easyxf('alignment: wrap on, vert center, horiz center')
		book = Workbook()
		sheet1 = book.add_sheet('gw significant regions')
		header = ['chrom', 'start', 'end', 'len (kb)', 'genes', 'pop']
		for icol in range(len(header)):
			sheet1.write(0, icol, header[icol], boldstyle)

		#get inclusive overlap
		irow = 0
		for iregionfilename in range(len(regionfiles)):
			regionfilename = regionfiles[iregionfilename]
			pop = takepops[iregionfilename]
			#for regionfilename in regionfiles:
			genes = BedTool(args.gene_bedfile) 
			regions = BedTool(regionfilename)
			intersect = regions.intersect(genes,  wa = True, wb = True)

			#narrow down
			geneDict = {}
			for item in intersect:
				selregion_chr = item[0]
				selregion_start, selregion_end = item[1], item[2]
				key = (selregion_chr, selregion_start, selregion_end)
				generegion_id = item[6]
				if key not in geneDict.keys():
					geneDict[key] = [generegion_id]
				else:
					geneDict[key].append(generegion_id)

			#for irow in range(len(regions)):
			for region in regions:
				totalselregions +=1
				irow +=1
				#region = regions[irow]
				chrom, start, end = region[0], region[1], region[2]
				key = (chrom, start, end)
				regionlen = int(end) - int(start)
				kb_regionlen=round(regionlen/1000)

				if key in geneDict.keys():
					genelist = geneDict[key]
					genes = set(genelist)
					genestring = ""
					for gene in genes:
						genestring += gene + ", "
					genestring = genestring[:-2]
					sheet1.write(irow+1, 0, chrom, wrapstyle)
					sheet1.write(irow+1, 1, int(start), wrapstyle)
					sheet1.write(irow+1, 2, int(end), wrapstyle)
					sheet1.write(irow+1, 3, kb_regionlen, wrapstyle)
					sheet1.write(irow+1, 4, genestring, wrapstyle)
					sheet1.write(irow+1, 5, pop, wrapstyle)
				else:
					sheet1.write(irow+1, 0, chrom, wrapstyle)
					sheet1.write(irow+1, 1, int(start), wrapstyle)
					sheet1.write(irow+1, 2, int(end), wrapstyle)
					sheet1.write(irow+1, 3, kb_regionlen,wrapstyle)			
					sheet1.write(irow+1, 5, pop,wrapstyle)

	colWidths = [10, 15, 15, 10, 25, 10]
	for icol in range(len(colWidths)):
		sheet1.col(icol).width = colWidths[icol] * 256 #~x char wide	

	book.save(args.savefile)
	book.save(TemporaryFile())
	print('wrote ' + str(totalselregions) + ' significant regions to: ' + args.savefile)
	return
def execute_manhattan(args):
	selpop = args.emppop
	model = args.model
	likessuffix = args.likessuffix
	savename = args.savefilename
	#nRep = args.nrep
	###############################
	### LOAD NEUTRAL SIM VALUES ###
	###############################
	if likessuffix == "neut":
		vsNeut = True
	elif likessuffix == "linked":
		vsNeut = False
	modelpops = {'YRI':1, 'CEU':2, 'CHB':3, 'BEB':4}
	pop = modelpops[selpop]

	#EXPERIMENTAL
	#This is not how Shari actually normalized for figure 1 2013.
	#if args.poolModels:
	#	all_neut_rep_scores = []
	#	for model in ['default_112115_825am', 'nulldefault', 'nulldefault_constantsize', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm']:
	#		neut_rep_scores = load_normed_simscores(model, pop, vsNeut=vsNeut, numRep=nRep)
	#		all_neut_rep_scores.extend(neut_rep_scores)
	#	neut_rep_scores = all_neut_rep_scores
	#else:
	#	neut_rep_scores = load_normed_simscores(model, pop, vsNeut=vsNeut, numRep=nRep)
	#print('loaded ' + str(len(neut_rep_scores)) + ' neut simscores...') 
	neut_rep_scores = []

	#############################
	### LOAD EMPIRICAL VALUES ###
	#############################
	all_emp_pos, all_emp_scores = [], []
	nSnps = 0
	for chrom in range(1,23):
		thesepos, thesescores = [], []
		emp_cms_filename = get_emp_cms_file(selpop, model, likessuffix, chrom, normed=True)
		print('loading chr ' + str(chrom) + ": " + emp_cms_filename)
		if not os.path.isfile(emp_cms_filename):
			print("missing: " + emp_cms_filename)
			break
		physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(emp_cms_filename)
		all_emp_pos.append(physpos)
		all_emp_scores.append(cms_normed)

	###########################
	### DRAW MANHATTAN PLOT ###
	###########################
	f, ax = plt.subplots(1)
	if args.zscores:
		calc_zscores = True
	else:
		calc_zscores = False

	plotManhattan(ax, neut_rep_scores, all_emp_scores, all_emp_pos, nSnps, zscores=calc_zscores, maxSkipVal = args.maxSkipVal)

	plt.savefig(savename)
	print('saved to: ' + savename)
	return
def execute_extended_manhattan(args):
	selpop = args.emppop
	model = args.model
	likessuffix = args.likessuffix
	savename = args.savefilename
	numChr = 22
	if likessuffix == "neut":
		vsNeut = True
	elif likessuffix == "linked":
		vsNeut = False
	modelpops = {'YRI':1, 'CEU':2, 'CHB':3, 'BEB':4}
	pop = modelpops[selpop]

	f, axarr = plt.subplots(numChr, 1, sharex = True)

	plt.xlabel('position')
	plt.ylabel('cms_gw normed score')
	#plt.tick_params(axis='y', left='off') only one subplot?

	all_emp_pos, all_emp_scores = [], []
	for chrom in range(1,numChr +1):
		emp_cms_filename = get_emp_cms_file(selpop, model, likessuffix, chrom, normed=True)
		print('loading chr ' + str(chrom) + ": " + emp_cms_filename)
		if not os.path.isfile(emp_cms_filename):
			print("missing: " + emp_cms_filename)
			break
		physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(emp_cms_filename)

		iax = chrom-1
		ax = axarr[iax]
		plotManhattan_extended(ax, cms_normed, physpos, chrom)
		all_emp_pos.append(physpos)
		all_emp_scores.append(cms_normed)
	################################
	## HILITE SIGNIFICANT REGIONS ##
	################################

	if args.regionsfile is not None:
		regionchrs, regionstarts, regionends = loadregions(args.regionsfile)

		print('loaded ' + str(len(regionchrs)) + ' significant regions from ' + args.regionsfile)
		for iregion in range(len(regionchrs)):
			regionchr, regionstart, regionend = regionchrs[iregion], regionstarts[iregion], regionends[iregion]
			this_chrom = int(regionchr.strip('chr'))
			ichrom = this_chrom-1
			chrompos, chromscores = all_emp_pos[ichrom], all_emp_scores[ichrom]
			zipped = zip(chrompos, chromscores)
			plotpos, plotvals = [], []
			for locus in zipped:
				if locus[0] >= regionstart:
					plotpos.append(locus[0])
					plotvals.append(locus[1])
				if locus[0] > regionend:
					break

			axarr[ichrom].plot(plotpos, plotvals, color="red")
			#axarr[ichrom].plot([regionstart, regionend], [0, 0], color="red")
	#axarr.get_yaxis().set_visible(False)
	plt.savefig(savename)
	print('saved to: ' + savename)
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_power()
	args = runparser.parse_args()

	# if called with no arguments, print help
	if len(sys.argv)==1:
		runparser.parse_args(['--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline

