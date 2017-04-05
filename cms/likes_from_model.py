#!/usr/bin/env python
## top-level script for generating probability distributions for component scores as part of CMS 2.0. 
## last updated: 04.04.2017 	vitti@broadinstitute.org

import matplotlib as mp 
mp.use('agg') 
from dists.freqbins_func import run_traj, get_bin_strings, get_bins, check_create_dir, check_create_file, write_bin_paramfile, execute, get_concat_files, get_info_from_tped_name
from dists.scores_func import calc_ihs, calc_delihh, calc_xpehh, calc_fst_deldaf, read_neut_normfile, norm_neut_ihs, norm_sel_ihs, norm_neut_xpehh, norm_sel_xpehh, get_sim_compscore_files, get_scores_from_files, get_compscores_from_files_flatten#get_compscores_from_files
from dists.likes_func import plot_pdf_comparison_from_scores, get_plot_pdf_params, quick_load_likes, quick_cl_from_likes, quick_clr_from_likes
import numpy as np
import argparse
import sys, os, subprocess
import matplotlib.pyplot as plt

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")
	generate_sel_bins_parser = subparsers.add_parser('generate_sel_bins', help="Pre-processing step: generate directories with parameter files divided by selDAF, according to specified bins") #nix or clean
	generate_sel_bins_parser.add_argument('outputDir', type=str, action='store', help='location to write cosi output')
	
	#####################################
	## SIMULATE SELECTION TRAJECTORIES ##
	#####################################
	get_sel_traj_parser = subparsers.add_parser('get_sel_traj', help='Run forward simulations of selection trajectories to populate selscenarios by final allele frequency before running coalescent simulations for entire sample.')
	get_sel_traj_parser.add_argument('traj_outputname', type=str, action='store', help="write to file")	
	get_sel_traj_parser.add_argument('--maxAttempts', type=int, action='store', help='maximum number of attempts to generate a trajectory before re-sampling selection coefficient / start time', default=100)
	
	##########################
	## RUN FULL SIMULATIONS ##
	##########################
	run_neut_sim_parser = subparsers.add_parser('run_neut_sim', help='run neutral simulations')
	run_sel_sim_parser = subparsers.add_parser('run_sel_sim', help='run sims with selection')
	run_sel_sim_parser.add_argument('traj_infilename', type=str, action='store', help="selection trajectory")

	####################################
	## GENERATE SCORES FROM SIMULATES ##
	####################################
	run_repscores_parser = subparsers.add_parser('run_repscores', help="run composite score calculations for simulated tped")
	run_repscores_parser.add_argument('score_basedir', type=str, action='store', help="parent directory in which to generate/populate folders for each composite score")
	run_repscores_parser.add_argument('inputTpedFile', type=str, action='store', help="TPED file with simulate data for putative selected population") 
	run_repscores_parser.add_argument('--simRecomFile', type=str, action="store", help="location of input recom file", default="/n/home08/jvitti/params/test_recom.recom")
	get_neut_norm_params_parser = subparsers.add_parser('get_neut_norm_params', help="jointly normalize scores from neutral replicates to get parameters with which to normalize all replicates")	
	norm_from_neut_params_parser = subparsers.add_parser('norm_from_neut_params', help="normalize component scores according to neutral distribution")
	norm_from_neut_params_parser.add_argument('--selbin', action='store', help="e.g. 0.10 -- if excluded, normalize neutral simulates")

	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################
	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='Collate scores from simulated data in order to generate component test probability distributions.')
	likes_from_scores_parser.add_argument('--nrep_neut', type=int, action="store", help="number of neutral replicates", default=1000)
	likes_from_scores_parser.add_argument('--nrep_sel', type=int, action="store", help="number of replicates per selection scenario bin", default=500)
	likes_from_scores_parser.add_argument('--binMerge', action="store_true", help="if included, generate higher-order groupings of selscenario frequency bins")
	likes_from_scores_parser.add_argument('--foldDist', action="store_true", help="if included, take absolute value of score values (fold distribution over y-axis)", default=False)
	likes_from_scores_parser.add_argument('--save_suffix', type=str, action="store", help="optional string to add to filenames (e.g. to avoid overwrite)", default=None)
	likes_from_scores_parser.add_argument('--save_dir', type=str, action="store", help="if included, specify alternate outpot location", default=None)
	likes_from_scores_parser.add_argument('--save_likelihoods', action="store_true", help="in addition to plotting, save data", default=False)
	if True:
		likes_from_scores_parser.add_argument('--ihs', action="store_true", help="visualize likelihoods for iHS", default=False)
		likes_from_scores_parser.add_argument('--delihh', action="store_true", help="visualize likelihoods for delIHH", default=False)
		likes_from_scores_parser.add_argument('--nsl', action="store_true", help="visualize likelihoods for nSL", default=False)
		likes_from_scores_parser.add_argument('--xpehh', action="store_true", help="visualize likelihoods for XP-EHH", default=False)
		likes_from_scores_parser.add_argument('--deldaf', action="store_true", help="visualize likelihoods for delDAF", default=False)
		likes_from_scores_parser.add_argument('--fst', action="store_true", help="visualize likelihoods for Fst", default=False)
	plot_likes_vs_scores_parser = subparsers.add_parser('plot_likes_vs_scores', help='useful QC step; visualize the function that converts score values into likelihood contributions towards composite scores')
	plot_likes_vs_scores_parser.add_argument('inputLikesPrefix', action="store", help="filename minus causal/linked/neutral.txt")
	plot_likes_vs_scores_parser.add_argument('--default_nregion_snps', type=int, action="store", help="presumptive number of SNPs in region (determines prior distributions for within-region CMS calculations)", default=5000)

	#################
	## SHARED ARGS ## 
	################# 
	for common_parser in [run_neut_sim_parser, run_sel_sim_parser, run_repscores_parser, get_neut_norm_params_parser, norm_from_neut_params_parser]:
		common_parser.add_argument('--checkOverwrite', action="store_true", default=True)
	for run_sim_parser in [run_neut_sim_parser, run_sel_sim_parser]:
		run_sim_parser.add_argument('writeBase', type=str, action='store', help="write prefix")
		run_sim_parser.add_argument('--dropSings', type=float, action='store',  help='randomly thin global singletons from output dataset to model ascertainment bias', default=.25)	
	for cosi_parser in [get_sel_traj_parser, run_neut_sim_parser, run_sel_sim_parser]:
		cosi_parser.add_argument('--cosiBuild', type=str, action='store', help='which version of cosi to run', default="coalescent")
		cosi_parser.add_argument('inputParamFile', type=str, action='store', help='file with model specifications for input')
		#cosi_parser.add_argument('--genmapRandomRegions', action='store_true', help='cosi option to sub-sample genetic map randomly from input')	
	for cms_preconda_parser in [run_repscores_parser, get_neut_norm_params_parser, norm_from_neut_params_parser]: #TEMPORARY; conda will obviate
		cms_preconda_parser.add_argument('--cmsdir', type=str, action='store', help="location of CMS scripts (TEMPORARY; conda will obviate)", default="/n/home08/jvitti/cms/cms/") 
	for interior_score_parser in [get_neut_norm_params_parser, likes_from_scores_parser]:
		interior_score_parser.add_argument('--edge', type=int, action="store", help="use interior of replicates; define per-end bp. (e.g. 1.5Mb -> 1Mb: 250000)", default=250000)
		interior_score_parser.add_argument('--chromlen', type=int, action="store", help="per bp (1.5mb = 1500000)", default=1500000)
	for norm_parser in [get_neut_norm_params_parser, norm_from_neut_params_parser]:
		norm_parser.add_argument('--score', type=str, action='store', default='ihs')
		norm_parser.add_argument('--simpop', type=int, action='store', default=1)
		norm_parser.add_argument('--altpop', type=int, action='store', default=2)
		norm_parser.add_argument('--nrep', type=int, action='store', default=100)
	for norm_sims_parser in [get_neut_norm_params_parser, norm_from_neut_params_parser, likes_from_scores_parser]:
		norm_sims_parser.add_argument('modeldir', type=str, action="store", help="location of component score folders for demographic scenario")
	for selbin_parser in [generate_sel_bins_parser, get_sel_traj_parser, likes_from_scores_parser]:
		selbin_parser.add_argument('--freqRange', type=str, help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
		selbin_parser.add_argument('--nBins', type=int, help="number of frequency bins", default=9)		
	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
### Run simuates from specified demographic
### model under various scenarios
def execute_generate_sel_bins(args):
	''' pre-processing step for sel_trajs(->sel_sim) ''' #hmmm nix this? or else, fix write_bin_paramfile: merge with quick_write_inclusive_sweep_paramfiles.py
	freqRange = args.freqRange
	nBins = args.nBins
	runDir = args.outputDir
	neutParamfile = args.inputParamFile
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(freqRange, nBins)
	for ibin in range(nBins):
		populateDir = runDir + "sel_" + bin_medians_str[ibin]
		binDir = check_create_dir(populateDir)
		bounds = bin_starts[ibin], bin_ends[ibin]
		paramfilename = populateDir + "/params"
		write_bin_paramfile(neutParamfile, paramfilename, bounds)
		print('wrote to: ' + paramfilename)	
	return
def execute_get_sel_traj(args):
	'''generate forward trajectories of simulated allele frequencies for demographic scenarios with selection'''
	traj_outputname = args.traj_outputname
	cosibuild = args.cosiBuild
	paramfilename = args.inputParamFile
	maxattempts = args.maxAttempts
	run_traj(traj_outputname, cosibuild, paramfilename, maxattempts)
	return
def execute_run_neut_sim(args):
	''' generates tped data for one neutral replicate from a demographic model parameter file using coalescent simulator cosi '''
	cosibuild = args.cosiBuild
	outbase = args.writeBase
	paramfilename = args.inputParamFile
	dropSing = args.dropSings
	cmd = cosibuild
	argstring = "-p " + paramfilename + " --genmapRandomRegions --drop-singletons " + str(dropsSing) + " --tped " + outbase + " --output-gen-map"
	#print(cmd + " " + argstring)	
	cosicreatefilename = outbase + "_0_1.tped"
	cosi_movedfilename = outbase + "_1.tped"	
	proceed = check_create_file(cosicreatefilename, args.checkOverwrite)
	if proceed:
	#if not os.path.isfile(cosicreatefilename):
		fullCmd = cmd + " " + argstring
		print(fullCmd)
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
def execute_run_sel_sim(args):
	''' generates tped data for one replicate from a demographic model parameter file and selection trajectory file using coalescent simulator cosi '''
	trajectory = args.traj_infilename
	cosibuild = args.cosiBuild
	outbase = args.writeBase
	paramfilename = args.inputParamFile
	dropSing = args.dropSings	
	cmd = "env COSI_NEWSIM=1 env COSI_LOAD_TRAJ=" + trajectory + " " + cosibuild 
	argstring = "-p " + paramfilename + " --genmapRandomRegions --drop-singletons " + str(dropsSing) + " --tped " + outbase + " --output-gen-map"
	#print(cmd + " " + argstring)
	cosicreatefilename = outbase + "_0_1.tped"
	cosi_movedfilename = outbase + "_1.tped"
	proceed = check_create_file(cosi_movedfilename, args.checkOverwrite)
	if proceed:
		fullCmd = cmd + " " + argstring
		print(fullCmd)
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

### Calculate component scores from 
### simulated data
def execute_run_repscores(args):
	''' for one simulated replicate (agnostic wrt neut/sel), generate all component scores '''
	basedir = args.score_basedir
	cmsdir = args.cmsdir 
	tped = args.inputTpedFile
	repNum, pop, tpeddir = get_info_from_tped_name(tped)
	checkOverwrite = args.checkOverwrite #nix/fix
	simRecomFile = args.simRecomFile
	assert os.path.isfile(tped)
	for scorefiledir in ['ihs', 'delihh', 'nsl', 'xpehh', 'fst_deldaf']:
		#dircmd = "mkdir -p " + basedir + scorefiledir + "/"
		#subprocess.check_output( dircmd.split() )
		check_create_dir(basedir + scorefiledir)

	####### Calculate per-population
	####### scores: iHS, delIHH, nSL
	ihs_commandstring = "python " + cmsdir + "scans.py selscan_ihs"
	ihs_outfileprefix = basedir + "ihs/rep" + str(repNum) + "_" + str(pop) 
	ihs_unnormedfile = ihs_outfileprefix + ".ihs.out"
	ihs_argstring = tped + " " + ihs_outfileprefix + " --threads 7 "
	ihs_fullcmd = ihs_commandstring + " " + ihs_argstring
	ihs_normedfile = ihs_unnormedfile + ".norm"
	proceed = check_create_file(ihs_normedfile, args.checkOverwrite)
	if proceed:
		print(ihs_fullcmd)
		execute(ihs_fullcmd)
	delihh_commandstring = "python " + cmsdir + "composite.py delihh_from_ihs"
	delihh_unnormedfile =  basedir + "delihh/rep" + str(repNum) + "_" + str(pop) + ".txt"
	delihh_argstring = ihs_unnormedfile + " "+ delihh_unnormedfile
	delihh_fullcmd = delihh_commandstring + " " + delihh_argstring 
	delihh_normedfile = delihh_unnormedfile + ".norm"
	proceed = check_create_file(delihh_normedfile, args.checkOverwrite)
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

		fstdeldaf_commandstring = "python " + cmsdir + "composite.py freqscores"
		fstdeldaf_outfilename = basedir + "fst_deldaf/rep" + str(repNum) + "_" + str(pop) + "_" + str(altpop)
		fstdeldaf_argumentstring = tped + " " + tped2 + " " + simRecomFile + " " + fstdeldaf_outfilename 
		fstdeldaf_fullcmd = fstdeldaf_commandstring + " " + fstdeldaf_argumentstring 
		proceed = check_create_file(fstdeldaf_outfilename, args.checkOverwrite)
		if proceed:
			print(fstdeldaf_fullcmd)
			execute(fstdeldaf_fullcmd)
	return
def execute_get_neut_norm_params(args):
	''' creates a concatenated file and saves bin output if run for the first time'''
	pop = args.simpop
	numReps = args.nrep
	basedir = args.modeldir 
	if basedir[-1] != "/":
		basedir += "/"
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop
	chrlen, edge = args.chromlen, args.edge
	startbound, endbound = int(edge), chrlen - int(edge) #define replicate interior (conservative strategy to avoid edge effects)

	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + "neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		altpop = args.altpop
		concatfilebase = basedir + "neut/concat_" + str(pop) + "_" + str(altpop) + "_"
	else:
		print('must call per composite score')
		sys.exit(0)

	concatfilename = concatfilebase + score + ".txt"
	binfilename = concatfilebase + score + ".bins"

	if not os.path.isfile(binfilename):
		if not os.path.isfile(concatfilename):
			repfiles = []
			for irep in range(1, numReps+1):
				if score == 'ihs':
					unnormedfile = basedir + "neut/ihs/rep" + str(irep) + "_" + str(pop) + ".ihs.out"
					physpos_ind = 1
				elif score == "delihh":
					unnormedfile = basedir + "neut/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
					physpos_ind = 1
				elif score == "nsl":
					unnormedfile = basedir + "neut/nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
					physpos_ind = 1
				elif score == "xpehh":
					unnormedfile = basedir + "neut/xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
					physpos_ind = 1
				elif score == "fst":
					unnormedfile = basedir + "neut/fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
					physpos_ind = 0
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
				for line in readfile:
					entries = line.split()
					physpos = int(entries[physpos_ind])
					if physpos >= startbound and physpos <= endbound:
						concatfile.write(line)
				readfile.close()
			concatfile.close()
			print('wrote to: ' + concatfilename)

		#already have concatfile
		infilename = concatfilename
		outfilename = concatfilename + ".norm" #or just write to /tmp ?
		argstring = infilename #+ " " + outfilename
		if score in ['ihs', 'delihh']:
			cmd = "python " + cmsdir + "scans.py selscan_norm_ihs"
		elif score in ['nsl']:
			cmd = "python " + cmsdir + "scans.py selscan_norm_nsl"
		elif score in ['xpehh']:
			cmd = "python " + cmsdir + "scans.py selscan_norm_xpehh"
		else:
			cmd = ""
		fullcmd = cmd + " " + argstring
	
		alreadyExists = False
		if args.checkOverwrite:
			if not os.path.isfile(binfilename): #check for overwrite
				alreadyExists = False
			else:
				alreadyExists = True				
		if alreadyExists == False:
			print(fullcmd)
			#execute(fullcmd)
			with open(binfilename, 'w') as outfile:
				subprocess.check_output( fullcmd.split(), stderr=outfile)
			outfile.close()
			print('wrote to: ' + binfilename)
	return
def execute_norm_from_neut_params(args): 
	''' using parameters from neutral simulates, normalizes component scores '''
	pop = args.simpop
	numReps = args.nrep
	basedir = args.modeldir
	if basedir[-1] != "/":
		basedir += "/"
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop
	if args.selbin is None:
		scenario_dir = "neut/"
	else:
		scenario_dir = "sel" + str(pop) + "/sel_" + str(args.selbin) + "/"
	concatfilename, binfilename = get_concat_files(pop, score, altpop, basedir=basedir) #I NEED TO MAKE THIS HANDLE SEL SITUATION
	print('loading normalization parameters from ' + binfilename + "...")
	###############
	## NORMALIZE ##
	###############
	for irep in range(1, numReps+1):
		if score in ['ihs']:
			unnormedfile = basedir + scenario_dir + "ihs/rep" + str(irep) + "_" + str(pop)  + ".ihs.out"
			norm_sel_ihs(unnormedfile, binfilename)
		elif score in ['delihh']:
			unnormedfile = basedir + scenario_dir + "delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
			norm_sel_ihs(unnormedfile, binfilename)
		elif score in ['nsl']:
			unnormedfile = basedir + scenario_dir + "nsl/rep" + str(irep) + "_" + str(pop)+ ".nsl.out"
			norm_sel_ihs(unnormedfile, binfilename)
		elif score in ['xpehh']:
			unnormedfile = basedir + scenario_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
			norm_sel_xpehh(unnormedfile, binfilename)
		else: #if score in ['fst']:
			print('currently handling this manually: rewrite_fst_bins.py')
			pass
		if irep%100 == 0:
			print("currently rep: " + str(irep))
	return	

### Define component score likelihood 
### functions as histograms of simulated scores
def execute_likes_from_scores(args): 
	''' define likelihood function for component score based on histograms of score values from simulated data.
	each run is per-score, per-scenario. does all pops / poppairs  '''
	modeldir, pops = args.modeldir, [1, 2, 3, 4]
	if modeldir[-1] != "/":
		modeldir += "/"
	modeldir_entries = modeldir.split('/')
	model = modeldir_entries[-2]
	nPerBin_Neut, nPerBin_Sel = args.nrep_neut, args.nrep_sel
	freqRange, nBins = args.freqRange, args.nBins #define selscenario bins
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(freqRange, nBins)
	chrlen, edge = args.chromlen, args.edge #define replicate interior (conservative strategy to avoid edge effects)
	startbound, endbound = int(edge), chrlen - int(edge) 
	#print(str([args.ihs, args.delihh, args.nsl, args.xpehh, args.deldaf, args.fst]))
	if [args.ihs, args.delihh, args.nsl, args.xpehh, args.deldaf, args.fst].count(True) != 1:
		print('must call one score at a time')
		sys.exit(0)
	#####################
	## LOAD NEUT FILES ##
	#####################
	neutdir = modeldir + "neut/"
	all_completed_neut = []
	for pop in pops:
		completed_neut = []
		for irep in range(1, nPerBin_Neut+1):
			loaded_neut_files = get_sim_compscore_files(pop, irep, neutdir) 
			all_present = [os.path.isfile(item) for item in loaded_neut_files]
			if (all_present.count(True)) == 9: #replicate done
				completed_neut.append(loaded_neut_files)
		all_completed_neut.append(completed_neut)
	print("loaded " + str(sum([len(item) for item in all_completed_neut])) + " neutral replicates with complete component scores from " + neutdir)
	####################
	## LOAD SEL FILES ## 
	####################
	binlabels = bin_medians_str 
	all_completed_sel, selcounter = [], 0
	for pop in pops:
		seldir = modeldir + "sel" + str(pop) + "/"
		completed_sel = []
		## generate likelihood distributions for selection simulates 
		## from multiple frequency bins
		if args.binMerge:  
			#need: input bins and name of cluster
			selbin_clusters = [["0.10", "0.20", "0.30"], ["0.40", "0.50", "0.60",], ["0.70", "0.80", "0.90",],["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]] 
			cluster_names = ["low", "mid", "hi","allfreq"]
			for icluster in range(len(selbin_clusters)):
				cluster, clustername = selbin_clusters[icluster], cluster_names[icluster]
				completed_cluster = []
				for selbin in cluster:
					bindir = seldir + "sel_" + str(selbin) + "/"				
					for irep in range(1, nPerBin_Sel+1):
						loaded_sel_files = get_sim_compscore_files(pop, irep, bindir)
						all_present = [os.path.isfile(item) for item in loaded_sel_files]
						if (all_present.count(True)) == 9: #replicate done
							completed_cluster.append(loaded_sel_files)
							selcounter +=1
				completed_sel.append(completed_cluster)
			all_completed_sel.append(completed_sel) #[ipop][iCLUSTER][irep][iscore]
			nChunks, chunk_labels = len(cluster_names), cluster_names
		
		## fine-grained (ie per-sel-bin) generation of 
		## score probability density functions
		else:
			for selbin in binlabels:
				bindir = seldir + "sel_" + str(selbin) + "/"
				completed_bin = []
				for irep in range(1, nPerBin_Sel+1):
					loaded_sel_files = get_sim_compscore_files(pop, irep, bindir)
					all_present = [os.path.isfile(item) for item in loaded_sel_files]
					if (all_present.count(True)) == 9: #replicate done
						completed_bin.append(loaded_sel_files)
						selcounter +=1
				completed_sel.append(completed_bin)
			all_completed_sel.append(completed_sel) #[ipop][ibin][irep][iscore]
			nChunks, chunk_labels = len(binlabels), bin_medians_str
	print("loaded " + str(selcounter) + " selection replicates with complete component scores")
	#################################
	## SORT SCORES INTO HISTOGRAMS ## 
	################################# 
	for ichunk in range(nChunks): #iterate over selfreqs
		chunkstring = chunk_labels[ichunk]
		####################
		## PER POP SCORES ##
		####################
		if (args.ihs or args.delihh or args.nsl): #per pop
			###########
			### SORT ## 
			###########
			if args.ihs:
				score = "ihs"
				print('binning iHS scores...')
				values = get_scores_from_files(all_completed_neut, all_completed_sel, 0, ichunk, startbound, endbound, foldDists = args.foldDist)
			if args.delihh:
				score = "delihh"
				print('binning delIHH scores...')
				values = get_scores_from_files(all_completed_neut, all_completed_sel, 1, ichunk, startbound, endbound, foldDists = args.foldDist)	
			if args.nsl:
				score = "nsl"
				print('binning nSL scores...')
				values = get_scores_from_files(all_completed_neut, all_completed_sel, 2, ichunk, startbound, endbound, foldDists = args.foldDist)
			if True:
				pop1vals, pop2vals, pop3vals, pop4vals = values
				neut_1, causal_1, linked_1 = pop1vals
				neut_2, causal_2, linked_2 = pop2vals
				neut_3, causal_3, linked_3 = pop3vals
				neut_4, causal_4, linked_4 = pop4vals
			##########
			## PLOT ##
			##########
			minVal, maxVal, nProbBins, annotate = get_plot_pdf_params(score)
			plot_title = "PDF for " + score + ", sel_" + str(chunkstring)
			output_dir = ""
			if args.save_dir is not None:
				output_dir = args.save_dir + model + "_"
			else:
				output_dir = modeldir + "likes/"

			savefilebase = output_dir + score + "_sel_" + chunkstring 
			if args.save_suffix is not None:
				savefilebase += "_" + args.save_suffix

			savefilename = savefilebase + ".png"
			#savefilebase = modeldir + "likes/" + score + "_sel_" + chunkstring
			likes_savebase_1 = output_dir + score + "_sel1_" + chunkstring + "_" 
			likes_savebase_2 = output_dir + score + "_sel2_" + chunkstring + "_" #+ args.save_suffix + "_"
			likes_savebase_3 = output_dir + score + "_sel3_" + chunkstring + "_" #+ args.save_suffix + "_"
			likes_savebase_4 = output_dir + score + "_sel4_" + chunkstring + "_" #+ args.save_suffix + "_"
			if args.save_suffix is not None:
				likes_savebase_1 += args.save_suffix + "_"
				likes_savebase_2 += args.save_suffix + "_"
				likes_savebase_3 += args.save_suffix + "_"
				likes_savebase_4 += args.save_suffix + "_"

			f1, (ax1, ax2, ax3, ax4)  = plt.subplots(4, sharex=True, sharey=True)
			plt.suptitle(plot_title)
			plot_pdf_comparison_from_scores(ax1, neut_1, causal_1, linked_1, minVal, maxVal, nProbBins, "1 (AFR)", likes_savebase_1, saveFiles = args.save_likelihoods, annotate = annotate)
			plot_pdf_comparison_from_scores(ax2, neut_2, causal_2, linked_2, minVal, maxVal, nProbBins, "2 (EUR)", likes_savebase_2, saveFiles = args.save_likelihoods, annotate = annotate)
			plot_pdf_comparison_from_scores(ax3, neut_3, causal_3, linked_3, minVal, maxVal, nProbBins, "3 (EAS)", likes_savebase_3, annotate = annotate, saveFiles = args.save_likelihoods)
			plot_pdf_comparison_from_scores(ax4, neut_4, causal_4, linked_4, minVal, maxVal, nProbBins, "4 (SAS)", likes_savebase_4, annotate = annotate, saveFiles = args.save_likelihoods)
			ax4.set_xlabel('score value')
			#f1.subplots_adjust(hspace=0)
			#f1.ylabel("p(score)")
			plt.savefig(savefilename)
			plt.close()
			print('saved to ' + savefilename)

		##################
		## PER POP COMP ##
		##################
		if (args.xpehh or args.deldaf or args.fst): #per pop-comp
			###########
			### SORT ##
			###########
			if args.xpehh:
				score = "xpehh"
				print('binning XP-EHH scores...')
				values = get_compscores_from_files_flatten(all_completed_neut, all_completed_sel, "xpehh", ichunk, startbound, endbound, foldDists = args.foldDist)
			if args.deldaf:
				score = "deldaf"
				print('binning delDAF scores...')
				values = get_compscores_from_files_flatten(all_completed_neut, all_completed_sel, "deldaf", ichunk, startbound, endbound, foldDists = args.foldDist)	
			if args.fst:
				score = "fst"
				print('binning Fst scores...')
				values = get_compscores_from_files_flatten(all_completed_neut, all_completed_sel, "fst", ichunk, startbound, endbound, foldDists = args.foldDist)
			if True:
				pop1vals, pop2vals, pop3vals, pop4vals = values
				neutvals1, causalvals1, linkedvals1 = pop1vals
				neutvals2, causalvals2, linkedvals2 = pop2vals
				neutvals3, causalvals3, linkedvals3 = pop3vals
				neutvals4, causalvals4, linkedvals4 = pop4vals

			##########
			## PLOT ##
			##########
			minVal, maxVal, nProbBins, annotate = get_plot_pdf_params(score)
			output_dir = ""
			if args.save_dir is not None:
				output_dir = args.save_dir + model + "_"
			else:
				output_dir = modeldir + "likes/"

			savefilebase = output_dir + score + "_sel_" + str(chunkstring)
			if args.save_suffix is not None:
				savefilebase += "_" + args.save_suffix
			savefilename = savefilebase + ".png"
	
			plot_title = "PDF for " + score + ", " + str(chunkstring)

			likes_savebase_1 = output_dir + score + "_sel1_" + chunkstring + "_" 
			likes_savebase_2 = output_dir + score + "_sel2_" + chunkstring + "_" #+ args.save_suffix + "_"
			likes_savebase_3 = output_dir + score + "_sel3_" + chunkstring + "_" #+ args.save_suffix + "_"
			likes_savebase_4 = output_dir + score + "_sel4_" + chunkstring + "_" #+ args.save_suffix + "_"
			if args.save_suffix is not None:
				likes_savebase_1 += args.save_suffix + "_"
				likes_savebase_2 += args.save_suffix + "_"
				likes_savebase_3 += args.save_suffix + "_"
				likes_savebase_4 += args.save_suffix + "_"
			f1, (ax1, ax2, ax3, ax4)  = plt.subplots(4, sharex=True, sharey=True)
			plt.suptitle(plot_title)
			plot_pdf_comparison_from_scores(ax1, neutvals1, causalvals1, linkedvals1, minVal, maxVal, nProbBins, "1 (AFR)", likes_savebase_1, saveFiles = args.save_likelihoods, annotate = annotate)
			plot_pdf_comparison_from_scores(ax2, neutvals2, causalvals2, linkedvals2, minVal, maxVal, nProbBins, "2 (EUR)", likes_savebase_2, saveFiles = args.save_likelihoods, annotate = annotate)
			plot_pdf_comparison_from_scores(ax3, neutvals3, causalvals3, linkedvals3, minVal, maxVal, nProbBins, "3 (EAS)", likes_savebase_3, annotate = annotate, saveFiles = args.save_likelihoods)
			plot_pdf_comparison_from_scores(ax4, neutvals4, causalvals4, linkedvals4, minVal, maxVal, nProbBins, "4 (SAS)", likes_savebase_4, annotate = annotate, saveFiles = args.save_likelihoods)
			ax4.set_xlabel('score value')
			#plot_pdf_comparison_from_scores(ax1, neuta, causala, linkeda, minVal, maxVal, nProbBins, ax_ylabela, savefilename_a, annotate = annotate, saveFiles = args.save_likelihoods)
			#ax3.set_xlabel('score value')
			#f1.subplots_adjust(hspace=0)
			plt.savefig(savefilename)
			plt.close()
			print('saved to ' + savefilename)
	return
def execute_plot_likes_vs_scores(args):
	''' useful QC step; visualize the function that converts score values into likelihood contributions towards composite scores '''
	inputprefix = args.inputLikesPrefix
	nregionsnps = args.default_nregion_snps
	causal_likes_filename, linked_likes_filename, neut_likes_filename = inputprefix + "_causal.txt", inputprefix + "_linked.txt", inputprefix + "_neut.txt"

	causal_likes = quick_load_likes(causal_likes_filename)
	linked_likes = quick_load_likes(linked_likes_filename)
	neut_likes = quick_load_likes(neut_likes_filename)

	comp_likes = quick_cl_from_likes(causal_likes, linked_likes, nSnp = nregionsnps, takeLn = True)
	comp_like_ratios = quick_clr_from_likes(causal_likes, neut_likes, takeLn = True)

	savefilename = inputprefix + "_clr_likes_vs_likesscores.png"
	fig, ax = plt.subplots()
	ax.scatter(np.arange(len(comp_like_ratios)), comp_like_ratios)
	plt.savefig(savefilename)
	print('saved to: ' + savefilename)
	plt.close()

	savefilename = inputprefix + "_cl_likes_vs_likesscores.png"
	fig, ax = plt.subplots()
	ax.scatter(np.arange(len(comp_likes)), comp_likes)
	plt.savefig(savefilename)
	print('saved to: ' + savefilename)
	plt.close()
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_likes_from_model()
	args = runparser.parse_args()

	# if called with no arguments, print help
	if len(sys.argv)==1:
		runparser.parse_args(['--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline
