#!/usr/bin/env python
## top-level script for generating probability distributions for component scores as part of CMS 2.0. 
## last updated: 02.27.2017 	vitti@broadinstitute.org #must update docstrings

import matplotlib as mp 
mp.use('agg') 
from power.power_func import write_master_likesfile 
from dists.likes_func import get_old_likes, read_likes_file, plot_likes, get_hist_bins, read_demographics_from_filename, define_axes
from dists.freqbins_func import get_bin_strings, get_bins, check_bin_filled, check_create_dir, check_create_file, write_bin_paramfile, execute, get_concat_files
from dists.scores_func import calc_ihs, calc_delihh, calc_xpehh, calc_fst_deldaf, read_neut_normfile, norm_neut_ihs, norm_sel_ihs, norm_neut_xpehh, norm_sel_xpehh, calc_hist_from_scores, write_hists_to_files, get_indices, load_vals_from_files, choose_vals_from_files
from util.parallel import slurm_array #nix this
import argparse
import sys, os
import matplotlib.pyplot as plt

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")

	###############################
	## RUN SELECTION SIMULATIONS ##
	###############################
	get_sel_trajs_parser = subparsers.add_parser('get_sel_trajs', help='Run forward simulations of selection trajectories to populate selscenarios by final allele frequency before running coalescent simulations for entire sample.')
	get_sel_trajs_parser.add_argument('--maxAttempts', action='store', help='maximum number of attempts to generate a trajectory before re-sampling selection coefficient / start time')
	run_neut_sims_parser = subparsers.add_parser('run_neut_sims', help='run additional neutral simulations to create test set')
	run_sel_sims_parser = subparsers.add_parser('run_sel_sims')
	for run_sims_parser in [run_neut_sims_parser, run_sel_sims_parser]:
		run_sims_parser.add_argument('--startrep', type=int, action="store", help='replicate number at which to start', default=1)
		run_sims_parser.add_argument('--paramfolder', type=str, action="store", help='location to read model parameters', default="/idi/sabeti-scratch/jvitti/params/")
	run_sel_sims_parser.add_argument('--trajdir', action="store")
	run_neut_repscores_parser = subparsers.add_parser('run_neut_repscores', help ='run per-replicate component scores')
	run_sel_repscores_parser = subparsers.add_parser('run_sel_repscores', help='run per-replicate component scores')
	for run_repscores_parser in [run_neut_repscores_parser, run_sel_repscores_parser]:
		run_repscores_parser.add_argument('--irep', type=int, action="store", help="replicate number", default=1)
		run_repscores_parser.add_argument('--tpedfolder', type=str, action="store", help="location of input tped data", default="/idi/sabeti-scratch/jvitti/clean/sims/")
		run_repscores_parser.add_argument('--simRecomFile', type=str, action="store", help="location of input recom file", default="/idi/sabeti-scratch/jvitti/params/test_recom.recom")

	run_norm_neut_repscores_parser = subparsers.add_parser('run_norm_neut_repscores')
	run_norm_neut_repscores_parser.add_argument('--edge', type=int, action="store", help="use interior of replicates; define per-end bp. (e.g. 1.5Mb -> 1Mb: 250000)")
	run_norm_neut_repscores_parser.add_argument('--chromlen', type=int, action="store", help="per bp (1.5mb = 1500000)")
	norm_from_binfile_parser = subparsers.add_parser('norm_from_binfile')
	sel_norm_from_binfile_parser = subparsers.add_parser('sel_norm_from_binfile')
	for norm_parser in [run_norm_neut_repscores_parser,norm_from_binfile_parser, sel_norm_from_binfile_parser]:
		norm_parser.add_argument('--score', action='store', default='')
		norm_parser.add_argument('--altpop', action='store', default='')

	for selbin_parser in [run_sel_sims_parser, sel_norm_from_binfile_parser]:
		selbin_parser.add_argument('--selbin', action="store")

	##########################
	### COSI - SHARED ARGS  ##
	##########################
	for cosi_parser in [get_sel_trajs_parser]:
		cosi_parser.add_argument('inputParamFile', action='store', help='file with model specifications for input')
		cosi_parser.add_argument('outputDir', action='store', help='location to write cosi output')
		cosi_parser.add_argument('--cosiBuild', action='store', help='which version of cosi to run', default="coalescent")
		cosi_parser.add_argument('--dropSings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')
		cosi_parser.add_argument('--genmapRandomRegions', action='store_true', help='cosi option to sub-sample genetic map randomly from input')
		cosi_parser.add_argument('n', action='store', type=int, help='num replicates to run') 

	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################
	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='Collate scores from simulated data in order to generate component test probability distributions.')
	if True:
		likes_from_scores_parser.add_argument('neutFile', action='store', help='file with scores for neutral scenarios (normalized if necessary). May be a file containing a list of files (e.g. multiple replicates) if suffix is <.list>.')
		likes_from_scores_parser.add_argument('selFile', action='store', help='file with scores for selected scenarios (normalized if necessary). May be a file containing a list of files (e.g. multiple replicates) if suffix is <.list>.')
		likes_from_scores_parser.add_argument('selPos', action='store', type=int, help='position of causal variant', default=500000)
		likes_from_scores_parser.add_argument('outPrefix', action='store', help='save file as')
		likes_from_scores_parser.add_argument('--thinToSize', action='store_true', help='subsample from simulated SNPs (since nSel << nLinked < nNeut)', default=False)	
		likes_from_scores_parser.add_argument('--ihs', action='store_true', help='define probability distribution for iHS')	
		likes_from_scores_parser.add_argument('--delihh', action='store_true', help='define probability distribution for delIHH')	
		likes_from_scores_parser.add_argument('--nsl', action='store_true', help='define probability distribution for nSl')	
		likes_from_scores_parser.add_argument('--xpehh', action='store_true', help='define probability distribution for XP-EHH')	
		likes_from_scores_parser.add_argument('--deldaf', action='store_true', help='define probability distribution for delDAF')	
		likes_from_scores_parser.add_argument('--fst', action='store_true', help='define probability distribution for Fst')	

		likes_from_scores_parser.add_argument('--edge', type=int, action="store", help="use interior of replicates; define per-end bp. (e.g. 1.5Mb -> 1Mb: 25000)", default=25000)
		likes_from_scores_parser.add_argument('--chromlen', type=int, action="store", help="per bp (1.5mb = 1500000)", default=1500000)

	##############
	## SEL BINS ##
	##############
	for sel_parser in [get_sel_trajs_parser, likes_from_scores_parser]:
		sel_parser.add_argument('--freqRange', type=str, help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
		sel_parser.add_argument('--nBins', type=int, help="number of frequency bins", default=9)

	visualize_likes_parser = subparsers.add_parser('visualize_likes', help='Visualize likelihood tables generated from simulated data.')
	visualize_likes_parser.add_argument('likesFiles', help='input files (comma-delimited, or passed as one file containing paths for each) with bins of probability density function to visualize.')	
	visualize_likes_parser.add_argument('--oldLikes', help='visualize likelihood tables from previous run')

	write_master_likes_parser = subparsers.add_parser('write_master_likes', help="write file for demographic scenario and population pointing to files for all score likelihoods")
	write_master_likes_parser.add_argument('--likes_basedir', default="/idi/sabeti-scratch/jvitti/likes_111516_b/", help="location of likelihood tables ")

	for likes_parser in [likes_from_scores_parser, visualize_likes_parser]:
		likes_parser.add_argument('--nLikesBins', action='store', type=int, help='number of bins to use for histogram to approximate probability density function', default=60)

	for common_parser in [run_neut_sims_parser, run_sel_sims_parser, run_neut_repscores_parser, run_sel_repscores_parser, run_norm_neut_repscores_parser, norm_from_binfile_parser, sel_norm_from_binfile_parser]:
		common_parser.add_argument('--writedir', type =str, help='where to write output', default = "/idi/sabeti-scratch/jvitti/")
		common_parser.add_argument('--checkOverwrite', action="store_true", default=False)
		common_parser.add_argument('--simpop', action='store', help='simulated population', default=1)
		common_parser.add_argument('--model', type=str, default="nulldefault")
		common_parser.add_argument('--nrep', type=int, default=1000)

	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
### Run simuates from specified demographic
### model under various scenarios
def execute_get_sel_trajs(args):
	'''wraps call to run_traj.py to generate forward trajectories of simulated allele frequencies for demographic scenarios with selection'''
	selTrajDir = args.outputDirs
	if selTrajDir[-1] != "/":
		selTrajDir += "/"
	selTrajDir += "sel_trajs/"
	runDir = check_create_dir(selTrajDir)
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)

	print("running " + str(args.n) + " selection trajectories per for each of " + str(args.nBins) + " frequency bins, using model: \n\t\t" + args.inputParamFile)
	print("outputting to " + runDir)

	for ibin in range(args.nBins):
		populateDir = runDir + "sel_" + bin_medians_str[ibin]
		binDir = check_create_dir(populateDir)
		bounds = bin_starts[ibin], bin_ends[ibin]
		paramfilename = populateDir + "/params"
		write_bin_paramfile(args.inputParamFile, paramfilename, bounds)		#rewrite paramfile here? to give rejection sampling 
		if binDir[-1] != "/":
			binDir += "/"
		traj_outputname = binDir + 'rep' + taskIndexStr + '.txt'
		traj_cmd = "python cms/run_traj.py " + traj_outputname + " " + args.cosiBuild + " " + paramfilename + " " + str(args.maxAttempts) +'\n'#
		if args.printOnly:
			dispatch = False
		else:
			dispatch = True
		slurm_array(binDir + "run_sel_traj.sh", traj_cmd, args.n, dispatch=dispatch)

	return #MAKE CLUSTER-INDEPENDENT
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
def execute_run_sel_sims(args):
	'''from run_additional_sims.py'''
	#cmd = "coalescent"
	model = args.model
	nrep = args.nrep
	nStart = args.startrep
	nEnd = nStart + nrep
	writefolder = args.writedir
	trajfolder = args.trajdir
	paramfolder = args.paramfolder
	selpop = args.simpop
	selbin = args.selbin
	
	for foldername in [writefolder, trajfolder, paramfolder]:
		assert(foldername != None)
		if foldername[-1] != "/":
			foldername += "/"
	paramfilename = paramfolder + model + ".par"

	for irep in range(nStart, nEnd + 1):
		outbase = writefolder + "rep" + str(irep)
		trajectory = trajfolder + "rep" + str(irep) + ".txt"
		cmd = "env COSI_NEWSIM=1 env COSI_LOAD_TRAJ=" + trajectory + " coalescent"
		argstring = "-p " + paramfilename + " --genmapRandomRegions --drop-singletons .25 --tped " + outbase + " --output-gen-map"
		#cosicreatefilename = outbase + "_0_1.tped"  	#previous implementation batch-ran sims then batch-renamed them.
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
def execute_run_neut_repscores(args):
	''' adapted from rerun_scores_onerep.py'''
	model = args.model
	pop = args.simpop
	repNum = args.irep
	basedir = args.writedir
	cmsdir = args.cmsdir
	tpeddir = args.tpedfolder
	simRecomFile = args.simRecomFile

	tped = tpeddir + "rep_" + str(repNum) + "_" + str(pop) + ".tped" #temp quick fix
	assert os.path.isfile(tped)
	#in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file  = get_sim_component_score_files(model, repNum, pop, altpop, scenario = "neut", filebase = basedir)
		
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
	proceed = check_create_file(ihs_unnormedfile, args.checkOverwrite)
	if proceed:
		print(ihs_fullcmd)
		execute(ihs_fullcmd)
	delihh_commandstring = "python " + cmsdir + "likes_from_model.py scores_from_sims --delIhh" #WAIT WHAT?
	delihh_unnormedfile =  basedir + "delihh/rep" + str(repNum) + "_" + str(pop) + ".txt"
	delihh_argstring = ihs_unnormedfile + " "+ delihh_unnormedfile
	delihh_fullcmd = delihh_commandstring + " " + delihh_argstring + " --cmsdir " + args.cmsdir
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
		tped2 = tpeddir + "rep_" + str(repNum) + "_" + str(altpop) + ".tped"
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
		fstdeldaf_argumentstring = tped + " --fst_deldaf " + tped2 + " " + fstdeldaf_outfilename + " --recomfile " + simRecomFile + " --cmsdir " + args.cmsdir
		fstdeldaf_fullcmd = fstdeldaf_commandstring + " " + fstdeldaf_argumentstring 
		proceed = check_create_file(fstdeldaf_outfilename, args.checkOverwrite)
		if proceed:
			print(fstdeldaf_fullcmd)
			execute(fstdeldaf_fullcmd)
	return
def execute_run_sel_repscores(args):
	''' adapted from rerun_scores_onerep.py'''
	model = args.model
	pop = args.simpop
	repNum = args.irep
	basedir = args.writedir
	cmsdir = args.cmsdir
	tpeddir = args.tpedfolder
	simRecomFile = args.simRecomFile
	tped = tpeddir + "rep_" + str(repNum) + "_" + str(pop) + ".tped"
	assert os.path.isfile(tped)

	####### Calculate per-population
	####### scores: iHS, delIHH, nSL
	ihs_commandstring = "python " + cmsdir + "scans.py selscan_ihs"
	ihs_outfileprefix = basedir + "ihs/rep" + str(repNum) + "_" + str(pop) 
	ihs_unnormedfile = ihs_outfileprefix + ".ihs.out"
	ihs_argstring = tped + " " + ihs_outfileprefix + " --threads 7 "
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
		tped2 = tpeddir + "rep_" + str(repNum) + "_" + str(altpop) + ".tped"
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
	chrlen, edge = args.chromlen, args.edge
	startbound, endbound = int(edge), chrlen - int(edge) #define replicate interior (conservative strategy to avoid edge effects)

	if score in ['ihs', 'delihh', 'nsl']:
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_"
	elif score in ['xpehh', 'fst']:
		altpop = args.altpop
		concatfilebase = basedir + model + "/neut/concat_" + str(pop) + "_" + str(altpop) + "_"
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
					unnormedfile = basedir + model + "/neut/ihs/rep" + str(irep) + "_" + str(pop) + ".ihs.out"
					physpos_ind = 1
				elif score == "delihh":
					unnormedfile = basedir + model + "/neut/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
					physpos_ind = 1
				elif score == "nsl":
					unnormedfile = basedir + model + "/neut/nsl/rep" + str(irep) + "_" + str(pop) + ".nsl.out"
					physpos_ind = 1
				elif score == "xpehh":
					unnormedfile = basedir + model + "/neut/xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
					physpos_ind = 1
				elif score == "fst":
					unnormedfile = basedir + model + "/neut/fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
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
	
		print(fullcmd)
		execute(fullcmd)
		
	return
def execute_norm_from_binfile(args):
	model = args.model
	pop = args.simpop
	numReps = args.nrep
	basedir = args.writedir
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop
	concatfilename, binfilename = get_concat_files(model, pop, score, altpop, basedir=basedir)
	
	for irep in range(numReps, 0, -1):
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
			normedfilename = unnormedfile + ".norm"
			argstring = "--normalizeXpehh --neutXpehhNormParams " + binfilename + " " + unnormedfile + " " + normedfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['fst']:
			print('currently handling this manually: rewrite_fst_bins.py')
			pass

		alreadyExists = False
		if args.checkOverwrite:
			if not os.path.isfile(normedfilename): #check for overwrite
				alreadyExists = False
			else:
				alreadyExists = True				
		if alreadyExists == False:	
			fullcmd = commandstring + " " + argstring
			print(fullcmd)
			execute(fullcmd)
	return
def execute_sel_norm_from_binfile(args):
	model = args.model
	pop = args.simpop
	numReps = args.nrep
	basedir = args.writedir
	cmsdir = args.cmsdir
	score = args.score
	altpop = args.altpop
	selbin = args.selbin
	concatfilename, binfilename = get_concat_files(model, pop, score, altpop, basedir=basedir)
	
	for irep in range(numReps, 0, -1):
		if score in ['ihs']:
			unnormedfile = basedir  + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/ihs/rep" + str(irep) + "_" + str(pop)  + ".ihs.out"
			normedfilename = unnormedfile + ".norm"
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['delihh']:
			unnormedfile = basedir  + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
			normedfilename = unnormedfile + ".norm"	
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['nsl']:
			unnormedfile = basedir  + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/nsl/rep" + str(irep) + "_" + str(pop)+ ".nsl.out"
			normedfilename = unnormedfile + ".norm"			
			argstring = "--normalizeIhs " + unnormedfile + " " + normedfilename + " --neutIhsNormParams " + binfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['xpehh']:
			unnormedfile = basedir  + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + ".xpehh.out"
			normedfilename = unnormedfile + ".norm"
			argstring = "--normalizeXpehh --neutXpehhNormParams " + binfilename + " " + unnormedfile + " " + normedfilename
			commandstring = "python " + args.cmsdir + "likes_from_model.py scores_from_sims"
		elif score in ['fst']:
			print('currently handling this manually: rewrite_fst_bins.py')
			pass
		fullcmd = commandstring + " " + argstring
		alreadyExists = False
		if args.checkOverwrite:
			if not os.path.isfile(normedfilename): #check for overwrite
				alreadyExists = False
			else:
				alreadyExists = True				
		if alreadyExists == False:
			print(fullcmd)
			execute(fullcmd)
	return	

### Define component score likelihood 
### functions as histograms of simulated scores
def execute_likes_from_scores(args):
	''' define likelihood function for component score based on histograms of score values from simulated data '''
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins) #selfreq bins
	numLikesBins = args.nLikesBins #histogram bins

	chrlen, edge = args.chromlen, args.edge
	startbound, endbound = int(edge), chrlen - int(edge) #define replicate interior (conservative strategy to avoid edge effects)

	datatypes = []
	for status in ['causal', 'linked', 'neutral']:
		for dpoint in ['positions', 'score_final']:
			datatypes.append(status + "_" + dpoint)

	#################
	## LOAD SCORES ##
	#################
	data = {}
	if args.ihs:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_neut = get_indices('ihs', "neut")
		expectedlen_sel, indices_sel = get_indices('ihs', "sel")
		histBins,scoreRange,yLims = get_hist_bins('ihs', numLikesBins)
	elif args.delihh:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_neut = get_indices('delihh', "neut")
		expectedlen_sel, indices_sel = get_indices('delihh', "sel")
		histBins,scoreRange,yLims = get_hist_bins('delihh', numLikesBins)
	elif args.xpehh:
		comparison, stripHeader = True, True
		expectedlen_neut, indices_neut = get_indices('xp', "neut")
		expectedlen_sel, indices_sel = get_indices('xp', "sel")
		histBins,scoreRange,yLims = get_hist_bins('xp', numLikesBins)
	elif args.deldaf:
		comparison, stripHeader = True, True
		expectedlen_neut, indices_neut = get_indices('deldaf', "neut")
		expectedlen_sel, indices_sel = get_indices('deldaf', "sel")
		histBins,scoreRange,yLims = get_hist_bins('deldaf', numLikesBins)
	elif args.fst:
		comparison, stripHeader = True, True
		expectedlen_neut, indices_neut = get_indices('fst', "neut")
		expectedlen_sel, indices_sel = get_indices('fst', "sel")
		histBins,scoreRange,yLims = get_hist_bins('fst', numLikesBins)
	elif args.nsl:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_neut = get_indices('nsl', "neut")
		expectedlen_sel, indices_sel = get_indices('nsl', "sel")
		histBins,scoreRange,yLims = get_hist_bins('nsl', numLikesBins)				
	else:
		print("Must specify score (e.g. --ihs)")
		sys.exit(0)

	xlims = scoreRange	

	val_array = load_vals_from_files(args.neutFile, expectedlen_neut, indices_neut, stripHeader)		
	neut_positions, neut_score_final, neut_anc_freq = val_array[0], val_array[1], val_array[2]

	val_array = load_vals_from_files(args.selFile, expectedlen_sel, indices_sel, stripHeader)		
	sel_positions, sel_score_final, sel_anc_freq = val_array[0], val_array[1], val_array[2]

	#if not args.xpehh and not args.deldaf and not args.fst: 
	#else:
	#	if args.xpehh:
	#		chooseMethod = "max" 
	#		comp = "xpehh"
	#	elif args.deldaf:
	#		chooseMethod = "daf"
	#		comp = "deldaf"
	#	elif args.fst:
	#		chooseMethod = "mean"
	#		comp = "fst"

	#	val_array = choose_vals_from_files(args.neutFile, expectedlen_neut, indices_neut, comp = comp, stripHeader = stripHeader, method=chooseMethod)
	#	neut_positions, neut_score_final = val_array[0], val_array[1]

	#	val_array = choose_vals_from_files(args.selFile, expectedlen_sel, indices_sel, comp = comp, stripHeader = stripHeader, method=chooseMethod)		
	#	sel_positions, sel_score_final = val_array[0], val_array[1]

	causal_indices = [i for i, x in enumerate(sel_positions) if x == args.selPos]
	linked_indices = [i for i, x in enumerate(sel_positions) if x != args.selPos and x > startbound and x < endbound] 

	print("loaded " + str(len(neut_positions)) + " neutral variants from : "  + args.neutFile)
	print(" and  "+str(len(causal_indices)) + " causal variants, and " + str(len(linked_indices)) + " linked variants from: " + args.selFile)

	causal_positions = [sel_positions[variant] for variant in causal_indices]
	causal_score_final = [sel_score_final[variant] for variant in causal_indices]

	linked_positions = [sel_positions[variant] for variant in linked_indices]
	linked_score_final = [sel_score_final[variant] for variant in linked_indices]

	n_causal, n_linked, n_neut, bin_causal, bins_linked, bins_neut = calc_hist_from_scores(causal_score_final, linked_score_final, neut_score_final, xlims, int(numLikesBins), args.thinToSize)
	write_hists_to_files(args.outPrefix, histBins, n_causal, n_linked, n_neut)
	print("wrote to " + args.outPrefix)
	return #REVISIT CHOOSE METHOD?
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
def execute_visualize_likes(args):
	''' view likelihood function (i.e. score histograms) of CMS component scores for all models/populations '''
	likesfilenames = args.likesFiles
	likesfilenames = likesfilenames.split(',')
	if len(likesfilenames) == 1:
		file_paths = []
		openfile=open(likesfilenames[0], 'r')
		for line in openfile:
			file_path = line.strip('\n')
			file_paths.append(file_path)
		openfile.close()
		likesfilenames = file_paths

	print('loading ' + str(len(likesfilenames)) + " likelihood tables..." )

	if args.oldLikes:
		likes_dict = get_old_likes()
	else:		
		likes_dict = {}
		for likesfilename in likesfilenames:
			if not os.path.isfile(likesfilename):
				print('must first run likes_from_scores: ' + likesfilename)
			else:
				starts, ends, vals = read_likes_file(likesfilename)
				key = read_demographics_from_filename(likesfilename) #key = (model, score, dist, pop)
				likes_dict[key] = [starts, ends, vals]
				for item in key:
					if str(item) not in likesfilename:
						print("ERROR: " + str(item) + " " + likesfilename)

	keys = likes_dict.keys()
	scores, pops, models, dists = [], [], [], []
	for key in keys:
		models.append(key[0])
		scores.append(key[1])
		dists.append(key[2])
		pops.append(key[3])

	scores = set(scores)
	scores = list(scores)
	models = set(models)
	models = list(models)
	pops = set(pops)
	pops = list(pops)
	dists = set(dists)
	dists = list(dists)

	print('loaded likes for ' + str(len(scores)) + " scores...")
	print('loaded likes for ' + str(len(models)) + " models...")
	print('loaded likes for ' + str(len(pops)) + " pops...")	
	print('loaded likes for ' + str(len(dists)) + " dists...")	

	colorDict = {'causal':'red', 'linked':'green', 'neut':'blue'}
		
	for score in scores: #each score gets its own figure
		bins, scorerange, ylims = get_hist_bins(score, args.nLikesBins)

		f, axes = plt.subplots(len(models), len(pops), sharex='col', sharey='row')
		f.suptitle("p( component score | demographic model + putative selpop )\n" + score)
		iAxis = 1

		for imodel in range(len(models)): #rows
			model = models[imodel]
			
			for ipop in range(len(pops)): #columns
				
				if len(models) == 1:
					ax = axes[ipop]
				else:
					ax = axes[imodel, ipop]
				pop = pops[ipop]

				causalkey = (model, score, 'causal', pop)
				linkedkey = (model, score, 'linked', pop)
				neutkey = (model, score, 'neut', pop)
				if causalkey in keys:
					starts_causal, ends_causal, vals_causal = likes_dict[causalkey]
					plot_likes(starts_causal, ends_causal, vals_causal, ax, xlims=scorerange, ylims=ylims, color=colorDict['causal'])
				if linkedkey in keys:
					starts_linked, ends_linked, vals_linked = likes_dict[linkedkey]
					plot_likes(starts_linked, ends_linked, vals_linked, ax, xlims=scorerange, ylims=ylims, color=colorDict['linked'])
				if neutkey in keys:
					starts_neut, ends_neut, vals_neut = likes_dict[neutkey]				
					plot_likes(starts_neut, ends_neut, vals_neut, ax, xlims=scorerange, ylims=ylims, color=colorDict['neut'])		
		
				if causalkey in keys and neutkey in keys:
					assert starts_causal == starts_neut
				if neutkey in keys and linkedkey in keys:
					assert starts_neut == starts_linked
				if causalkey in keys and linkedkey in keys:
					assert starts_causal == starts_linked

				iAxis +=1
				if imodel == (len(models)-1):
					ax.set_xlabel(pop)
				if ipop == (len(pops)-1):
					entries = model.split("_")
					label = ""
					for item in entries:
						label += item +"\n"
					label.strip('\n')
					ax.set_ylabel(label, fontsize=7)
					ax.yaxis.set_label_position("right")
		plt.show()	
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
