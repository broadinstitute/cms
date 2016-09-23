#!/usr/bin/env python
## top-level script for generating probability distributions for component scores as part of CMS 2.0. 
## last updated: 09.23.16 vitti@broadinstitute.org

from dists.likes_func import get_old_likes, read_likes_file, plot_likes, get_hist_bins
from dists.freqbins_func import get_bin_strings, get_bins, check_bin_filled, check_make_dir, write_bin_paramfile
from dists.scores_func import calc_ihs, calc_delihh, calc_xpehh, calc_fst_deldaf, read_neut_normfile, norm_neut_ihs, norm_sel_ihs, norm_neut_xpehh, norm_sel_xpehh, calc_hist_from_scores, write_hists_to_files, get_indices, load_vals_from_files
from util.parallel import slurm_array
import argparse
import subprocess
import sys, os
import matplotlib.pyplot as plt

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")
	
	#############################
	## RUN NEUTRAL SIMULATIONS ##
	#############################
	run_neut_sims_parser = subparsers.add_parser('run_neut_sims', help='Run neutral simulations from a demographic model.')
	
	###############################
	## RUN SELECTION SIMULATIONS ##
	###############################
	get_sel_trajs_parser = subparsers.add_parser('get_sel_trajs', help='Run forward simulations of selection trajectories to populate selscenarios by final allele frequency before running coalescent simulations for entire sample.')
	get_sel_trajs_parser.add_argument('--maxAttempts', action='store', help='maximum number of attempts to generate a trajectory before re-sampling selection coefficient / start time')
	run_sel_sims_parser = subparsers.add_parser('run_sel_sims', help='Run simulations of selection from demographic model and selection trajectories.')	
	run_sel_sims_parser.add_argument('trajDir', action='store', help='location of simulated trajectories (i.e. outputDir from get_sel_trajs)')

	##########################
	### COSI - SHARED ARGS  ##
	##########################
	for cosi_parser in [run_neut_sims_parser, get_sel_trajs_parser, run_sel_sims_parser]:
		cosi_parser.add_argument('inputParamFile', action='store', help='file with model specifications for input')
		cosi_parser.add_argument('outputDir', action='store', help='location to write cosi output')
		cosi_parser.add_argument('--cosiBuild', action='store', help='which version of cosi to run', default="coalescent")
		cosi_parser.add_argument('--dropSings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')
		cosi_parser.add_argument('--genmapRandomRegions', action='store_true', help='cosi option to sub-sample genetic map randomly from input')
		cosi_parser.add_argument('n', action='store', type=int, help='num replicates to run') 

	################################
	## CALCULATE SCORES FROM SIMS ## 
	################################
	scores_from_sims_parser = subparsers.add_parser('scores_from_sims', help='Calculate scores from simulated data.')
	if True:
		scores_from_sims_parser.add_argument('inputFilename', action='store', help='tped from which to calculate score / iHS from which to calculate delIhh / unnormalized to normalize.')
		scores_from_sims_parser.add_argument('outputFilename', action='store', help='where to write scorefile')		
		#PER POP:			
		scores_from_sims_parser.add_argument('--ihs', action='store_true', help="calculate iHS from simulates")
		scores_from_sims_parser.add_argument('--delIhh', action='store_true', help="calculate delIHH from iHS from simulates")
		#POP COMPARISONS:
		scores_from_sims_parser.add_argument('--xpehh', action='store', help="inputTped for altpop")
		scores_from_sims_parser.add_argument('--fst_deldaf', action='store', help="inputTped for altpop")
		scores_from_sims_parser.add_argument('--recomfile', action='store', help="input recomfile for sims")
		#NORMALIZE: 
		#scores_from_sims_parser.add_argument('--concatNeutScores', action='store', help="directory from which to concatenate neutral scores (in order to normalize)")
		scores_from_sims_parser.add_argument('--normalizeIhs', action='store_true')
		scores_from_sims_parser.add_argument('--neutIhsNormParams', action='store', help="filename for parameters to normalize to; if flag not given then will by default normalize file to its own global dist")
		scores_from_sims_parser.add_argument('--normalizeDelIhh', action='store_true')	
		scores_from_sims_parser.add_argument('--neutDelIhhNormParams', action='store', help="filename for parameters to normalize to; if flag not given then will by default normalize file to its own global dist")
		scores_from_sims_parser.add_argument('--normalizeXpehh', action='store_true')
		scores_from_sims_parser.add_argument('--neutXpehhNormParams', help="filename for parameters to normalize to; if flag not given then will by default normalize file to its own global dist")

	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################
	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='Collate scores from simulated data in order to generate component test probability distributions.')
	if True:
		#likes_from_scores_parser.add_argument('--write', action='store_true', help='once scores have been calculated and normalized, bin values and write probability distributions to file')
		#likes_from_scores_parser.add_argument('--plot', action='store_true', help='once probability distributions written to file, visualize')	
		likes_from_scores_parser.add_argument('neutFile', action='store', help='file with scores for neutral scenarios (normalized if necessary)')
		likes_from_scores_parser.add_argument('selFile', action='store', help='file with scores for selected scenarios (normalized if necessary)')
		likes_from_scores_parser.add_argument('selPos', action='store', type=int, help='position of causal variant', default=500000)
		likes_from_scores_parser.add_argument('nLikesBins', action='store', type=int, help='number of bins to use for histogram to approximate probability density function', default=60)
		likes_from_scores_parser.add_argument('outPrefix', action='store', help='save file as')
		likes_from_scores_parser.add_argument('--thinToSize', action='store_true', help='subsample from simulated SNPs (since nSel << nLinked < nNeut)')	
		likes_from_scores_parser.add_argument('--ihs', action='store_true', help='define probability distribution for iHS')	
		likes_from_scores_parser.add_argument('--delihh', action='store_true', help='define probability distribution for delIHH')	
		likes_from_scores_parser.add_argument('--xp', action='store_true', help='define probability distribution for XP-EHH')	
		likes_from_scores_parser.add_argument('--deldaf', action='store_true', help='define probability distribution for delDAF')	
		likes_from_scores_parser.add_argument('--fst', action='store_true', help='define probability distribution for Fst')	

	##############
	## SEL BINS ##
	##############
	for sel_parser in [get_sel_trajs_parser, run_sel_sims_parser, scores_from_sims_parser, likes_from_scores_parser]:
		sel_parser.add_argument('--freqRange', type=str, help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
		sel_parser.add_argument('--nBins', type=int, help="number of frequency bins", default=9)

	visualize_likes_parser = subparsers.add_parser('visualize_likes', help='Visualize likelihood tables generated from simulated data.')
	#visualize_likes_parser.add_argument('')

	for common_parser in [run_neut_sims_parser]:#, get_sel_trajs_parser, run_sel_sims_parser, scores_from_sims_parser, likes_from_scores_parser, visualize_likes_parser]:
		common_parser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')

	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
def execute_run_neut_sims(args):
	'''run neutral simulations from model'''
	neutRunDir = args.outputDir
	if args.outputDir[-1] != "/":
		neutRunDir += "/"
	neutRunDir += "neut_sims"
	runDir = check_make_dir(neutRunDir)
	runDir += "/"
	print("running " + str(args.n) + " neutral simulates from model: " + args.inputParamFile)
	print("writing to: " + runDir)
	neutSimCommand = args.cosiBuild + " -p " + args.inputParamFile + " --output-gen-map " 
	if args.genmapRandomRegions:
		neutSimCommand += " --genmapRandomRegions"
	if args.dropSings is not None:
		neutSimCommand += " --drop-singletons " + str(args.dropSings)
	neutSimCommand += " -n "  + str(args.n) + " --tped " + runDir + "rep"
	if args.printOnly:
		print(neutSimCommand)
	else:
		subprocess.check_output(neutSimCommand.split())
	return
def execute_get_sel_trajs(args):
	'''generate selection trajectories as needed'''
	selTrajDir = args.outputDir
	if selTrajDir[-1] != "/":
		selTrajDir += "/"
	selTrajDir += "sel_trajs/"
	runDir = check_make_dir(selTrajDir)
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)

	print("running " + str(args.n) + " selection trajectories per for each of " + str(args.nBins) + " frequency bins, using model: \n\t\t" + args.inputParamFile)
	print("outputting to " + runDir)

	for ibin in range(args.nBins):
		populateDir = runDir + "sel_" + bin_medians_str[ibin]
		binDir = check_make_dir(populateDir)
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

	return
def execute_run_sel_sims(args):
	'''after get_sel_trajs has been run, use trajectories to generate simulated tped data.'''
	selSimDir = args.outputDir
	if selSimDir[-1] != "/":
		selSimDir += "/"
	selSimDir += "sel_sims/"
	runDir = check_make_dir(selSimDir)
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)

	print("loading trajectories from " + args.trajDir)
	print("outputting to " + runDir)

	for ibin in range(args.nBins):
		binTrajDir = args.trajDir + "sel_" + bin_medians_str[ibin]
		populateDir = runDir + "sel_" + bin_medians_str[ibin]
		populateDir = check_make_dir(populateDir)
		trajectoryFilenames = os.listdir(binTrajDir)
		paramfilename = binTrajDir + "/params"
		assert os.path.isfile(paramfilename)			
		if binTrajDir [-1] != "/":
			binTrajDir += "/"
		sim_cmd = "env COSI_NEWSIM=1 env COSI_LOAD_TRAJ=" + binTrajDir + "/rep" + taskIndexStr + ".txt " + args.cosiBuild, + " -p " + paramfilename + " --output-gen-map " 
		if args.genmapRandomRegions:
			sim_cmd += " --genmapRandomRegions"
		if args.dropSings is not None:
			sim_cmd += " --drop-singletons " + str(args.dropSings)
		sim_cmd += " --tped " + populateDir + 'rep' + taskIndexStr
		if args.printOnly:
			dispatch = False
		else:
			dispatch = True
		slurm_array(runDir + "run_sel_sims.sh", sim_cmd, args.n, dispatch=dispatch)
	return
def execute_scores_from_sims(args):
	''' adapted from JV scores_from_tped_vers.py. functions point to scans.py'''
	inputFilename, outputFilename = args.inputFilename, args.outputFilename

	#################
	## CALC SCORES ##
	#################
	if args.ihs:
		calc_ihs(inputFilename, outputFilename)
	if args.delIhh:
		ihsfilename = inputFilename
		calc_delihh(ihsfilename, outputFilename)
	if args.xpehh is not None:
		altinputTped = args.xpehh
		calc_xpehh(inputFilename, altinputTped, outputFilename)		
	if args.fst_deldaf is not None:
		altinputTped = args.fst_deldaf
		calc_fst_deldaf(inputFilename, altinputTped, args.recomfile, outputFilename, 'cms/cms/model/')

	###############
	## NORMALIZE ##
	###############
	#if args.concatNeutScores is not None:
	#	neutSimDir = args.concatNeutScores
	#	neutSimDir_contents = os.listdir(neutSimDir)
	#	neutSim_scores = [item for item in neutSimDir_contents if ".out" in item]
	#	print('loading neutral simulates from ' str(len(neutSim_scores)) + " files from " + neutSimDir + ' ...')

	#	neut_xpehh = [item for item in neutSim_scores if "xpehh" in item]
	#	neut_ihs = [item for item in neutSim_scores if "ihs" in item]

		#sort by poppairs
		#if len(neut_xpehh) > 0:
		#	pass

		#sort by pops
		#if len(neut_ihs) > 0:
		#	pass

	if args.normalizeIhs:
		if args.neutIhsNormParams is not None:
			fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
			print("loading normalization parameters from " + args.neutIhsNormParams + " ...")
			norm_sel_ihs(args.inputFilename, args.neutIhsNormParams)
		else:
			norm_neut_ihs(args.inputFilename, args.outputFilename)

	if args.normalizeDelIhh:
		if args.neutDelIhhNormParams is not None:
			fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
			print("loading normalization parameters from " + args.neutDelIhhNormParams + " ...")		
			norm_sel_ihs(args.inputFilename, args.neutDelIhhNormParams)
		else:
			norm_neut_ihs(args.inputFilename, args.outputFilename)

	if args.normalizeXpehh:
		if args.normalizeXpehh:
			if args.neutXpehhNormParams is not None:
				norm_sel_xpehh(args.inputFilename, args.neutXpehhNormParams)
			else:
				norm_neut_xpehh(args.inputFilename, args.outputFilename)
		
	return
def execute_likes_from_scores(args):
	'''adapted from likes_from_scores_vers.py'''
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins) #selfreq bins
	numLikesBins = args.nLikesBins #histogram bins
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
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "ihs_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "ihs_sel")
		histBins,scoreRange,yLims = get_hist_bins('ihs', numLikesBins)
	if args.delihh:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "delihh_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "delihh_sel")
		histBins,scoreRange,yLims = get_hist_bins('delihh', numLikesBins)
	if args.xp:
		comparison, stripHeader = True, True
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "xp_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "xp_sel")
		histBins,scoreRange,yLims = get_hist_bins('xp', numLikesBins)
	if args.deldaf:
		comparison, stripHeader = True, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "deldaf_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "deldaf_sel")
		histBins,scoreRange,yLims = get_hist_bins('deldaf', numLikesBins)
	if args.fst:
		comparison, stripHeader = True, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "fst_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "fst_sel")
		histBins,scoreRange,yLims = get_hist_bins('fst', numLikesBins)
		
	val_array = load_vals_from_files(args.neutFile, expectedlen_neut, indices_neut, stripHeader, verbose)		
	neut_positions, neut_score_final, neut_anc_freq = val_array[0], val_array[1], val_array[2]
	val_array = load_vals_from_files(args.selFile, expectedlen_sel, indices_sel, stripHeader, verbose)		
	sel_positions, sel_score_final, sel_anc_freq = val_array[0], val_array[1], val_array[2]

	neut_der_freq = [1. - float(x) for x in neut_anc_freq]
	sel_der_freq = [1. - float(x) for x in sel_anc_freq]
	for ibin in range(len(bin_starts)):
		causal_indices = [i for i, x in enumerate(sel_positions) if x == args.selPos and sel_der_freq[i]>=bin_starts[ibin] and sel_der_freq[i]<bin_ends[ibin]]
		linked_indices = [i for i, x in enumerate(sel_positions) if x != args.selPos and sel_der_freq[i]>=bin_starts[ibin] and sel_der_freq[i]<bin_ends[ibin]] 
		neutral_indices = [i for i, x in enumerate(neut_positions) if neut_der_freq[i]>=bin_starts[ibin] and neut_der_freq[i]<bin_ends[ibin]] 

		causal_positions = [sel_positions[variant] for variant in causal_indices]
		causal_score_final = [sel_score_final[variant] for variant in causal_indices]
		
		linked_positions = [sel_positions[variant] for variant in linked_indices]
		linked_score_final = [sel_score_final[variant] for variant in linked_indices]

		neutral_positions = [neut_positions[variant] for variant in neutral_indices]
		neutral_score_final = [neut_score_final[variant] for variant in neutral_indices]
		
		for datatype in datatypes:
			key = (bin_medians_str[ibin], datatype) #BINS!
			data[key] = eval(datatype)
	
	#################
	## BIN SCORES ##
	#################		
	#build in option to merge bins here?
	for ibin in range(len(bin_starts)):
		xlims = scoreRange
		causal_scores, linked_scores, neut_scores = data[(bin_medians_str[ibin], 'causal_scores_final')], data[(bin_medians_str[ibin], 'linked_scores_final')], data[(bin_medians_str[ibin], 'neut_scores_final')]
		n_causal, n_linked, n_neut, bin_causal, bins_linked, bins_neut = calc_hist_from_scores(causal_scores, linked_scores, neut_scores, xlims, args.thinToSize)
		write_hists_to_files(args.outPrefix +"_" + bin_medians_str[ibin], histBins, n_causal, n_linked, n_neut)
	return
def execute_visualize_likes(args):
	'''currently: view all'''
	old_likes = get_old_likes()
	scores = ['ihs', 'delihh', 'fst', 'xp', 'deldaf']
	pops = range(1,5)
	models = ['default_112115_825am', 'default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault', 'nulldefault_constantsize']
	dists = ['causal', 'linked', 'neut']
	colorDict = {'causal':'red', 'linked':'green', 'neut':'blue'}
	
	for score in scores: #each score gets its own figure
		#fig = plt.figure()
		#fig, axarr = plt.subplots(len(models),len(pops), sharex = True, sharey =False)#len(models), len(pops), sharex=False, sharey=True)
		bins, scorerange, ylims = get_hist_bins(score, 50)
		f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7 ,ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16), (ax17, ax18,ax19, ax20)) = plt.subplots(5, 4, sharex='col', sharey='row')
		f.suptitle("p( component score | demographic model + putative selpop )\n" + score)
		iAxis = 1

		for imodel in range(len(models)): #rows
			model = models[imodel]
			for ipop in range(len(pops)): #columns
				ax = eval('ax' + str(iAxis))
				pop = pops[ipop]
				#ax = fig.add_subplot(111)#111)	
				#print(str(imodel)  + "\t" + str(ipop))
				#ax = axarr[imodel, ipop]			
				starts_causal, ends_causal, vals_causal = old_likes[(model, score, 'causal', pop)]
				starts_linked, ends_linked, vals_linked = old_likes[(model, score, 'linked', pop)]
				starts_neut, ends_neut, vals_neut = old_likes[(model, score, 'neut', pop)]				
				assert starts_causal == starts_neut and starts_neut == starts_linked
				plot_likes(starts_causal, ends_causal, vals_causal, ax, xlims=scorerange, ylims=ylims, color=colorDict['causal'])
				plot_likes(starts_linked, ends_linked, vals_linked, ax, xlims=scorerange, ylims=ylims, color=colorDict['linked'])
				plot_likes(starts_neut, ends_neut, vals_neut, ax, xlims=scorerange, ylims=ylims, color=colorDict['neut'])		
				iAxis +=1
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
