## top-level script for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.03.16 vitti@broadinstitute.org

prefixstring = "{CMS2.0}>>\t\t" #for stderr (make global?)
from dists.func import get_bins, check_make_dir, give_command, check_bin_filled
import argparse
import sys, os

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")

	#############################
	## RUN NEUTRAL SIMULATIONS ##
	#############################
	run_neut_sims_parser = subparsers.add_parser('run_neut_sims', help='run neutral simulations')
	run_neut_sims_parser.add_argument('n', action='store', type=int, help='num replicates to run')
	
	###############################
	## RUN SELECTION SIMULATIONS ##
	###############################
	get_sel_trajs_parser = subparsers.add_parser('get_sel_trajs', help='run forward simulations of selection trajectories and perform rejection sampling to populate selscenarios by final allele frequency before running coalescent simulations for entire sample')
	get_sel_trajs_parser.add_argument('nSimsPerBin', type=int, help="number of selection trajectories to generate per allele frequency bin")	
	
	run_sel_sims_parser = subparsers.add_parser('run_sel_sims', help='run sel. simulations')	
	run_sel_sims_parser.add_argument('n', action='store', type=int, help='num replicates to run per sel scenario')
	run_sel_sims_parser.add_argument('trajDir', action='store', type=int, help='location of simulated trajectories (i.e. outputDir from get_sel_trajs)')
	
	##########################
	### COSI - SHARED ARGS  ##
	##########################
	for cosi_parser in [run_neut_sims_parser, get_sel_trajs_parser, run_sel_sims_parser]:
		cosi_parser.add_argument('inputParamFile', action='store', help='file with model specifications for input')
		cosi_parser.add_argument('outputDir', action='store', help='location to write cosi output')
		cosi_parser.add_argument('--cosiBuild', action='store', help='which version of cosi to run? (*automate installation)', default="/Users/vitti/Desktop/COSI_DEBUG_TEST/cosi-2.0/coalescent")
		cosi_parser.add_argument('--dropSings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')

	################################
	## CALCULATE SCORES FROM SIMS ## I'm not sure if we even want to include this. It just wraps infrastructure that already exists elsewhere. leave as exercise for the user?
	################################
	scores_from_sims_parser = subparsers.add_parser('scores_from_sims', help='get scores from simulations')
	#PER POP:
	scores_from_sims_parser.add_argument('inputTped', help='tped from which to calculate score')
	scores_from_sims_parser.add_argument('outputFilename', help='where to write scorefile')		
	scores_from_sims_parser.add_argument('--ihs', action='store_true')
	scores_from_sims_parser.add_argument('--delIhh', action='store_true')
	#POP COMPARISONS:
	scores_from_sims_parser.add_argument('--xpehh', action='store', help="inputTped for altpop")
	scores_from_sims_parser.add_argument('--fst_deldaf', action='store', help="inputTped for altpop")
	#NORMALIZE: 
	scores_from_sims_parser.add_argument('--normalizeIhs', action='store', help="filename for parameters to normalize to; if not given then will by default normalize file to its own global dist")	
	scores_from_sims_parser.add_argument('--normalizeDelIhh', action='store', help="filename for parameters to normalize to; if not given then will by default normalize file to its own global dist")
	scores_from_sims_parser.add_argument('--normalizeXpehh', action='store', help="filename for parameters to normalize to; if not given then will by default normalize file to its own global dist")


	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################
	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='get component score probability distributions from scores')
	likes_from_scores_parser.add_argument('--write', action='store_true', help='once scores have been calculated and normalized, bin values and write probability distributions to file')
	likes_from_scores_parser.add_argument('--plot', action='store_true', help='once probability distributions written to file, visualize')	
	likes_from_scores_parser.add_argument('neutFile', action='store', help='file with scores for neutral scenarios (normalized if necessary)')
	likes_from_scores_parser.add_argument('selFile', action='store', help='')#THIS GONNA HAVE TO BE A DIRECTORY for all the selfreqs? #file with scores for selected scenarios (normalized if necessary)')
	likes_from_scores_parser.add_argument('selPos', action='store', type=int, help='position of causal variant', default=500000)
	likes_from_scores_parser.add_argument('--ihs', action='store_true', help='')	
	likes_from_scores_parser.add_argument('--delihh', action='store_true', help='')	
	likes_from_scores_parser.add_argument('--xp', action='store_true', help='')	
	likes_from_scores_parser.add_argument('--deldaf', action='store_true', help='')	
	likes_from_scores_parser.add_argument('--fst', action='store_true', help='')	

	##############
	## SEL BINS ##
	##############
	for sel_parser in [get_sel_trajs_parser, run_sel_sims_parser, likes_from_scores_parser]:
		sel_parser.add_argument('--freqRange', type=str, help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
		sel_parser.add_argument('--nBins', type=int, help="number of frequency bins", default=9)
	return parser

#####################
## AUX. FUNCTIONS ###
######################
def execute_run_neut_sims(args):
	'''adapted from JV run_sims_from_model_vers.py; previously wrote to ms and ran via UGER taskarrays. adjust per new cosi tped output option'''
	neutRunDir = args.outputDir
	if args.outputDir[-1] != "/":
		neutRunDir += "/"
	neutRunDir += "run_neut_sims"
	runDir = check_make_dir(neutRunDir)
	runDir += "/"
	print prefixstring + "running " + str(args.n) + " neutral simulates from model: " + args.inputParamFile
	print prefixstring + "writing to: " + runDir
	neutSimCommand = args.cosiBuild + " -p " + args.inputParamFile +  " --genmapRandomRegions --output-gen-map " #allow user to pass these?
	if args.dropSings is not None:
		neutSimCommand += " --drop-singletons " + str(args.dropSings)
	neutSimCommand += " -n "  + str(args.n) + " -o " + runDir + "rep"
	print prefixstring + neutSimCommand
	give_command(neutSimCommand)
	return
def execute_get_sel_trajs(args):
	'''adapted from JV run_sims_from_model_vers.py'''
	selTrajDir = args.outputDir
	if selTrajDir[-1] != "/":
		selTrajDir += "/"
	selTrajDir += "sel_trajs/"
	runDir = check_make_dir(selTrajDir)
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
	print prefixstring + "running " + str(args.nSimsPerBin) + " selection trajectories per for each of " + str(args.nBins) + " frequency bins, using model: \n\t\t" + args.inputParamFile
	print prefixstring + "outputting to " + runDir
	selSimCommand_base = args.cosiBuild + " -p " + args.inputParamFile +  " --genmapRandomRegions --output-gen-map " #allow user to pass these?
	if args.dropSings is not None:
		selSimCommand_base += " --drop-singletons " + args.dropSings
	for ibin in range(args.nBins):
		populateDir = runDir + "seltraj_" + bin_medians_str[ibin]
		check_make_dir(populateDir)
		unfilled = True 
		iRep = 1
		while unfilled == True:
			trajectoryname = populateDir + "/rep" + str(iRep) + ".traj"
			commandstring = "env COSI_NEWSIM=1 COSI_SAVE_TRAJ=" + trajectoryname
			fullcommand = commandstring + " " + selSimCommand_base 
			print prefixstring + fullcommand
			give_command(fullcommand)
			#give this command until we get a succesful one?
			#unfilled = False #FOR DEBUG; just loop once
			print "JV: must optimize use of cosi here for rejection sampling!"
			#automatically catch [cosi::tag_error_msg*] = no trajectory found within given number of attempts
			if check_bin_filled(populateDir, args.nSimsPerBin):
				unfilled = False
	return
def execute_run_sel_sims(args):
	'''after get_sel_trajs has been run, use trajectories to generate simulated tped data.'''
	selSimDir = args.outputDir
	if selSimDir[-1] != "/":
		selSimDir += "/"
	selSimDir += "sel_sims/"
	runDir = check_make_dir(selSimDir)
	print prefixstring + "loading trajectories from " + args.trajDir
	print prefixstring + "outputting to " + runDir
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
	for ibin in range(args.nBins):
		trajDir = runDir + "seltraj_" + bin_medians_str[ibin]
		trajectoryFilenames = os.listdir(trajDir)
		print "JV: loaded " + str(len(trajectoryfilenames)) + " files ..."
		print "JV: how best to parallelize this step? use same infrastructure as execute_run_neut?"
	return
def execute_scores_from_sims(args):
	''' adapted from JV scores_from_tped_vers.py. points to scans.py '''
	#inputs = os.listdir(args.inputTpedDir)
	#print prefixstring + "calculating from " + str(len(inputs)) + " simulates: " + args.inputTpedDir
	#print prefixstring + "writing scores to: " + args.outputScoreDir
	inputTped, outputFilename = args.inputTped, args.outputFilename
	if args.ihs:
		calc_ihs(inputTped, outputFilename)
		#for inputTped in inputs:
			#outfilename = args.outputScoreDir + "_ihs"
			#calc_ihs(inputTped, outfilename) #parallelize? 
	if args.delIhh:
		#get ihs file!
		print "get concatenated ihs file..."
		ihsfilename = ""
		calc_delihh(ihsfilename, outputFilename)
	if args.xpehh is not None:
		altinputTped = args.xpehh
		calc_xpehh(inputTped, altinputTped, outputFilename)		
	if args.fst_deldaf is not None:
		altinputTped = args.fst_deldaf
		calc_fst_deldaf(inputTped, altinputTped, outputFilename)
	if args.normalizeIhs is not None:
		#if the normargfile exists, use it. otherwise, just pipe to scans.py
		pass
	if args.normalizeDelIhh is not None:
		pass
	if args.normalizeXpehh is not None:
		pass
	return
def execute_likes_from_scores(args):
	'''adapted from likes_from_scores_vers.py'''
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
	data = {}
	if args.ihs:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "ihs_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "ihs_sel")
	if args.delihh:
		comparison, stripHeader = False, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "delihh_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "delihh_sel")
		pass
	if args.xp:
		comparison, stripHeader = True, True
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "xp_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "xp_sel")
	if args.deldaf:
		comparison, stripHeader = True, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "deldaf_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "deldaf_sel")
	if args.fst:
		comparison, stripHeader = True, False
		expectedlen_neut, indices_sel = get_indices(args.neutFile, "fst_neut")
		expectedlen_sel, indices_sel = get_indices(args.selFile, "fst_sel")
		pass
				#MUST DO THIS FOR EACH FREQBIN
		# OR DO I? CAN I ASSUME USER HAS CONCATENATED IN ORDER TO NORMALIZE...?

	val_array = load_vals_from_files(args.neutFile, expectedlen_neut, indices_neut, stripHeader, verbose)		
	neut_positions, neut_score_final = val_array[0], val_array[1]
	val_array = load_vals_from_files(args.selFile, expectedlen_sel, indices_sel, stripHeader, verbose)		
	sel_positions, sel_score_final = val_array[0], val_array[1]
	# BINS! 

	causal_indices = [i for i, x in enumerate(sel_positions) if x == args.selPos]
	linked_indices = [i for i, x in enumerate(sel_positions) if x != args.selPos] 

	causal_positions = [sel_positions[variant] for variant in causal_indices]
	causal_score_final = [sel_score_final[variant] for variant in causal_indices]
	
	linked_positions = [sel_positions[variant] for variant in linked_indices]
	linked_score_final = [sel_score_final[variant] for variant in linked_indices]
	datatypes = []
	for status in ['causal', 'linked', 'neut']:
		for dpoint in ['positions', 'score_final']:
			datatypes.append(status + "_" + dpoint)
	for datatype in datatypes:
		key = (score, dem_scenario, selpop, altpop, datatype, selfreq) #BINS!
		data[key] = eval(datatype)
	#if args.write:
	#	pass
	if args.plot:
		if comparison: #fst, xpehh, deldaf
			pass
		else: #ihs, delihh
			pass
	print prefixstring + "MUST CONNECT likes_from_model.py to likes_from_scores_vers.py"
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_likes_from_model()
	args = runparser.parse_args()
	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline