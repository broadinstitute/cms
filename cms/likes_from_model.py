## top-level script for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.03.16 vitti@broadinstitute.org

prefixstring = "{CMS2.0}>>\t\t" #for stderr (make global?)
from dists.func import get_bins, check_make_dir, give_command, check_bin_filled
import argparse
import sys

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
	
	##########################
	### COSI - SHARED ARGS  ##
	##########################
	for cosi_parser in [run_neut_sims_parser, get_sel_trajs_parser, run_sel_sims_parser]:
		cosi_parser.add_argument('inputParamFile', action='store', help='file with model specifications for input')
		cosi_parser.add_argument('outputDir', action='store', help='location to write cosi output')
		cosi_parser.add_argument('--cosiBuild', action='store', help='which version of cosi to run? (*automate installation)', default="/Users/vitti/Desktop/COSI_DEBUG_TEST/cosi-2.0/coalescent")
		cosi_parser.add_argument('--dropSings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')

	##############
	## SEL BINS ##
	##############
	for sel_parser in [get_sel_trajs_parser, run_sel_sims_parser]:
		sel_parser.add_argument('--freqRange', type=str, help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
		sel_parser.add_argument('--nBins', type=int, help="number of frequency bins", default=9)
	
	################################
	## CALCULATE SCORES FROM SIMS ##
	################################
	scores_from_sims_parser = subparsers.add_parser('scores_from_sims', help='get scores from simulations')
	scores_from_sims_parser.add_argument('inputTpedDir', help='directory containing all simulated tpeds from which to calculate')
	scores_from_sims_parser.add_argument('outputScoreDir', help='directory to write output score files')		
	scores_from_sims_parser.add_argument('--ihs', action='store_true')	#includes norm? MAYBE: run ihs is one thing (user specifies neut/sel), norm is another thing, norm_sel is another thing
	scores_from_sims_parser.add_argument('--normalizeIhs', action='store_true')	

	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################

	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='get component score probability distributions from scores')
	likes_from_scores_parser.add_argument('--write', action='store_true', help='once scores have been calculated and normalized, bin values and write probability distributions to file')
	likes_from_scores_parser.add_argument('--plot', action='store_true', help='once probability distributions written to file, visualize')	
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
	if args.outputDir[-1] != "/":
		selTrajDir += "/"
	runDir = check_make_dir(selTrajDir)
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqRange, args.nBins)
	print prefixstring + "running " + str(args.nSimsPerBin) + " selection simulates per for each of " + str(args.nBins) + " frequency bins, using model: \n\t\t" + args.inputParamFile
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
	print prefixstring + "running " + str(args.n) + " sel simulates per scenario from model: " + args.inputParamFile
	return
def execute_scores_from_sims(args):
	''' adapted from JV scores_from_tped_vers.py '''
	inputs = os.listdir(args.inputTpedDir)
	print prefixstring + "calculating from " + str(len(inputs)) + " simulates: " + args.inputTpedDir
	print prefixstring + "writing scores to: " + args.outputDir
	if args.neutIhs is not None:
		inputs = os.listdir(args.inputTpedDir)

		pass
	return
def execute_likes_from_scores(args):
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