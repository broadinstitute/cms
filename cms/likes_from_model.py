## top-level script for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 07.01.16 vitti@broadinstitute.org

import subprocess
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
	run_neut_sims_parser.add_argument('inputparamfile', action='store', help='file with model specifications for input')
	run_neut_sims_parser.add_argument('outputdir', action='store', help='location to write cosi tpeds')
	run_neut_sims_parser.add_argument('cosibuild', action='store', help='which version of cosi to run? (*automate installation)')
	run_neut_sims_parser.add_argument('--dropsings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')

	###############################
	## RUN SELECTION SIMULATIONS ##
	###############################
	get_sel_trajs_parser = subparsers.add_parser('get_sel_trajs', help='run forward simulations of selection trajectories and perform rejection sampling to populate selscenarios by final allele frequency before running coalescent simulations for entire sample')
	get_sel_trajs_parser.add_argument('inputparamfile', action='store', help='file with model specifications for input, including sweep parameters')
	get_sel_trajs_parser.add_argument('--freqrange', help="range of final selected allele frequencies to simulate, e.g. .05-.95", default='.05-.95')
	get_sel_trajs_parser.add_argument('--nbins', type=int, help="number of frequency bins", default=8)
	get_sel_trajs_parser.add_argument('--dropsings', action='store', type=float, help='randomly thin global singletons from output dataset to model ascertainment bias')
	get_sel_trajs_parser.add_argument('nSimsPerBin', type=int, help="number of selection trajectories to generate per allele frequency bin")	
	get_sel_trajs_parser.add_argument('outputdir', action='store', help='location to write selection trajectories')

	run_sel_sims_parser = subparsers.add_parser('run_sel_sims', help='run sel. simulations')	
	run_sel_sims_parser.add_argument('n', action='store', type=int, help='num replicates to run per sel scenario')
	run_sel_sims_parser.add_argument('inputparamfile', action='store', help='file with specifications for input, including sweep')
	#selTRAJ VARIABLES... sel bins... model,

	################################
	## CALCULATE SCORES FROM SIMS ##
	################################
	scores_from_sims_parser = subparsers.add_parser('scores_from_sims', help='get scores from simulations')
	scores_from_sims_parser.add_argument('inputTpedDir', help='directory containing all simulated tpeds from which to calculate')
	scores_from_sims_parser.add_argument('outputDir', help='directory to write output')		
	scores_from_sims_parser.add_argument('--neutIhs', action='store_true')	#includes norm?
	scores_from_sims_parser.add_argument('--selIhs', action='store_true')	
		
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
def get_bins(freqrange, nbins, nSimsPerBin):
	fullrange = [float (x) for x in freqrange.split(',')]
	bin_starts = [fullrange[0] + binlen*i for i in range(numbins+1)]
	bin_ends = [fullrange[0] + binlen*i for i in range(1,numbins+2)]
	bin_medians = [fullrange[0] + (binlen/2) * i*.1 for i in range(numbins+1)]
	bin_medians_str = [str(x)[:3] for x in bin_medians]
	return fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str	
def execute_run_neut_sims(args):
	'''adapted from JV run_sims_from_model_vers.py; previously wrote to ms and ran via UGER taskarrays. adjust per new cosi tped output option'''
	print "running " + str(args.n) + " neutral simulates from model: " + args.inputparamfile
	neutSimCommand = args.cosibuild + " -p " + args.inputparamfile +  " --genmapRandomRegions --output-gen-map " #allow user to pass these?
	if hasattr(args, dropsings):
		neutSimCommand += " --drop-singletons " + args.dropsings
	neutSimCommand += " -n "  + args.n + " -o " + args.outputdir + "/neut"
	print neutSimCommand
	return
def execute_get_sel_trajs(args):
	'''adapted from JV run_sims_from_model_vers.py'''
	fullrange, bin_starts, bin_ends, bin_medians, bin_medians_str = get_bins(args.freqrange, args.nbins, args.nSimsPerBin)
	print "running " + str(args.nSimsPerBin) + " selection simulates per for each of " + str(args.nbins) + " frequency bins, using model: " + args.inputparamfile
	selSimCommand_base = args.cosibuild + " -p " + args.inputparamfile +  " --genmapRandomRegions --output-gen-map " #allow user to pass these?
	if hasattr(args, dropsings):
		selSimCommand_base += " --drop-singletons " + args.dropsings
	print "previously: built in UGER here and used task arrags to index trajectories and perform rejection sampling, must now determine how best to use cos"		
	for ibin in range(args.nbins):
		commandstring = "env COSI_NEWSIM=1 COSI_SAVE_TRAJ={trajectoryname}..."
		unfilled = False#True 
		#loop to populate bin. should make it clear to user that the input for sel coefficients and start times should be distribution?
		#should automate creation of folders? 
		while unfilled == True:
			if thereAreEnoughSimsThatHaveBeenMadeSuccessfully:
				unfilled = False
	return
def execute_run_sel_sims(args):
	print "running " + str(args.n) + " sel simulates per scenario from model: " + args.inputparamfile

	return
def execute_scores_from_sims(args):
	''' adapted from JV scores_from_tped_vers.py '''
	print "calculating from simulates: " + args.inputTpedDir
	print "writing scores to: " + args.outputDir
	if hasattr(args, 'neutIhs'):
		print "calculating ihs..."
	return
def execute_likes_from_scores(args):
	print "MUST CONNECT likes_from_model.py to likes_from_scores_vers.py"
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