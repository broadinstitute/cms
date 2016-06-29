## top-level script for generating probability distributions for component scores as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_likes_from_model():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for generating probability distributions for component scores from pre-specified demographic model(s).")
	subparsers = parser.add_subparsers(help="sub-commands")

	################################
	## RUN COALESCENT SIMULATIONS ##
	################################
	run_neut_sims_parser = subparsers.add_parser('run_neut_sims', help='run neutral simulations')
	run_neut_sims_parser.add_argument('n', action='store', type=int, help='num replicates to run')
	run_neut_sims_parser.add_argument('inputparamfile', action='store', help='file with model specifications for input')

	run_sel_sims_parser = subparsers.add_parser('run_sel_sims', help='run sel. simulations')	
	run_sel_sims_parser.add_argument('n', action='store', type=int, help='num replicates to run per sel scenario')
	run_sel_sims_parser.add_argument('inputparamfile', action='store', help='file with specifications for input')
	#selTRAJ VARIABLES... sel bins... model,

	################################
	## CALCULATE SCORES FROM SIMS ##
	################################
	scores_from_sims_parser = subparsers.add_parser('scores_from_sims', help='get scores from simulations')
	#INPUT: score + inherited from above

	##################################################
	## GATHER SCORES AND CALCULATE LIKELIHOOD TABLE ##
	##################################################

	likes_from_scores_parser = subparsers.add_parser('likes_from_scores', help='get component score probability distributions from scores')
	#input: inherited
	return parser

#####################
## AUX. FUNCTIONS ###
######################
def execute_run_neut_sims(args):
	print "running " + str(args.n) + " neutral simulates from model: " + args.inputparamfile
	print "MUST CONNECT likes_from_model.py to run_sims_from_model_vers.py"
	return
def execute_run_sel_sims(args):
	print "running " + str(args.n) + " sel simulates per scenario from model: " + args.inputparamfile
	print "MUST CONNECT likes_from_model.py to run_sims_from_model_vers.py"
	return
def execute_scores_from_sims(args):
	print "running " + str(args.n) + " sel simulates per scenario from model: " + args.inputparamfile
	print "MUST CONNECT likes_from_model.py to scores_from_tped_vers.py"
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