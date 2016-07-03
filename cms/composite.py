## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 06.28.16 vitti@broadinstitute.org

import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for combining component statistics -- i.e., the final step of the CMS 2.0 pipeline.")
	subparsers = parser.add_subparsers(help="sub-commands")

	###############
	## POP PAIRS ##
	###############
	poppair_parser = subparser.add_parser('poppair', help='collate all component statistics for a given population pair (as a prerequisite to more sophisticated group comparisons')
	poppair_parser.add_argument('in_ihs_file1', help="file with normalized iHS values for putative selpop", action="store")
	poppair_parser.add_argument('in_ihs_file2', help="file with normalized iHS values for altpop", action="store")
	poppair_parser.add_argument('in_ihs_file1', help="file with normalized iHS values for putative selpop", action="store")
	poppair_parser.add_argument('in_delihh_file2', help="file with normalized delIhh values for altpop", action="store")
	poppair_parser.add_argument('in_delihh_file1', help="file with normalized delIhh values for putative selpop", action="store")	
	poppair_parser.add_argument('in_xp_file', help="file with normalized XP-EHH values", action="store") #reversed? 0T 1F
	poppair_parser.add_argument('in_fst_deldaf_file', help="file with Fst and delDAF values", action="store") #reversed? 0T 1F

	#################
	## GENOME-WIDE ##
	#################
	bayesian_gw_parser = subparsers.add_parser('bayesian_gw', help='default algorithm and weighting, genome-wide')
	bayesian_gw_parser.add_argument('inputparamfile', action='store', help='file with specifications for input')

	###################
	## WITHIN REGION ##
	###################
	bayesian_region_parser = subparsers.add_parser('bayesian_region', help='default algorithm and weighting, within-region')
	ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)')

	return parser


#####################
## AUX. FUNCTIONS ###
#####################
def execute_poppair(args):
	cmd = "./combine_cms_scores_poppairs "
	#<chrom> <selpop> <otherpop> <ihs1filename> <delihh1filename> <xpfilename> <xp reversed? 0T 1F> <fst_deldaffilename> <deldaf reversed? 0T 1F>\n");
	print "must finalize connection to combine_cms_scores_poppairs"
	return
def execute_bayesian_gw(args):
	print "must connect composite.py to combine_cms framework"
	return
def execute_bayesian_region(args):
	print "must connect composite.py to combine_cms_regions framework"
	return
def execute_ml_region(args):
	print "must connect composite.py to cms_ml framework"
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_composite()
	args = runparser.parse_args()
	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline