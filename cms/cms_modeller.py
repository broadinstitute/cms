## top-level script for demographic modeling as part of CMS 2.0. (assumes C programs have been compiled in same directory)
## last updated: 06.28.16 vitti@broadinstitute.org

import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################

def full_parser_cms_modeller():
	"""parser execution/documentation for this script"""
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for exploratory fitting of demographic models to population genetic data.")
	subparsers = parser.add_subparsers(help="sub-commands")

	######################
	## CALCULATE TARGET ##
	#######################
	target_stats_parser = subparsers.add_parser('target_stats', help='perform per-site(/per-site-pair) calculations of population summary statistics for model target values')
	target_stats_parser.add_argument('tpeds', action='store', help='comma-delimited list of input tped files (i.e., for all populations to jointly model)')
	target_stats_parser.add_argument('recom', action='store', help='recombination map') #check consistent with used to create tped
	target_stats_parser.add_argument('regions', action='store', help='putative neutral regions') #make this optional and/or TSV
	target_stats_parser.add_argument('out', action='store', help='outfile prefix') 
	
	bootstrap_parser = subparsers.add_parser('bootstrap', help='perform bootstrap estimates of population summary statistics in order to finalize model target values')
	#bootstrap_parser.add_argument('n', action='store', type=int, help='number of bootstraps to perform in order to estimate standard error of the dataset (should converge for reasonably small n)')
	#bootstrap_parser.add_argument('in', action='store', help='infile with per-site(/per-site-pair) calculations') 
	#bootstrap_parser.add_argument('out', action='store', help='outfile prefix') 
	#INPUT: nbootstrapreps, outputfile, -g-w- inputFILES from above, currently get_neutral_targetstats_from_bootstrap.py includes func to point to outfiles from above. rely on user to concatenate them manually?

	######################
	## VISUALIZE MODEL  ##
	######################

	point_parser = subparsers.add_parser('point', help='run simulates of a point in parameter-space')
	#point_parser.add_argument('modelfile', help='demographic model in cosi format')
	#point_parser.add_argument('recomfile', help='recombination map')
	#point_parser.add_argument('targetvalsfile', help='targetvalsfile for model')	
	#INPUT: modelfile, targetvalsfile, writelocation for graphics

	#########################
	## FIT MODEL TO TARGET ##
	#########################

	grid_parser = subparsers.add_parser('grid', help='run grid search')
	#grid_parser.add_argument('inputparamfile', action='store', help='file with specifications of grid search')
	#do from command line? dimensions to vary etc?

	#optimize_parser = subparsers.add_parser('optimize', help='run optimization algorithm to fit model parameters')
	#optimize_parser.add_argument('inputparamfile', action='store', help='file with specifications of optimization search to run')

	return parser

#####################
## AUX. FUNCTIONS ###
######################

def execute_target_stats(args):
	npops = len(args.tpeds)
	print "calculating summary statistics for " +  str(npops) + " populations..."
	allCmds = []
	for ipop in range(npops):
		inputtped = args.tpeds[ipop]
		freqCmd = ['bootstrap_freq_popstats_regions ', inputtped, args.recom, args.regions, args.out + "_freqs_" + str(ipop)]
		ldCmd = ['bootstrap_ld_popstats_regions ', inputtped, args.recom, args.regions, args.out + "_ld_" + str(ipop)]
		allCmds.extend([freqCmd, ldCmd])
		for jpop in range(ipop+1, npops):
			inputtped2 = args.tpeds[jpop]
			fstCmd = ['bootstrap_fst_popstats_regions ', inputtped1, inputtped2, args.recom, args.regions, args.out + "_fst_" + str(ipop) + "_" + str(jpop)]
			allCmds.append(fstCmd)
	for command in allCmds:
		command = [str(x) for x in command]
		subprocess.check_call( command )
	return
def execute_bootstrap(args):
	return
def execute_point(args):
	return
def execute_grid(args):
	return

##########
## MAIN ##
##########

if __name__ == '__main__':
#def main():
	print "halleloo"
	runparser = full_parser_cms_modeller()
	#args = runparser.parse_args()
	subcommand = sys.argv[1]
	#eval('execute_' + subcommand + "(args)") #points to functions defined above, which wrap other programs in the pipeline
	
#main()