## top-level script for demographic modeling as part of CMS 2.0. 
#(CURRENTLY: assumes C programs have been compiled in same directory -- JV must provide Makefile; assumes python in $PATH)
## last updated: 06.29.16 vitti@broadinstitute.org

import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_cms_modeller():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for exploratory fitting of demographic models to population genetic data.")
	subparsers = parser.add_subparsers(help="sub-commands")

	######################
	## CALCULATE TARGET ##
	######################
	target_stats_parser = subparsers.add_parser('target_stats', help='perform per-site(/per-site-pair) calculations of population summary statistics for model target values')
	target_stats_parser.add_argument('tpeds', action='store', help='comma-delimited list of input tped files (only one file per pop being modelled; must run chroms separately or concatenate)',type=list)
	target_stats_parser.add_argument('recom', action='store', help='recombination map') #check consistent with used to create tped
	target_stats_parser.add_argument('--regions', action='store', help='putative neutral regions') #make this optional and/or TSV
	target_stats_parser.add_argument('--freqs', action='store_true', help='calculate summary statistics from within-population allele frequencies') 
	target_stats_parser.add_argument('--ld', action='store_true', help='calculate summary statistics from within-population linkage disequilibrium') 
	target_stats_parser.add_argument('--fst', action='store_true', help='calculate summary statistics from population comparison using allele frequencies') 
	target_stats_parser.add_argument('out', action='store', help='outfile prefix') 
	
	bootstrap_parser = subparsers.add_parser('bootstrap', help='perform bootstrap estimates of population summary statistics in order to finalize model target values')
	bootstrap_parser.add_argument('n', action='store', type=int, help='number of bootstraps to perform in order to estimate standard error of the dataset (should converge for reasonably small n)')

	bootstrap_parser.add_argument('--in-freqs', action='store', help='comma-delimited list of infiles with per-site calculations for population. For bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 
	bootstrap_parser.add_argument('--in-ld', action='store', help='comma-delimited list of infiles with per-site-pair calculations for population. For bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 
	bootstrap_parser.add_argument('--in-fst', action='store', help='comma-delimited list of infiles with per-site calculations for population pair. For bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 	
	bootstrap_parser.add_argument('out', action='store', help='outfile prefix') 
	
	######################
	## VISUALIZE MODEL  ##
	######################

	point_parser = subparsers.add_parser('point', help='run simulates of a point in parameter-space')
	point_parser.add_argument('n', help='num reps', type=int)
	point_parser.add_argument('modelfile', help='demographic model in cosi format')
	point_parser.add_argument('recomfile', help='recombination map')
	point_parser.add_argument('targetvalsfile', help='targetvalsfile for model')	
	
	#########################
	## FIT MODEL TO TARGET ##
	#########################

	grid_parser = subparsers.add_parser('grid', help='run grid search')
	grid_parser.add_argument('grid_inputdimensionsfile', action='store', help='file with specifications of grid search') #must be defined for each search 

	optimize_parser = subparsers.add_parser('optimize', help='run optimization algorithm to fit model parameters')
	optimize_parser.add_argument('optimize_inputdimensionsfile', action='store', help='file with specifications of optimization search to run')

	return parser

#####################
## AUX. FUNCTIONS ###
######################
def execute_target_stats(args):
	'''calls bootstrap_*_popstats_regions to get per-snp/per-snp-pair values; these programs currently have hard-coded arg input -- JV consider switching to getopt'''
	inputtpedstring = ''.join(args.tpeds)
	inputtpeds = inputtpedstring.split(',')
	npops = len(inputtpeds)
	print "calculating summary statistics for " +  str(npops) + " populations..."
	allCmds = []
	for ipop in range(npops):
		inputtped = inputtpeds[ipop]
		if args.freqs:
			freqCmd = ['bootstrap_freq_popstats_regions ', inputtped, args.recom, args.regions, args.out + "_freqs_" + str(ipop)]
			allCmds.append(freqCmd)
		if args.ld:
			ldCmd = ['bootstrap_ld_popstats_regions ', inputtped, args.recom, args.regions, args.out + "_ld_" + str(ipop)]
			allCmds.append(ldCmd)
		if args.fst:
			for jpop in range(ipop+1, npops):
				inputtped2 = inputtpeds[jpop]
				fstCmd = ['bootstrap_fst_popstats_regions ', inputtped, inputtped2, args.recom, args.regions, args.out + "_fst_" + str(ipop) + "_" + str(jpop)]
				allCmds.append(fstCmd)
	for command in allCmds:
		command = [str(x) for x in command]
		#subprocess.check_call( command )
		print command
	return
def execute_bootstrap(args):
	'''pulls all per-snp/per-snp-pair values to get genome-wide bootstrap estimates. adapted from JV experimental: get_neutral_targetstats_from_bootstrap.py'''
	inputestimatefilenames = ''.join(args.ins)
	inputfilenames = inputestimatefilenames.split(',')
	npops = len(inputfilenames)
	print "estimating bootstrap values of summary statistics for " + str(npops) + " populations..."
	allCmds = []
	for ipop in range(npops):
		inputfilename = inputfilenames[ipop]
		bootstrapCmd = 'python get_neutral_targetstats_from_bootstrap.py '

	print "MUST CONNECT PIPELINE from cms_modeller.py TO get_neutral_targetstats_from_bootstrap.py"
	return
def execute_point(args):
	'''runs simulates of a point in parameter-space, comparing to specified target. adapted from JV experimental: grid_point.py'''
	print "generating " + str(n) + " simulations from model: " + args.modelfile
	print "MUST CONNECT PIPELINE from cms_modeller.py TO grid_point.py"
	return
def execute_grid(args):
	''' '''
	print "loading dimensions of grid to search from: " + args.grid_inputdimensionsfile
	print "MUST CONNECT PIPELINE from cms_modeller.py TO grid.py"
	return
def execute_optimize(args):
	print "loading dimensions to search from: " + args.optimize_inputdimensionsfile
	print "MUST CONNECT PIPELINE from cms_modeller.py TO optimize.py (; optimize_out.py...)"
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_cms_modeller()
	args = runparser.parse_args()
	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline