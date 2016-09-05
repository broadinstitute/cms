#!/usr/bin/env python
## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 09.05.16 vitti@broadinstitute.org

from combine.recalc_func import write_delIHH_file, interpolate_haps, windows, interpolate_from_windows
from combine.likes_func import get_likesfiles_frommaster
from dists.scores_func import calc_fst_deldaf 
import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for manipulating and combining component statistics")
	subparsers = parser.add_subparsers(help="sub-commands")

	freqscores_parser = subparsers.add_parser('freqscores', help="calculate Fst and delDAF") 
	freqscores_parser.add_argument('inTped1', type=str, action="store", help="input tped 1")			
	freqscores_parser.add_argument('inTped2', type=str, action="store", help="input tped 2")
	freqscores_parser.add_argument('recomFile', type=str, action="store", help="input recombination file")			 
	freqscores_parser.add_argument('outfile', type=str, action="store", help="file to write")

	win_haps_parser = subparsers.add_parser('win_haps', help='perform window-based calculation of haplotype scores')
	win_haps_parser.add_argument('infilename', type=str, action='store', help='file containing per-site scores')
	win_haps_parser.add_argument('writefilename', type=str, action='store', help='file to write')
	win_haps_parser.add_argument('--windowsize', default=30, type=int, help='number of SNPs per window')
	win_haps_parser.add_argument('--jumplen', default=15, type=int, help='number of SNPs to distance windows')
	#win_haps_parser.add_argument('--ihs', action='store_true')
	#win_haps_parser.add_argument('--delihh', action='store_true')

	interpolate_hapscores_parser = subparsers.add_parser('interpolate_hapscores', help="fill haplotype vals at low-freq sites based on sliding window averages")
	interpolate_hapscores_parser.add_argument('intpedfilename', type=str, action='store', help="input tped")
	interpolate_hapscores_parser.add_argument('inihsfilename', type=str, action='store', help="input per-site score file")
	interpolate_hapscores_parser.add_argument('inwinihsfilename', type=str, action='store', help="input window score file")
	interpolate_hapscores_parser.add_argument('outfilename', type=str, action='store', help="file to write")			

	delihh_from_ihs_parser = subparsers.add_parser('delihh_from_ihs')
	delihh_from_ihs_parser.add_argument('readfile', type=str, action='store', help='input ihs file')
	delihh_from_ihs_parser.add_argument('writefile', type=str, action='store', help='delihh file to write')

	xp_from_ihh_parser = subparsers.add_parser('xp_from_ihh', help="calculate XP-EHH based on two per-pop iHH files (ie for computational efficiency)")
	xp_from_ihh_parser.add_argument('inIhh1', type=str, action='store', help="input ihh file 1")
	xp_from_ihh_parser.add_argument('inIhh2', type=str, action='store', help="input ihh file 2")
	xp_from_ihh_parser.add_argument('outfilename', type=str, action='store', help="write to file")

	###############
	## POP PAIRS ##
	###############
	poppair_parser = subparsers.add_parser('poppair', help='collate all component statistics for a given population pair (as a prerequisite to more sophisticated group comparisons')
	poppair_parser.add_argument('in_ihs_file', type=str, action='store', help="file with normalized iHS values for putative selpop")
	poppair_parser.add_argument('in_delihh_file', type=str, action='store', help="file with normalized delIhh values for putative selpop")	
	poppair_parser.add_argument('in_xp_file', type=str, action='store', help="file with normalized XP-EHH values")
	poppair_parser.add_argument('in_fst_deldaf_file', type=str, action='store', help="file with Fst, delDaf values for poppair")
	poppair_parser.add_argument('--xp_reverse_pops', action="store_true", help="include if the putative selpop for outcome is the altpop in XPEHH (and vice versa)")	
	poppair_parser.add_argument('--fst_deldaf_reverse_pops', action="store_true", help="include if the putative selpop for outcome is the altpop in delDAF (and vice versa)") #reversed? 0T 1F
	poppair_parser.add_argument('outfile', type=str, action='store', help="file to write with collated scores") 

	########################
	## LARGER COMPARISONS ##
	########################
	outgroups_parser = subparsers.add_parser('outgroups', help='combine scores from comparisons of a putative selected pop to 2+ outgroups.')
	outgroups_parser.add_argument('infiles', type=str, action="store", help="comma-delimited set of pop-pair comparisons")
	outgroups_parser.add_argument('likesfile', type=str, action="store", help="text file where probability distributions are specified for component scores")
	outgroups_parser.add_argument('outfile', type=str, action="store", help="file to write with finalized scores") 
	outgroups_parser.add_argument('--region', action="store_true", help="for within-region (rather than genome-wide) CMS") 
	outgroups_parser.add_argument('--chrom', type=str, action="store", help="chromosome containing region") #FOR WITHIN-REGION CMS
	outgroups_parser.add_argument('--startBp', type=int, action="store", help="start location of region in basepairs")
	outgroups_parser.add_argument('--endBp', type=int, action="store", help="end location of region in basepairs")

	for common_parser in [xp_from_ihh_parser, poppair_parser, outgroups_parser]:
		common_parser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')

	#ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)')
	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
def execute_freqscores(args):
	calc_fst_deldaf(args.inTped1, args.inTped2, args.recomFile, args.outfile)
	return
def execute_win_haps(args):
	windowsize = args.windowsize
	jumplen = args.jumplen
	infilename = args.infilename
	writefilename = args.writefilename
	windows(windowsize, jumplen, infilename, writefilename)
	return
def execute_interpolate_hapscores(args):
	inputTpedFilename = args.intpedfilename
	inputIhsFilename = args.inihsfilename
	inputWinihsFilename = args.inwinihsfilename
	outputFilename = args.outfilename
	interpolate_from_windows(inputTpedFilename, inputIhsFilename, inputWinihsFilename, outputFilename)
	return
def execute_delihh_from_ihs(args):
	write_delIHH_file(args.readfile, args.writefile)
	return
def execute_xp_from_ihh(args):
	inputtped1 = args.inIhh1
	inputtped2 = args.inIhh2
	outfilename = args.outfilename
	cmd = "/dists/write_xpehh_fromihh"
	argstring = inputtped1 + " " + inputtped2 + " " + outfilename
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(command)
		else:
			subprocess.check_call( cmdstring )	
	return
def execute_poppair(args):
	cmd = "/combine/combine_scores_poppair"
	if args.xp_reverse_pops:
		xp_reversed = 0
	else:
		xp_reversed = 1
	if args.fst_deldaf_reverse_pops:
		deldaf_reversed = 0
	else:
		deldaf_reversed = 1
	argstring = args.in_ihs_file + " " + args.in_delihh_file + " " + args.in_xp_file + " " + str(xp_reversed) + " " + args.in_fst_deldaf_file + " " + str(deldaf_reversed) + " " + args.outfile 
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(command)
		else:
			subprocess.check_call( cmdstring )	
	return
def execute_outgroups(args):
	delihh_hit_filename, delihh_miss_filename, ihs_hit_filename, ihs_miss_filename, xpehh_hit_filename, xpehh_miss_filename, fst_hit_filename, fst_miss_filename, deldaf_hit_filename, deldaf_miss_filename = get_likesfiles_frommaster(args.likesfile)
	if not args.region: 	#GENOME-WIDE
		cmd = "/combine/combine_scores_multiplepops"
		argstring = args.outfile + " " + delihh_hit_filename + " " + delihh_miss_filename + " " + ihs_hit_filename + " " + ihs_miss_filename + " " + xpehh_hit_filename + " " + xpehh_miss_filename + " " + fst_hit_filename + " " + fst_miss_filename + " " + deldaf_hit_filename + " " + deldaf_miss_filename 
	else:	#WITHIN REGION
		cmd = "/combine/combine_scores_multiplepops_region"
		argstring = str(args.startBp) + " " + str(args.endBp) + " " + args.outfile + " " + delihh_hit_filename + " " + delihh_miss_filename + " " + ihs_hit_filename + " " + ihs_miss_filename + " " + xpehh_hit_filename + " " + xpehh_miss_filename + " " + fst_hit_filename + " " + fst_miss_filename + " " + deldaf_hit_filename + " " + deldaf_miss_filename 
	for pairfile in args.infiles.split(','):
		argstring += " " + pairfile
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(command)
		else:
			subprocess.check_call( cmdstring )	
	return
def execute_ml_region(args):
	chrom, startBp, endBp = args.chrom, args.startBp, args.endBp
	print("must connect composite.py to combine_cms_regions framework")
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_composite()
	args = runparser.parse_args()

	# if called with no arguments, print help
	if len(sys.argv)==1:
		runparser.parse_args(['--help'])
	elif len(sys.argv)==2 and (len(commands)>1 or commands[0][0]!=None):
		runparser.parse_args([sys.argv[1], '--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline
