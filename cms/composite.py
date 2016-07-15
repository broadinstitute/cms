## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 07.15.16 vitti@broadinstitute.org

prefixstring = "{CMS2.0}>>\t\t" #for stderr (make global?)
from combine.likes_func import get_likesfiles_frommaster
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
	poppair_parser = subparsers.add_parser('poppair', help='collate all component statistics for a given population pair (as a prerequisite to more sophisticated group comparisons')
	poppair_parser.add_argument('in_ihs_file', help="file with normalized iHS values for putative selpop", action="store")
	poppair_parser.add_argument('in_delihh_file', help="file with normalized delIhh values for putative selpop", action="store")	
	poppair_parser.add_argument('in_xp_file', help="file with normalized XP-EHH values", action="store")
	poppair_parser.add_argument('in_fst_deldaf_file', help="file with Fst, delDaf values for poppair", action="store")
	poppair_parser.add_argument('--xp_reverse_pops', help="include if the putative selpop for outcome is the altpop in XPEHH (and vice versa)", action="store_true")	
	poppair_parser.add_argument('--deldaf_reverse_pops', help="finclude if the putative selpop for outcome is the altpop in delDAF (and vice versa)", action="store_true") #reversed? 0T 1F
	poppair_parser.add_argument('outfile', help="file to write with collated scores", action="store") 

	########################
	## LARGER COMPARISONS ##
	########################
	outgroups_parser = subparsers.add_parser('outgroups', help='combine scores from comparisons of a putative selected pop to 2+ outgroups.')
	outgroups_parser.add_argument('infiles', help="comma-delimited set of pop-pair comparisons", action="store")
	outgroups_parser.add_argument('likesfile', help="text file where probability distributions are specified for component scores", action="store")
	outgroups_parser.add_argument('outfile', help="file to write with finalized scores", action="store") 
	outgroups_parser.add_argument('--region', help="for within-region (rather than genome-wide) CMS", action="store_true") 
	outgroups_parser.add_argument('--chrom', type=str, help="chromosome containing region", action="store") #FOR WITHIN-REGION CMS
	outgroups_parser.add_argument('--startBp', type=int, help="start location of region in basepairs", action="store")
	outgroups_parser.add_argument('--endBp', type=int, help="end location of region in basepairs", action="store")

	#ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)')
	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
def execute_poppair(args):
	cmd = "/combine/combine_scores_poppair"
	if args.xp_reverse_pops is not None:
		xp_reversed = 0
	else:
		xp_reversed = 1
	if args.deldaf_reverse_pops is not None:
		deldaf_reversed = 0
	else:
		deldaf_reversed = 1
	argstring = args.in_ihs_file + " " + args.in_delihh_file + " " + args.in_xp_file + " " + str(xp_reversed) + " " + args.in_fst_deldaf_file + " " + str(deldaf_reversed) + " " + args.outfile 
	cmdstring = cmd + " " + argstring
	print(cmdstring)
	#subprocess.check_output(argstring.split())
	return
def execute_outgroups(args):
	delihh_hit_filename, delihh_miss_filename, ihs_hit_filename, ihs_miss_filename, xpehh_hit_filename, xpehh_miss_filename, fst_hit_filename, fst_miss_filename, deldaf_hit_filename, deldaf_miss_filename = get_likesfiles_frommaster(args.likesfile)
	if not args.region: 	#GENOME-WIDE
		cmd = "/combine/combine_scores_multiplepops"
		argstring = args.outfile + " " + delihh_hit_filename + " " + delihh_miss_filename + " " + ihs_hit_filename + " " + ihs_miss_filename + " " + xpehh_hit_filename + " " + xpehh_miss_filename + " " + fst_hit_filename + " " + fst_miss_filename + " " + deldaf_hit_filename + " " + deldaf_miss_filename 
	else:	#WITHIN REGION
		cmd = "/combine/combine_scores_multiplepops_region"
		argstring = str(args.startBp) + " " + str(args.endBp) + " " + 	argstring = args.outfile + " " + delihh_hit_filename + " " + delihh_miss_filename + " " + ihs_hit_filename + " " + ihs_miss_filename + " " + xpehh_hit_filename + " " + xpehh_miss_filename + " " + fst_hit_filename + " " + fst_miss_filename + " " + deldaf_hit_filename + " " + deldaf_miss_filename 
	for pairfile in args.infiles.split(','):
		argstring += " " + pairfile
	cmdstring = cmd + " " + argstring
	print(cmdstring)
	#subprocess.check_output(argstring.split())
	return
def execute_ml_region(args):
	chrom, startBp, endBp = args.chrom, args.startBp, args.endBp
	print(prefixstring + "must connect composite.py to combine_cms_regions framework")
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_composite()
	args = runparser.parse_args()
	if len(sys.argv) < 2:
		print(prefixstring + "{composite.py}>>\t\t Run with flag -h to view script options.")
		sys.exit()
	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline
