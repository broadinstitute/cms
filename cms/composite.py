#!/usr/bin/env python
## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 11.24.16 vitti@broadinstitute.org

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from combine.recalc_func import interpolate_haps, windows, interpolate_from_windows
from combine.likes_func import get_likesfiles_frommaster
from combine.viz_func import hapSort_coreallele, hapSort, hapViz, readAnnotations, find_snp_index, pullRegion, load_from_hap
from dists.scores_func import calc_fst_deldaf, calc_delihh
import subprocess
import argparse
import gzip
import sys
import os


#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for manipulating and combining component statistics")
	subparsers = parser.add_subparsers(help="sub-commands")

	freqscores_parser = subparsers.add_parser('freqscores', help="Calculate allele frequency-based scores (Fst and delDAF) for a pair of populations.") 
	freqscores_parser.add_argument('inTped1', type=str, action="store", help="input tped 1")			
	freqscores_parser.add_argument('inTped2', type=str, action="store", help="input tped 2")
	freqscores_parser.add_argument('recomFile', type=str, action="store", help="input recombination file")			 
	freqscores_parser.add_argument('outfile', type=str, action="store", help="file to write")
	freqscores_parser.add_argument('--modelpath', action='store', type=str, default='cms/cms/model/', help="path to model directory containing executables")

	hapviz_parser = subparsers.add_parser('hapviz', help="Visualize haplotypes for region")
	hapviz_parser.add_argument('inputfile', type=str, action="store", help="input tped")
	hapviz_parser.add_argument('--startpos',type=int, help="define physical bounds of region")
	hapviz_parser.add_argument('--endpos',type=int, help="define physical bounds of region")
	hapviz_parser.add_argument('out',type=str,default=None, help="save image as file")
	hapviz_parser.add_argument('--corepos',type=int,default=-1, help="partition haplotypes based on allele status at this position")
	hapviz_parser.add_argument('--title', type=str, default=None, help="title to give to plot")
	hapviz_parser.add_argument('--annotate', type=str, default=None, help="tab-delimited file where each line gives <chr.pos>\t<annotation>")			
	hapviz_parser.add_argument('--maf', type=str, default=None, help="filter on minor allele frequency (e.g. .01, .05)")			
	hapviz_parser.add_argument('--dpi', type=str, default=300, help="image resolution")			


	if True:

		win_haps_parser = subparsers.add_parser('win_haps', help='Perform window-based calculation of haplotype scores. ')
		win_haps_parser.add_argument('infilename', type=str, action='store', help='file containing per-site scores')
		win_haps_parser.add_argument('writefilename', type=str, action='store', help='file to write')
		win_haps_parser.add_argument('--windowsize', default=30, type=int, help='number of SNPs per window')
		win_haps_parser.add_argument('--jumplen', default=15, type=int, help='number of SNPs to distance windows')

		interpolate_hapscores_parser = subparsers.add_parser('interpolate_hapscores', help="Fill (otherwise omitted) iHS values at low-freq sites based on sliding window averages.")
		interpolate_hapscores_parser.add_argument('intpedfilename', type=str, action='store', help="input tped")
		interpolate_hapscores_parser.add_argument('inihsfilename', type=str, action='store', help="input per-site score file")
		interpolate_hapscores_parser.add_argument('inwinihsfilename', type=str, action='store', help="input window score file")
		interpolate_hapscores_parser.add_argument('outfilename', type=str, action='store', help="file to write")			

		delihh_from_ihs_parser = subparsers.add_parser('delihh_from_ihs', help="Calculate delIHH values from iHS output files.")
		delihh_from_ihs_parser.add_argument('readfile', type=str, action='store', help='input ihs file')
		delihh_from_ihs_parser.add_argument('writefile', type=str, action='store', help='delihh file to write')

		xp_from_ihh_parser = subparsers.add_parser('xp_from_ihh', help="Calculate XP-EHH based on two per-pop iHH files.")
		xp_from_ihh_parser.add_argument('inIhh1', type=str, action='store', help="input ihh file 1")
		xp_from_ihh_parser.add_argument('inIhh2', type=str, action='store', help="input ihh file 2")
		xp_from_ihh_parser.add_argument('outfilename', type=str, action='store', help="write to file")

		###############
		## POP PAIRS ##
		###############
		poppair_parser = subparsers.add_parser('poppair', help='Collate all component statistics for a given population pair (as a prerequisite to more sophisticated group comparisons).')
		poppair_parser.add_argument('in_ihs_file', type=str, action='store', help="file with normalized iHS values for putative selpop")
		poppair_parser.add_argument('in_nsl_file', type=str, action='store', help="file with normalized nSL values for putative selpop")
		poppair_parser.add_argument('in_delihh_file', type=str, action='store', help="file with normalized delIhh values for putative selpop")	
		poppair_parser.add_argument('in_xp_file', type=str, action='store', help="file with normalized XP-EHH values")
		poppair_parser.add_argument('in_fst_deldaf_file', type=str, action='store', help="file with Fst, delDaf values for poppair")
		poppair_parser.add_argument('--xp_reverse_pops', action="store_true", help="include if the putative selpop for outcome is the altpop in XPEHH (and vice versa)")	
		poppair_parser.add_argument('--fst_deldaf_reverse_pops', action="store_true", help="include if the putative selpop for outcome is the altpop in delDAF (and vice versa)") #reversed? 0T 1F
		poppair_parser.add_argument('outfile', type=str, action='store', help="file to write with collated scores") 

		########################
		## LARGER COMPARISONS ##
		########################
		outgroups_parser = subparsers.add_parser('outgroups', help='Combine scores from comparisons of a putative selected pop to 2+ outgroups.')
		outgroups_parser.add_argument('infiles', type=str, action="store", help="comma-delimited set of pop-pair comparisons")
		outgroups_parser.add_argument('likesfile', type=str, action="store", help="text file where probability distributions are specified for component scores")
		outgroups_parser.add_argument('--likesfile_low', type=str, action="store", help="text file where probability distributions are specified for component scores")
		outgroups_parser.add_argument('--likesfile_mid', type=str, action="store", help="text file where probability distributions are specified for component scores")
		outgroups_parser.add_argument('selpop_likes', type=int, action="store", help="model id for selected pop [1, 2, 3, 4]")


		outgroups_parser.add_argument('outfile', type=str, action="store", help="file to write with finalized scores") 
		outgroups_parser.add_argument('--region', action="store_true", help="for within-region (rather than genome-wide) CMS") 
		outgroups_parser.add_argument('--chrom', type=str, action="store", help="chromosome containing region") #FOR WITHIN-REGION CMS
		outgroups_parser.add_argument('--startBp', type=int, action="store", help="start location of region in basepairs")
		outgroups_parser.add_argument('--endBp', type=int, action="store", help="end location of region in basepairs")

		for common_parser in [xp_from_ihh_parser, poppair_parser, outgroups_parser]:
			common_parser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')

		ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)')
		
		ucsc_viz_parser = subparsers.add_parser('ucsc_viz', help="Generate trackfiles of CMS scores for visualization in the UCSC genome browser.")
		ucsc_viz_parser.add_argument('infile_prefix', type=str, action="store", help="prefix of file containing scores to be reformatted (e.g. 'score_chr' for files named scores_chr#.txt)")
		ucsc_viz_parser.add_argument('outfile', type=str, action="store", help="file to write")
		ucsc_viz_parser.add_argument('--posIndex', type=int, action="store", default=1, help="index for column of datafile containing physical position (zero-indexed)")
		ucsc_viz_parser.add_argument('--scoreIndex', type=str, action="store", default=-2, help="index for column of datafile containing score (zero-indexed)")
		ucsc_viz_parser.add_argument('--strip_header', action="store_true", help="if input files include header line")

	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
def execute_freqscores(args):
	calc_fst_deldaf(args.inTped1, args.inTped2, args.recomFile, args.outfile, args.modelpath)
	return
def execute_hapviz(args):
	"""original Shervin Tabrizi update Joe Vitti"""
	##############
	## LOAD DATA##
	##############
	inputfilename = args.inputfile

	if ".hap" in inputfilename:
		haplotypes, coreindex, physpositions = load_from_hap(inputfilename, args.maf, corePos = args.corepos)
	else:
		if args.startpos is None or args.endpos is None:
			print("must provide bounds with --startpos and --endpos")
			sys.exit(0)
		else:
			startpos = int(args.startpos)
			endpos = int(args.endpos)
			haplotypes, coreindex, physpositions = pullRegion(inputfilename, startpos, endpos, args.maf, corePos = args.corepos)

	print("loaded genotypes for " + str(len(haplotypes[0])) + " sites... ")

	########################
	## SORT BY SIMILARITY ##
	########################
	if args.corepos is not -1:
		hap = hapSort_coreallele(haplotypes, coreindex)
	else:
		hap = hapSort(haplotypes) 

	##########
	## PLOT ##
	##########
	fig = plt.figure()
	ax = fig.add_subplot(111)
	hapViz(ax, hap[0], args.out)

	if args.annotate is not None:
		positions, annotations = readAnnotations(args.annotate)
		ylim = ax.axis()[-1]
		for i_snppos in range(len(positions)):
			snppos = positions[i_snppos]
			annotation = annotations[i_snppos]
			if int(snppos) in physpositions:
				foundindex = physpositions.index(int(snppos))
				ax.plot(foundindex, ylim, "v", color="black", markersize=1)
				ax.plot(foundindex, -5, "^", color="black", markersize=1)
				ax.text(foundindex, -35, str(snppos) +"\n" + annotation, fontsize=2, horizontalalignment='center')

	if args.title is not None:
		plt.title(args.title, fontsize=5)

	plt.tight_layout()
	plt.savefig(args.out, dpi=float(args.dpi))
	plt.close()
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
	calc_delihh(args.readfile, args.writefile)
	return
def execute_xp_from_ihh(args):
	inputtped1 = args.inIhh1
	inputtped2 = args.inIhh2
	outfilename = args.outfilename
	cmd = "dists/write_xpehh_fromihh"
	argstring = inputtped1 + " " + inputtped2 + " " + outfilename
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(cmdstring)
	else:
		subprocess.check_call( cmdstring.split() )	
	return
def execute_poppair(args):
	cmd = "combine/combine_scores_poppair"
	if args.xp_reverse_pops:
		xp_reversed = 0
	else:
		xp_reversed = 1
	if args.fst_deldaf_reverse_pops:
		deldaf_reversed = 0
	else:
		deldaf_reversed = 1
	argstring = args.in_ihs_file + " "  + args.in_nsl_file + " " + args.in_delihh_file + " " + args.in_xp_file + " " + str(xp_reversed) + " " + args.in_fst_deldaf_file + " " + str(deldaf_reversed) + " " + args.outfile 
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(cmdstring)
	else:
		subprocess.check_call( cmdstring.split() )	
	return
def execute_outgroups(args):
	ihs_hit_hi_filename, ihs_miss_hi_filename, nsl_hit_hi_filename, nsl_miss_hi_filename, delihh_hit_hi_filename, delihh_miss_hi_filename, xpehh_hit_hi_filename, xpehh_miss_hi_filename, fst_hit_hi_filename, fst_miss_hi_filename, deldaf_hit_hi_filename, deldaf_miss_hi_filename = get_likesfiles_frommaster(args.likesfile, args.selpop_likes)

	usefreqs = []
	if args.likesfile_low is not None:
		ihs_hit_low_filename, ihs_miss_low_filename, nsl_hit_low_filename, nsl_miss_low_filename, delihh_hit_low_filename, delihh_miss_low_filename, xpehh_hit_low_filename, xpehh_miss_low_filename, fst_hit_low_filename, fst_miss_low_filename, deldaf_hit_low_filename, deldaf_miss_low_filename = get_likesfiles_frommaster(args.likesfile_low, args.selpop_likes)
		usefreqs.append('low')
	if args.likesfile_mid is not None:
		ihs_hit_mid_filename, ihs_miss_mid_filename, nsl_hit_mid_filename, nsl_miss_mid_filename, delihh_hit_mid_filename, delihh_miss_mid_filename,  xpehh_hit_mid_filename, xpehh_miss_mid_filename, fst_hit_mid_filename, fst_miss_mid_filename, deldaf_hit_mid_filename, deldaf_miss_mid_filename = get_likesfiles_frommaster(args.likesfile_mid, args.selpop_likes)
		usefreqs.append('mid')
	while len(usefreqs) < 3:	#HI-FREQ by default
		usefreqs.append('hi')
		
	if not args.region: 	#GENOME-WIDE
		cmd = "combine/combine_scores_multiplepops"
		argstring = args.outfile 
	else:	#WITHIN REGION
		cmd = "combine/combine_scores_multiplepops_region"
		argstring = args.outfile + " " + str(args.startBp) + " " + str(args.endBp) + " "  
	
	for score in ['ihs', 'nsl', 'delihh', 'xpehh', 'fst', 'deldaf']:
		for dist_type in ['hit', 'miss']:
			for freq in usefreqs: #['low', 'mid', 'hi']:
				argument = eval(score + "_" + dist_type + "_" + freq + "_filename")
				argstring += " " + argument 

	for pairfile in args.infiles.split(','):
		argstring += " " + pairfile
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(cmdstring)
	else:
		subprocess.check_call( cmdstring.split() )	
	return
def execute_ml_region(args):
	chrom, startBp, endBp = args.chrom, args.startBp, args.endBp
	#IN PROGRESS
	return
def execute_ucsc_viz(args):
	#convertBedGraph
	inprefix = args.infile_prefix
	outfilename = args.outfile
	outfile = open(outfilename, 'w')
	for chrom in [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 3, 4, 5, 6, 7, 8, 9]: 	#BEDGRAPH MUST BE CASE-SENSITIVE SORTED
		chromfile = inprefix + ".chr" + str(chrom) + ".txt.norm"
		assert os.path.isfile(chromfile)
		infile = open(chromfile, 'r')
		if args.strip_header:
			infile.readline()
		for line in infile:
			entries = line.strip('\n').split()
			startPos = int(entries[int(args.posIndex)])
			score = float(entries[int(args.scoreIndex)]) #use normalized value
			writestring = "chr" + str(chrom) + "\t" + str(startPos) + "\t" + str(startPos + 1) + "\t" + str(score) + "\n"
			outfile.write(writestring)
		infile.close()
	outfile.close()
	print("wrote to: " + outfilename)
	#convertBedGraphtoBigWig:
	print("for large datasets, convert to BigWig format, e.g.: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig\n")
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
	elif len(sys.argv)==2:
		runparser.parse_args([sys.argv[1], '--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline
