#!/usr/bin/env python
## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 02.27.2017 	vitti@broadinstitute.org

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from combine.input_func import get_likesfiles_frommaster, write_perpop_ihh_from_xp
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
	freqscores_parser.add_argument('recomFile', type=str, action="store", help="input recombination file") #should work around this		 
	freqscores_parser.add_argument('outfile', type=str, action="store", help="file to write")
	freqscores_parser.add_argument('--modelpath', action='store', type=str, default='cms/cms/model/', help="path to model directory containing executables") #will become redundant with conda

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

	delihh_from_ihs_parser = subparsers.add_parser('delihh_from_ihs', help="Calculate delIHH values from iHS output files.")
	delihh_from_ihs_parser.add_argument('readfile', type=str, action='store', help='input ihs file')
	delihh_from_ihs_parser.add_argument('writefile', type=str, action='store', help='delihh file to write')

	ihh_from_xp_parser = subparser.add_parser('ihh_from_xp', help="extract per-pop iHH values from XP-EHH and write to individual files to facilitate arbitrary population comparisons ")
	ihh_from_xp_parser.add_argument('inXpehh', type=str, help="input xpehh file")
	ihh_from_xp_parser.add_argument('outIhh', type=str, help="write to file")
	ihh_from_xp_parser.add_argument('takePop', default=1, type=int, help="write for first (1) or second (2) pop in XP-EHH file?")

	xp_from_ihh_parser = subparsers.add_parser('xp_from_ihh', help="Calculate XP-EHH based on two per-pop iHH files.")
	xp_from_ihh_parser.add_argument('inIhh1', type=str, action='store', help="input ihh file 1")
	xp_from_ihh_parser.add_argument('inIhh2', type=str, action='store', help="input ihh file 2")
	xp_from_ihh_parser.add_argument('outfilename', type=str, action='store', help="write to file")
	xp_from_ihh_parser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')
	xp_from_ihh_parser.add_argument('--cmsdir', type=str, action='store', help="TEMP; will become redundant with conda packaging", default="/n/home08/jvitti/cms/cms/") 

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
	''' python wrapper for program to calculate Fst and delDAF '''
	calc_fst_deldaf(args.inTped1, args.inTped2, args.recomFile, args.outfile, args.modelpath) #should obviate need for recom file. 
	return
def execute_hapviz(args):
	''' view haplotype data as a colored grid. original Shervin Tabrizi update Joe Vitti ''' 
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
	print("plotted to: " + args.out)
	plt.close()
	return	
def execute_delihh_from_ihs(args):
	''' python wrapper for program to calculate delIHH from iHS file (containing iHH0, iHH1) '''
	calc_delihh(args.readfile, args.writefile)
	return
def execute_ihh_from_xp(args):
	''' extract per-pop iHH scores from XP-EHH file and write to individual files to facilitate arbitrary population comparions '''
	takefile = args.inXpehh
	writefile = args.outIhh
	popNum = args.takePop
	write_perpop_ihh_from_xp(takefile, writefile, popNum)
	return
def execute_xp_from_ihh(args):
	''' python wrapper for program to (re)calculate XP-EHH from per-pop iHH values '''
	inputtped1 = args.inIhh1
	inputtped2 = args.inIhh2
	outfilename = args.outfilename
	cmd = args.cmsdir + "combine/write_xpehh_from_ihh" #will become redundant with conda
	argstring = inputtped1 + " " + inputtped2 + " " + outfilename
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(cmdstring)
	else:
		subprocess.check_call( cmdstring.split() )	
	return
def execute_ml_region(args):
	''' perform within-region localization using machine learning algorithm '''
	#chrom, startBp, endBp = args.chrom, args.startBp, args.endBp
	#IN PROGRESS
	return
def execute_ucsc_viz(args):
	''' write score/position data to file for visualization in UCSC genome browser '''
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
