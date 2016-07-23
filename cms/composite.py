## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 07.23.16 vitti@broadinstitute.org

prefixstring = "{CMS2.0}>>\t\t" #for stderr (make global?)
from combine.recalc_func import write_delIHH_file, interpolate_haps, calc_winIhs
from combine.likes_func import get_likesfiles_frommaster
from dists.scores_func import calc_fst_deldaf
import subprocess
import argparse
import sys

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_composite():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for combining component statistics -- i.e., the final step of the CMS 2.0 pipeline.")
	subparsers = parser.add_subparsers(help="sub-commands")

	freqscores_parser = subparsers.add_parser('freqscores', help="calculate Fst and delDAF") #arguably, this (and some of the below?) should live in scans.py.
	freqscores_parser.add_argument('inTped1', help="input tped 1", action="store")			 #trying to avoid desecrating CT's code for now. 
	freqscores_parser.add_argument('inTped2', help="input tped 2", action="store")
	freqscores_parser.add_argument('recomFile', help="input recombination file", action="store")			 #trying to avoid desecrating CT's code for now. 
	freqscores_parser.add_argument('outfile', help="file to write", action="store")

	win_haps_parser = subparsers.add_parser('win_haps', help='perform window-based calculation of haplotype scores')
	win_haps_parser.add_argument('infilename', action='store')
	win_haps_parser.add_argument('writefilename', action='store')
	win_haps_parser.add_argument('--windowsize', default=30, type=int)
	win_haps_parser.add_argument('--jumplen', default=15, type=int)
	#win_haps_parser.add_argument('--ihs', action='store_true')
	#win_haps_parser.add_argument('--delihh', action='store_true')

	interpolate_hapscores_parser = subparsers.add_parser('interpolate_hapscores', help="fill haplotype vals at low-freq sites based on sliding window averages")
	interpolate_hapscores_parser.add_argument('intpedfilename', action='store')
	interpolate_hapscores_parser.add_argument('inihsfilename', action='store')
	interpolate_hapscores_parser.add_argument('inwinihsfilename', action='store')
	interpolate_hapscores_parser.add_argument('outfilename', action='store')			

	delihh_from_ihs_parser = subparsers.add_parser('delihh_from_ihs')
	delihh_from_ihs_parser.add_argument('readfile', help='input ihs file', action='store')
	delihh_from_ihs_parser.add_argument('writefile', help='delihh file to write', action='store')

	xp_from_ihh_parser = subparsers.add_parser('xp_from_ihh', help="calculate XP-EHH based on two per-pop iHH files (ie for computational efficiency)")
	xp_from_ihh_parser.add_argument('inIhh1', help="input ihh file 1", action="store")
	xp_from_ihh_parser.add_argument('inIhh2', help="input ihh file 2", action="store")
	xp_from_ihh_parser.add_argument('outfilename', help="write to file", action="store")

	###############
	## POP PAIRS ##
	###############
	poppair_parser = subparsers.add_parser('poppair', help='collate all component statistics for a given population pair (as a prerequisite to more sophisticated group comparisons')
	poppair_parser.add_argument('in_ihs_file', help="file with normalized iHS values for putative selpop", action="store")
	poppair_parser.add_argument('in_delihh_file', help="file with normalized delIhh values for putative selpop", action="store")	
	poppair_parser.add_argument('in_xp_file', help="file with normalized XP-EHH values", action="store")
	poppair_parser.add_argument('in_fst_deldaf_file', help="file with Fst, delDaf values for poppair", action="store")
	poppair_parser.add_argument('--xp_reverse_pops', help="include if the putative selpop for outcome is the altpop in XPEHH (and vice versa)", action="store_true")	
	poppair_parser.add_argument('--fst_deldaf_reverse_pops', help="finclude if the putative selpop for outcome is the altpop in delDAF (and vice versa)", action="store_true") #reversed? 0T 1F
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
def execute_freqscores(args):
	#cmd = "/model/calc_fst_deldaf"
	#argstring = inputTped1 + " " + inputTped2 + " " + recomFile + " " + outputFile
	#cmdstring = cmd + " " + argstring
	#print(cmdstring)
	#subprocess.check_output(argstring.split())
	calc_fst_deldaf(args.inTped1, args.inTped2, args.recomFile, args.outfile)
	return
def execute_win_haps(args):
	'''from winiHS.py (/windelihh?)'''
	#if args.ihs:	
	#if args.delihh:
	posindex, scoreindex = 1, 6
	ihh1_index, ihh0_index, freq1_index = 3,4,2

	windowsize = args.windowsize
	jumplen = args.jumplen
	infilename = args.infilename
	writefilename = args.writefilename

	positions, scores = [], []
	ihh1s, ihh0s, freq1s = [], [], []
	infile = open(infilename, 'r')
	for line in infile:
		entries = line.split()
		position, score = int(entries[posindex]), float(entries[scoreindex])
		positions.append(position)
		scores.append(score)

		ihh1, ihh0, freq1 = float(entries[ihh1_index]), float(entries[ihh0_index]), float(entries[freq1_index])
		ihh1s.append(ihh1)
		ihh0s.append(ihh0)
		freq1s.append(freq1) #not sure that this even counts for anything. 

	infile.close()

	###################
	## DEFINE WINDOWS #
	#this is what I interpret 
	#Schlebusch to be doing (?)

	numSnps = len(scores)
	print(str(numSnps) + " scores available")
	all_phys_windows, all_score_windows = [], []
	all_ihh1_windows, all_ihh0_windows, all_freq1_windows = [], [], []
	iSnp = 0
	while iSnp < numSnps:
		window_phys, window_scores = [], []
		window_ihh1, window_ihh0, window_freq1 = [], [], []
		while len(window_scores) < windowsize:
			if iSnp < numSnps: #?
				window_phys.append(positions[iSnp])
				window_scores.append(scores[iSnp])
				window_ihh1.append(ihh1s[iSnp])
				window_ihh0.append(ihh0s[iSnp])
				window_freq1.append(freq1s[iSnp])
			else:
				break
			iSnp +=1
		all_phys_windows.append(window_phys)
		all_score_windows.append(window_scores)
		all_ihh1_windows.append(window_ihh1)
		all_ihh0_windows.append(window_ihh0)
		all_freq1_windows.append(window_freq1)
		iSnp += jumplen
		#print(str(len(all_score_windows)))
	numWindows = len(all_score_windows)
	print("chunked " + str(numSnps) + " SNPs with available haplotype scores into " + str(numWindows) + " windows of size " + str(windowsize) + " SNPs with jump length " + str(jumplen) + ".")

	################################
	## CALC AVE AND WRITE TO FILE  #
	################################

	writefile = open(writefilename, 'w')
	for iWindow in range(numWindows):
		positions = all_phys_windows[iWindow]
		startPos, endPos = positions[0], positions[1]
		scores = all_score_windows[iWindow]
		winiHS = calc_winIhs(scores)

		ihh1_scores = all_ihh1_windows[iWindow] #since these are all positive,
		ihh0_scores = all_ihh0_windows[iWindow] #we can just reuse calc_winIhs func
		freq1_vals = all_freq1_windows[iWindow] #to get averages.

		winihh1 = calc_winIhs(ihh1_scores)
		winihh0 = calc_winIhs(ihh0_scores)
		freq_ave = calc_winIhs(freq1_vals)
		

		#writeline = str(startPos) + "\t" + str(endPos) + "\t" + str(winiHS) +'\n'
		writeline = str(startPos) + "\t" + str(endPos) + "\t" + str(winiHS) +'\t' + str(winihh1) + "\t" + str(winihh0) + "\t" + str(freq_ave) + "\n"
		writefile.write(writeline)
	writefile.close()
	print("wrote to file: " + writefilename)

	return
def execute_interpolate_hapscores(args):
	'''from interpolate.py'''
	inputTpedFilename = args.intpedfilename
	inputIhsFilename = args.inihsfilename
	inputWinihsFilename = args.inwinihsfilename
	outputFilename = args.outfilename
	all_snps = []
	openfile = open(inputTpedFilename, 'r')
	for line in openfile:
		entries = line.split()
		pos = int(entries[1])
		all_snps.append(pos)
	openfile.close()
	nsnps = len(all_snps)
	print('loaded ' + str(nsnps) + ' snps from' + inputTpedFilename)

	# load scores and windows from winIhs 
	openfile = open(inputWinihsFilename, 'r')
	starts, ends, scores = [], [], []
	ihh1s, ihh0s = [], []
	freqs = []

	for line in openfile:
		entries = line.split()
		startBp, endBp, winIhs = int(entries[0]), int(entries[1]), float(entries[2])
		
		ihh1, ihh0 = float(entries[3]), float(entries[4])
		freq_win = float(entries[5])

		starts.append(startBp)
		ends.append(endBp)
		scores.append(winIhs)
		ihh1s.append(ihh1)
		ihh0s.append(ihh0)
		freqs.append(freq_win)
	openfile.close()

	readfile = open(inputIhsFilename, 'r')
	writefile = open(outputFilename, 'w')
	readline = readfile.readline()
	entries = readline.split()
	thisScorePos = int(entries[0])

	for iSnp in range(nsnps):
		thisPos = all_snps[iSnp]
		if thisPos < thisScorePos: #interpolate until we advance to the first one for which we have a score
			interpolated_ihs, interpolated_ihh1, interpolated_ihh0, interpol_freq = interpolate(starts, ends, scores, ihh1s, ihh0s, thisPos, freqs)

			unnormed_ihs= 0 #dummy to pass to selscan so norm will function
			writeline = str(thisPos) + "\t" + str(thisPos) + "\t" + str(interpol_freq) + "\t" + str(interpolated_ihh1) + "\t" + str(interpolated_ihh0) + "\t" + str(unnormed_ihs) + "\t" + str(interpolated_ihs) + "\t9\n" #index interpolated in outputfile, in case it's not already obvious?
			writefile.write(writeline)
		elif thisPos == thisScorePos: #found a match
			writefile.write(readline) #propagate to writefile
			readline = readfile.readline() #advance to next SNP for which we have a full iHS calc
			entries = readline.split()
			if len(entries) < 1:
				break
			else:
				thisScorePos = int(entries[0])
		else: #thisPos > thisScorePos
			print("ERROR: does your input iHS score file contain SNPs missing from your input TPED? Shame.")

	readfile.close()
	writefile.close()
	print('wrote to ' + outputFilename)


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
	print cmdstring
	#subproces.check_output(cmdstring.split())
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
		argstring = str(args.startBp) + " " + str(args.endBp) + " " + args.outfile + " " + delihh_hit_filename + " " + delihh_miss_filename + " " + ihs_hit_filename + " " + ihs_miss_filename + " " + xpehh_hit_filename + " " + xpehh_miss_filename + " " + fst_hit_filename + " " + fst_miss_filename + " " + deldaf_hit_filename + " " + deldaf_miss_filename 
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
