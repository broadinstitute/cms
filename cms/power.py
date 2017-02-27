##	top-level script to manipulate and analyze empirical/simulated CMS output
##	last updated 02.27.2017	vitti@broadinstitute.org

import matplotlib as mp 
mp.use('agg')
import matplotlib.pyplot as plt
from power.power_func import merge_windows, get_window, check_outliers, check_rep_windows, calc_pr, get_pval, plotManhattan, \
						plotManhattan_extended, quick_plot, get_causal_rank, get_cdf_from_causal_ranks, plot_dist
from power.parse_func import get_neut_repfile_name, get_sel_repfile_name, get_emp_cms_file, read_cms_repfile, \
						read_pr, read_vals_lastcol, get_pr_filesnames, load_regions
from tempfile import TemporaryFile
#from xlwt import Workbook, easyxf #add to cms-venv 
#from pybedtools import BedTool 
import numpy as np
import argparse
import sys
import os

####################
## DEFINE PARSER ###
####################
def full_parser_power():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for calculating CMS 2.0 power from simulated data and significance for CMS scores from empirical data.")
	subparsers = parser.add_subparsers(help="sub-commands")
	
	#######################
	## VISUALIZE OUTPUT ###
	#######################
	repviz_parser = subparsers.add_parser('repviz', help="visualize component and combined scores")
	repviz_parser.add_argument('--vizRep', type=int, help="rep number", default=1, action="store")
	distviz_parser = subparsers.add_parser('distviz', help="visualize distribution of CMS scores")
	distviz_parser.add_argument('--oneFile', action="store", default=None, help="quick viz score dist from a single file")
	distviz_parser.add_argument('--simsDist', action="store_true", default=False, help="visualize from all simulated replicates")
	distviz_parser.add_argument('--empDist', action="store_true", default=False, help="visualize from empirical scores")
	distviz_parser.add_argument('--normed_cms', action="store_true", help="normalized values", default = False)

	#####################
	## QUANTIFY POWER ###
	#####################
	cdf_parser = subparsers.add_parser('cdf', help = 'plot cumulative density function of causal rank')
	fpr_parser = subparsers.add_parser('fpr', help='calculate false positive rate for CMS_gw based on neutral simulations')
	tpr_parser = subparsers.add_parser('tpr', help='calculate false positive rate for CMS_gw based on simulations with selection')
	roc_parser = subparsers.add_parser('roc', help="calculate receiving operator characteristic curve given false and true positive rates")
	roc_parser.add_argument('--plot_curve', action="store_true", default=False)
	roc_parser.add_argument('--find_opt', action="store_true", default=False)
	roc_parser.add_argument('--maxFPR', type=float, action="store", default=.001)
	cdf_parser.add_argument('--selPos', type=int, action='store', default=750000, help="position of the causal allele in simulates")
	find_cutoff_parser = subparsers.add_parser('find_cutoff', help="get best TPR for a given FPR and return threshhold cutoffs for region detection")
	find_cutoff_parser.add_argument('--maxFPR', type=float, action="store", default=".05")
	find_cutoff_parser.add_argument('fprloc', type=str, action="store", help="specific to model/pop")
	find_cutoff_parser.add_argument('tprloc', type=str, action="store", help="specific to model/pop")

	#############################
	## EMPIRICAL SIGNIFICANCE ###
	#############################
	gw_regions_parser = subparsers.add_parser('gw_regions', help="pull designated significant regions from genome-wide normalized results")
	gw_regions_parser.add_argument('--geneFile')
	regionlog_parser = subparsers.add_parser('regionlog', help='write regions to excel sheet with gene overlap')
	regionlog_parser.add_argument('--gene_bedfile', help="name of file", type = str, action='store', default = "/idi/sabeti-scratch/jvitti/knownGenes_110116.txt")
	manhattan_parser = subparsers.add_parser('manhattan', help='generate manhattan plot of p-values of empirical results.')	
	manhattan_parser.add_argument('--zscores', action = 'store_true', help="plot -log10(p-values) estimated from neutral simulation") #nix
	manhattan_parser.add_argument('--maxSkipVal', help="expedite plotting by ignoring anything obviously insignificant", default=-10e10)
	manhattan_parser.add_argument('--poolModels', help="experimental - combine neut distributions from multiple demographic models")
	extended_manhattan_parser = subparsers.add_parser('extended_manhattan', help = "generate per-chrom plots as one fig")
	extended_manhattan_parser.add_argument('--plotscore', help="string label for score to plot: {seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed}", type=str, default="cms_normed")
	extended_manhattan_parser.add_argument('--regionsfile', help="optional; input file of regions designated as above threshhold")
	extended_manhattan_parser.add_argument('--percentile', help="percentile to hilite")
	extended_manhattan_parser.add_argument('--titlestring', help="title for plot")
	extended_manhattan_parser.add_argument('--dpi', help="resolution for matplotlib", type=int, default=100)

	##################
	## SHARED ARGS ###
	##################

	for write_parser in [fpr_parser, tpr_parser, roc_parser, repviz_parser, cdf_parser]:
		write_parser.add_argument('--writedir', type =str, help='where to write output', default = "/idi/sabeti-scratch/jvitti/")
		write_parser.add_argument('--checkOverwrite', action="store_true", default=False)
		write_parser.add_argument('--simpop', action='store', help='simulated population', default=1)

	for regions_parser in [fpr_parser, gw_regions_parser, regionlog_parser, tpr_parser]:
		regions_parser.add_argument('regionlen', type = int, action='store', help='length of region to query', default="100000")
		regions_parser.add_argument('thresshold', type = float, action='store', help='percentage of region to exceed cutoff', default="30")
		regions_parser.add_argument('cutoff', type = float, action='store', help='minimum significant value for region definition', default="3.0")
		regions_parser.add_argument('--saveLog', type =str, help="save results as text file", )

	for emp_parser in [manhattan_parser, extended_manhattan_parser, gw_regions_parser, distviz_parser]:
		emp_parser.add_argument('--emppop', action='store', help='empirical population', default="YRI")

	for suffixed_parser in [fpr_parser, tpr_parser, cdf_parser, manhattan_parser, extended_manhattan_parser, gw_regions_parser, distviz_parser]:
		suffixed_parser.add_argument('--suffix', type= str, action='store', default='')

	for plot_parser in [repviz_parser, distviz_parser, manhattan_parser, extended_manhattan_parser, cdf_parser]:
		plot_parser.add_argument('--savefilename', action='store', help='path of file to save', default="/web/personal/vitti/test.png")

	for commonparser in [repviz_parser, distviz_parser, fpr_parser, gw_regions_parser, manhattan_parser,
						regionlog_parser, cdf_parser, tpr_parser, extended_manhattan_parser]:
		commonparser.add_argument('--model', type=str, default="nulldefault")
		commonparser.add_argument('--nrep', type=int, default=1000)

	return parser

#################################
## DEFINE EXECUTIVE FUNCTIONS ###
#################################
########	Visualize composite score 
########	output for a given CMS run.
def execute_repviz(args):
	'''visualize CMS scores for simulated data as individual replicates (w/ component and combined scores vs. pos)'''
	causal_index = -1
	model, pop = args.model, args.simpop
	irep = args.vizRep
	savefilename = args.savefilename
	writedir = args.writedir
	#CURRENTLY JUST NEUT
	scenars = ['neut']#'neut', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
	for scenar in scenars:
		if scenar == "neut":
			cmsfilename = get_neut_repfile_name(model, irep, pop, normed = True, basedir=writedir)
		else:
			cmsfilename = get_sel_repfile_name(model, irep, pop, scenar, normed =True, vsNeut = True, basedir=writedir)
			#print(cmsfilename)
		if os.path.isfile(cmsfilename):
			print('loading from... ' + cmsfilename)
			physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(cmsfilename)

			f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8, sharex = True)
			quick_plot(ax1, physpos, ihs_normed, "ihs_normed", causal_index)
			quick_plot(ax2, physpos, delihh_normed, "delihh_normed", causal_index)
			quick_plot(ax3, physpos, nsl_normed, "nsl_normed", causal_index)
			quick_plot(ax4, physpos, xpehh_normed, "xpehh_normed", causal_index)
			quick_plot(ax5, physpos, fst, "fst", causal_index)
			quick_plot(ax6, physpos, deldaf, "deldaf", causal_index)
			quick_plot(ax7, physpos, cms_unnormed, "cms_unnormed", causal_index)
			quick_plot(ax8, physpos, cms_normed, "cms_normed", causal_index)				
			savefilename = args.savefilename + "_"+ model + "_" + str(pop) + "_" + str(scenar) + "_" + str(irep) + ".png"

			plt.savefig(savefilename)
			#if irep == 1:
			print("plotted to " + savefilename)
			plt.close()
	return
def execute_distviz(args):
	"""from explore_p3_dists.py"""
	reps = args.nrep

	############################
	## VIEW DIST FROM 1 FILE  ##
	############################
	if args.oneFile is not None:
		savefilename = args.savefilename
		infilename = args.oneFile
		if os.path.isfile(infilename):
			values = read_vals_lastcol(infilename)
			plot_dist(values, savefilename)

	#########################
	## VIEW DIST FROM SIMS ##
	#########################
	if args.simsDist:
		#### SHOULD USE DISPATCH TO PASS MODEL AND MODELPPOP
		#	argstring = "--simpop " + str(simpop) + " --model " + model
		#models = [ 'gradient_101915_treebase_6_best', 'default_default_101715_12pm', 'default_112115_825am','nulldefault_constantsize', 'nulldefault', ]
		#modelpops = [1, 2, 3, 4]
		models = [args.model]
		modelpops = [args.simpop]

		for model in models:		
			for pop in modelpops:
				allvals = []
				#NEUT REPS
				for irep in range(1, reps+1):
					filebase = "/idi/sabeti-scratch/jvitti/clean/scores/"
					infilename = filebase + model + "/neut/composite/rep" + str(irep) + "_" + str(pop) + ".cms.out"
					if os.path.isfile(infilename):
						values = read_vals_lastcol(infilename)
						allvals.extend(values)
				plot_dist(allvals, "/web/personal/vitti/sim_dists/" + model +"_" + str(pop) + "_justneut.png")

				#SEL REPS
				#for selbin in selbins:	
				#	for irep in range(1,501):
				#		infilename = filebase + model + "/sel" + str(pop) + "/sel_" + str(selbin) + "/rep" + str(irep) + "_" + str(pop) + "_vsneut.cms.out"
				#		if os.path.isfile(infilename):
				#			values = read_vals_lastcol(infilename)
				#			allvals.extend(values)
			
				#plot_dist(allvals, savefilename)

	#############################
	## VIEW DIST FROM EMP DATA ##
	#############################
	elif args.empDist:
		for model in models:
			for pop in emppops:
				savefilename = "/web/personal/vitti/" + model +"_" + str(pop) + ".png"
				allvals = []
				for chrom in chroms:
					infilename = filebase + pop + "_" + model + "_" + "neut" + ".chr" + str(chrom) + ".txt"
					if os.path.isfile(infilename):
						values = read_vals_lastcol(infilename)
						allvals.extend(values)
				plot_dist(allvals, savefilename)
	return
def execute_manhattan(args):
	selpop = args.emppop
	model = args.model
	savename = args.savefilename
	suffix = args.suffix
	#nRep = args.nrep
	###############################
	### LOAD NEUTRAL SIM VALUES ###
	###############################
	modelpops = {'YRI':1, 'CEU':2, 'CHB':3, 'BEB':4}
	pop = modelpops[selpop]

	#EXPERIMENTAL
	#This is not how Shari actually normalized for figure 1 2013.
	#if args.poolModels:
	#	all_neut_rep_scores = []
	#	for model in ['default_112115_825am', 'nulldefault', 'nulldefault_constantsize', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm']:
	#		neut_rep_scores = load_normed_simscores(model, pop, vsNeut=vsNeut, numRep=nRep)
	#		all_neut_rep_scores.extend(neut_rep_scores)
	#	neut_rep_scores = all_neut_rep_scores
	#else:
	#	neut_rep_scores = load_normed_simscores(model, pop, vsNeut=vsNeut, numRep=nRep)
	#print('loaded ' + str(len(neut_rep_scores)) + ' neut simscores...') 
	neut_rep_scores = []

	#############################
	### LOAD EMPIRICAL VALUES ###
	#############################
	all_emp_pos, all_emp_scores = [], []
	nSnps = 0
	for chrom in range(1,23):
		thesepos, thesescores = [], []
		emp_cms_filename = get_emp_cms_file(selpop, model, chrom, normed=True, suffix=suffix)
		print('loading chr ' + str(chrom) + ": " + emp_cms_filename)
		if not os.path.isfile(emp_cms_filename):
			print("missing: " + emp_cms_filename)
			break
		physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(cmsfilename)
		#physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(emp_cms_filename)
		all_emp_pos.append(physpos)
		all_emp_scores.append(cms_normed)

	###########################
	### DRAW MANHATTAN PLOT ###
	###########################
	f, ax = plt.subplots(1)
	if args.zscores:
		calc_zscores = True
	else:
		calc_zscores = False

	plotManhattan(ax, neut_rep_scores, all_emp_scores, all_emp_pos, nSnps, zscores=calc_zscores, maxSkipVal = args.maxSkipVal)

	plt.savefig(savename)
	print('saved to: ' + savename)
	return
def execute_extended_manhattan(args):
	plotscore = args.plotscore
	selpop = args.emppop
	model = args.model
	suffix = args.suffix
	savename = args.savefilename
	dpi = args.dpi
	numChr = 22
	titlestring = args.titlestring

	modelpops = {'YRI':1, 'CEU':2, 'CHB':3, 'BEB':4}
	pop = modelpops[selpop]
	#colorDict = {1:'#FFB933', 2:'#0EBFF0', 3:'#ADCD00', 4:'#8B08B0'} #1000 Genomes group color scheme
	colorDict = {1:'#cec627', 2:'#0EBFF0', 3:'#65ff00', 4:'#8B08B0'} #make it pop-!
	

	f, axarr = plt.subplots(numChr, 1, sharex = True, sharey=True, dpi=dpi, figsize=(7, 10))
	plt.suptitle(titlestring, fontsize=10)

	plt.xlabel('position')
	plt.ylabel('cms_gw normed score')

	all_emp_pos, all_emp_scores = [], []
	for chrom in range(1,numChr +1):
		emp_cms_filename = get_emp_cms_file(selpop, model, chrom, normed=True, suffix=suffix)
		print('loading chr ' + str(chrom) + ": " + emp_cms_filename)
		if not os.path.isfile(emp_cms_filename):
			print("missing: " + emp_cms_filename)
			break
		physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(emp_cms_filename)

		iax = chrom-1
		ax = axarr[iax]
		#ax.grid()
		plot_data = eval(plotscore)
		plotManhattan_extended(ax, plot_data, physpos, chrom)
		all_emp_pos.append(physpos)
		all_emp_scores.append(plot_data)

	################################
	## HILITE SIGNIFICANT REGIONS ##
	################################

	if args.regionsfile is not None:
		regionchrs, regionstarts, regionends = load_regions(args.regionsfile)
		print('loaded ' + str(len(regionchrs)) + ' significant regions from ' + args.regionsfile)
		for iregion in range(len(regionchrs)):
			regionchr, regionstart, regionend = regionchrs[iregion], regionstarts[iregion], regionends[iregion]
			this_chrom = int(regionchr.strip('chr'))
			ichrom = this_chrom-1
			chrompos, chromscores = all_emp_pos[ichrom], all_emp_scores[ichrom]
			zipped = zip(chrompos, chromscores)
			plotpos, plotvals = [], []
			for locus in zipped:
				if locus[0] >= regionstart:
					plotpos.append(locus[0])
					plotvals.append(locus[1])
				if locus[0] > regionend:
					break
			axarr[ichrom].plot(plotpos, plotvals, color=colorDict[pop], markersize=1)

	if args.percentile is not None:
		percentile = float(args.percentile)
		print('plotting data with heuristic cutoff for ' + str(percentile) + " percentile...")
		flat_emp_scores = [item for sublist in all_emp_scores for item in sublist if not np.isnan(item)]
		score_cutoff = float(np.percentile(flat_emp_scores, percentile))
		print("score cutoff: " + str(score_cutoff))
		for chrom in range(1,numChr +1):
			iax = chrom-1
			ax = axarr[iax]
			maximumVal = ax.get_xlim()[1]
			xpoints = np.array([0, maximumVal])
			ypoints = np.array([score_cutoff, score_cutoff])
			ax.plot(xpoints, ypoints ,linestyle = "dotted", color="red", markersize=.3)

			#get empirical scores and positions for pass threshhold and plot them as above with color
			these_scores, these_pos = all_emp_scores[iax], all_emp_pos[iax]
			zipped =  zip(these_scores, these_pos)
			significant = [item for item in zipped if item[0] >= score_cutoff]
			signif_vals = [item[0] for item in significant]
			signif_pos = [item[1] for item in significant]
			ax.plot(signif_pos, signif_vals, color=colorDict[pop], linestyle='None', marker=".", markersize=.3)#, markersize=1)

	plt.savefig(savename)
	print('saved to: ' + savename)
	return

########	Quantify and visualize power
########	across significance cutoffs.
def execute_cdf(args):
	reps = args.nrep
	savefilename = args.savefilename
	writedir = args.writedir
	scenars = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
	model = args.model
	causalPos = args.selPos
	suffix = args.suffix
	#causal_ranks_all = []
	causal_ranks_1, causal_ranks_2, causal_ranks_3, causal_ranks_4 = [], [], [], []
	for pop in [1, 2, 3, 4]:
		for scenar in scenars:
			for irep in range(1, reps+1):
				cmsfilename = get_sel_repfile_name(model, irep, pop, scenar, normed =True, basedir=writedir, suffix=suffix)
			
				if os.path.isfile(cmsfilename):
					physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(cmsfilename)
					if causalPos in physpos:
						causal_index = physpos.index(causalPos)
						causal_unnormed = cms_unnormed[causal_index]
						causal_rank = get_causal_rank(cms_unnormed, causal_unnormed)
						#print('causal rank: ' + str(causal_rank)) 
						#causal_ranks.append(causal_rank)
						this_array = eval('causal_ranks_' + str(pop))
						this_array.append(causal_rank)
				#else:
				#	print("missing; " + cmsfilename)
	print("for pop 1, loaded " + str(len(causal_ranks_1)) + " replicates.")
	print("for pop 2, loaded " + str(len(causal_ranks_2)) + " replicates.")
	print("for pop 3, loaded " + str(len(causal_ranks_3)) + " replicates.")
	print("for pop 4, loaded " + str(len(causal_ranks_4)) + " replicates.")

	cdf_fig, cdf_ax = plt.subplots()
	if len(causal_ranks_1) > 0:
		cdf_bins1, cdf1 = get_cdf_from_causal_ranks(causal_ranks_1)
		cdf_ax.plot(cdf_bins1[1:], cdf1, color="yellow")
	if len(causal_ranks_2) > 0:
		cdf_bins2, cdf2 = get_cdf_from_causal_ranks(causal_ranks_2)
		cdf_ax.plot(cdf_bins2[1:], cdf2, color="blue")
	if len(causal_ranks_3) > 0:
		cdf_bins3, cdf3 = get_cdf_from_causal_ranks(causal_ranks_3)
		cdf_ax.plot(cdf_bins3[1:], cdf3, color="green")
	if len(causal_ranks_4) > 0:
		cdf_bins4, cdf4 = get_cdf_from_causal_ranks(causal_ranks_4)			
		cdf_ax.plot(cdf_bins4[1:], cdf4, color="purple")
	cdf_ax.set_xlim([0, 50])
	plt.title(model) #+ ", " + str(len(causal_ranks)) + " selection replicates")
	plt.ylabel('probability that the causal variant is captured')
	plt.xlabel('significance thresshold (i.e., examining the top x variants)')
	plt.savefig(savefilename)
	plt.close()
	print('plotted to ' + savefilename)
	return
def execute_fpr(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	numReps = args.nrep
	pop = args.simpop
	suffix = args.suffix
	writedir = args.writedir

	all_scores = []
	all_percentages = []
	
	if True:
		for irep in range(1, numReps + 1):
			repfilename = get_neut_repfile_name(model, irep, pop, normed=True, suffix=suffix, basedir=writedir)
			if (irep==1):
				print(repfilename)
			physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
			#physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
			if len(cms_normed) > 0:
				all_scores.append(cms_normed)
				rep_percentages = check_rep_windows(physpos, cms_normed, regionlen, cutoff = cutoff)
				all_percentages.append(rep_percentages)		
				#FOR DEBUG
				#print(str(rep_percentages) + "\t" + repfilename)
				if len(rep_percentages) > 0:
					if max(rep_percentages) > thresshold:
						print("false positive: " + repfilename)

	print('loaded ' + str(len(all_scores)) + " replicates populations for model " + model + "...")
	fpr = calc_pr(all_percentages, thresshold)
	print('false positive rate: ' + str(fpr) + "\n")

	if args.saveLog	is not None:
		writefilename = args.saveLog 
		writefile = open(writefilename, 'w')
		writefile.write(str(fpr)+'\n')

		writefile.write(model + "\t" + str(regionlen) + "\t" + str(thresshold) + '\t' + str(cutoff) + '\n')
		writefile.close()
		print('wrote to :  ' + str(writefilename))
	return
def execute_tpr(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	numReps = args.nrep
	pop = args.simpop
	suffix = args.suffix
	writedir = args.writedir

	all_scores = []
	all_percentages = []
	
	#if args.saveLog	is not None:
	#	writefilename = args.saveLog
	#	if os.path.isfile(writefilename):
	#		print(writefilename + " already exists; aborting.")
	#		sys.exit(0)

	#per seldaf
	dafbins = [['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90'], ['0.10', '0.20', '0.30'], ['0.40', '0.50', '0.60'], ['0.70', '0.80', '0.90'], ['0.90']]
	daflabels = ['all', 'lo', 'mid', 'hi','highest']
	for ibin in [1, 2, 3, 4]:#range(1):
		thesebins, thislabel = dafbins[ibin], daflabels[ibin]
		allrepfilenames = []
		for selbin in thesebins:
			for irep in range(1, numReps + 1):
				repfilename = get_sel_repfile_name(model, irep, pop, selbin, normed=True, suffix=suffix, basedir=writedir)
				if (irep==1):
					print(repfilename)
				if os.path.isfile(repfilename):
					allrepfilenames.append(repfilename)
		print('loaded ' + str(len(allrepfilenames)) + " replicates...")
		#numToTake = min(500, len(allrepfilenames))
		#chosen = np.random.choice(allrepfilenames, numToTake, replace=False) #take random sample	
		chosen = allrepfilenames #this was just to expedite, no?
		for repfilename in chosen:
			physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
			#physpos, genpos, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(repfilename)
			if len(cms_normed) > 0:
				all_scores.append(cms_normed)
				rep_percentages = check_rep_windows(physpos, cms_normed, regionlen, cutoff = cutoff)
				all_percentages.append(rep_percentages)		

		print('loaded ' + str(len(all_scores)) + " replicates populations for model " + model + "...")
		tpr = calc_pr(all_percentages, thresshold)
		print('true positive rate: ' + str(tpr) + "\n")

		if args.saveLog	is not None:
			writefilename = args.saveLog +"_" + thislabel
			writefile = open(writefilename, 'w')
			writefile.write(str(tpr)+'\n')

			writefile.write(model + "\t" + str(regionlen) + "\t" + str(thresshold) + '\t' + str(cutoff) + '\n')
			writefile.close()
			print('wrote to :  ' + str(writefilename))
	return	
def execute_roc(args):
	''' from quick_roc.py '''
	plot_roc = args.plot_curve
	find_opt = args.find_opt
	maxFPR = args.maxFPR
	writedir = args.writedir

	regionlens = [10000, 25000, 50000, 100000]
	thressholds = [25, 30, 35, 40, 45, 50]
	cutoffs = [2, 2.5, 3, 3.5, 4, 4.5, 5]
	models = ['default_112115_825am', 'nulldefault', 'nulldefault_constantsize', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm']#
	pops = [1, 2, 3, 4]

	###############
	## LOAD DATA ##
	###############
	allfpr = {}
	alltpr = {}
	for regionlen in regionlens:
		for percentage in thressholds:
			for cutoff in cutoffs:
				for model in models:
					allpops_fpr, allpops_tpr = [], []
					for pop in pops:
						this_key = (regionlen, percentage, cutoff, model, pop)
						fprfile, tprfile = get_pr_filesnames(this_key, writedir)
						if os.path.isfile(fprfile) and os.path.isfile(tprfile) and os.path.getsize(fprfile) > 0 and os.path.getsize(tprfile) > 0:
							fpr = read_pr(fprfile)
							tpr = read_pr(tprfile)
							allfpr[this_key] = fpr
							alltpr[this_key] = tpr
							allpops_fpr.append(fpr)
							allpops_tpr.append(tpr)
						else:
							print("missing " + tprfile + "\t" + fprfile)
					#assert len(allpops_fpr) == 4
					ave_fpr = np.average(allpops_fpr)
					#assert len(allpops_tpr) == 4
					ave_tpr = np.average(allpops_tpr)

					popave_key = (regionlen, percentage, cutoff, model, "ave")
					allfpr[popave_key] = ave_fpr
					alltpr[popave_key] = ave_tpr

	###############
	## PLOT DATA ##
	###############
	if plot_roc:
		fig, ax = plt.subplots(1)
		for model in models:
			#need to plotfpr from dict
			plotfpr, plottpr = [], []
			for regionlen in regionlens:
				for percentage in thressholds:
					for cutoff in cutoffs:
						this_key = (regionlen, percentage, cutoff, model, "ave")
						plotfpr.append(allfpr[this_key])
						plottpr.append(alltpr[this_key])

			if (len(allfpr)) > 0:
				#ax.scatter(plotfpr, plottpr, label=model)
				plotfpr, plottpr = zip(*sorted(zip(plotfpr, plottpr)))
				ax.plot(plotfpr, plottpr, label=model)

		ax.set_xlabel('FPR')
		ax.set_ylabel('TPR')
		ax.set_xlim([0,.1])
		ax.set_ylim([0,1])
		savefilename = "/web/personal/vitti/roc/" +"compare_maf.png"#+ model + ".png"
		plt.legend()
		plt.savefig(savefilename)
		print("plotted to " + savefilename)
		plt.close()
			
	#####################
	## COMPARE CUTOFFS ##
	#####################
	if find_opt:
		for model in models:
			print(model)
			for regionlen in regionlens:
				for percentage in thressholds:
					for cutoff in cutoffs:
						avekey = (regionlen, percentage, cutoff, model, "ave")
						#FPR passes thresshold. What's the best TPR we can get?
						if allfpr[avekey] <= maxFPR:
							print(avekey)
							print(str(alltpr[avekey]))

	

	return
def execute_find_cutoff(args):
	''' a little cleaner and simpler '''

	maxFPR = args.maxFPR
	fpr_loc = args.fprloc
	tpr_loc = args.tprloc 
	################################
	## LOAD WHATEVER DATA WE HAVE ##
	################################
	all_fpr, all_tpr = {}, {}
	freq_classes = ['lo', 'mid', 'hi', 'highest']
	fpr_files = os.listdir(fpr_loc)
	tpr_files = os.listdir(tpr_loc)
	for fprfile in fpr_files:
		if "alt" not in fprfile:
			fpr = read_pr(fpr_loc + fprfile)
			entries = fprfile.split("_")
			regionlen, thresshold, cutoff = int(entries[1]), int(entries[2]), int(entries[3])
			fprkey = (regionlen, thresshold, cutoff)
			all_fpr[fprkey] = fpr
			for freq_class in freq_classes:
				tprfile = tpr_loc + "tpr_" + str(regionlen) + "_" + str(thresshold) + "_" + str(cutoff) + "_" + str(freq_class)
				if os.path.isfile(tprfile):			
					tpr = read_pr(tprfile)
					tprkey = (regionlen, thresshold, cutoff, freq_class)
					all_tpr[tprkey] = tpr
	#############################
	## CHOOSE OPT MEETING CRIT ##
	#############################
	best_tpr, best_fpr = 0, 0
	best_cutoff = 0
	fpr_keys = all_fpr.keys()
	tpr_keys = all_tpr.keys()
	for freq_class in freq_classes:
		print("Now finding optimal thresshold for " + freq_class + " with a maximum FPR of " + str(maxFPR))
		for key in fpr_keys:
			if all_fpr[key] <= maxFPR:
				tprkey = (key[0], key[1], key[2], freq_class)
				if tprkey in tpr_keys:
					tpr = all_tpr[tprkey]
					if tpr > best_tpr:
						best_tpr = tpr
						best_cutoff = tprkey
						best_fpr = all_fpr[key]
		print(best_cutoff)
		print("FPR: " + str(best_fpr))
		print("TPR: " + str(best_tpr))
	return

########	Apply significance cutoffs
########	to empirical results.
def execute_gw_regions(args):
	model = args.model
	regionlen = args.regionlen
	thresshold = args.thresshold
	cutoff = args.cutoff
	pop = args.emppop
	windowlen = args.regionlen
	suffix = args.suffix

	chroms = range(1,23)
	signif_windows = []
	####################
	## LOOP OVER CHRS ##
	####################
	for chrom in chroms:
		chrom_signif = []
		normedempfilename = get_emp_cms_file(pop, model, chrom, normed=True, suffix = suffix)
		if not os.path.isfile(normedempfilename):
			print("missing: " + normedempfilename)
		else:
			physpos, genpos, seldaf, ihs_normed, delihh_normed, nsl_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(normedempfilename)
			for iPos in range(len(physpos)):
				##################
				## CHECK REGION ##
				##################
				window_scores = get_window(iPos, physpos, cms_normed, windowlen)
				percentage = check_outliers(window_scores, cutoff)
				if percentage > thresshold:
					chrom_signif.append(physpos[iPos])
		signif_windows.append(chrom_signif)

	##############################
	## MERGE CONTIGUOUS WINDOWS ##
	##############################
	final_starts = []
	final_ends = []
	print('merging regions')
	for chrom_signif in signif_windows:
		starts, ends =  merge_windows(chrom_signif, windowlen)
		final_starts.append(starts)
		final_ends.append(ends)
	#print(final_starts)
	#print(final_ends)

	###################
	## WRITE TO FILE ##
	###################
	if args.saveLog is not None:
		writefilename = args.saveLog 
		writefile = open(writefilename, 'w')
		for ichrom in range(len(final_starts)):
			chromnum = ichrom + 1

			starts = final_starts[ichrom]
			ends = final_ends[ichrom]

			for iregion in range(len(starts)-1):
				writeline = "chr" + str(chromnum) + "\t" + str(starts[iregion]) + "\t" + str(ends[iregion]) + '\n'
				writefile.write(writeline)
		writefile.close()
		print('wrote to ' + writefilename)	
	return
def execute_regionlog(args):
	''' writes an excel file '''
	model = args.model	
	regionfiles = []
	pops = ['YRI', 'CEU', 'CHB', 'BEB']
	takepops = []
	for pop in pops: #too tired to soft-code rn
		regionfilename = "/n/regal/sabeti_lab/jvitti/clear-synth/1kg_regions/" + model + "/" + pop + "_" + str(int(args.regionlen)) + "_" + str(int(args.thresshold)) + "_" + str(int(args.cutoff)) + ".txt" #or _alt
		print(regionfilename)
		if os.path.isfile(regionfilename):
			regionfiles.append(regionfilename)
			takepops.append(pop)
	if len(regionfiles) == 0:
		print("found no regions")
		return
	else:
		totalselregions = 0
		print('loaded regions from ' + str(len(regionfiles)) + " files...")
		boldstyle = easyxf('font: bold 1;')
		wrapstyle = easyxf('alignment: wrap on, vert center, horiz center')
		book = Workbook()
		sheet1 = book.add_sheet('gw significant regions')
		header = ['chrom', 'start', 'end', 'len (kb)', 'genes', 'pop']
		for icol in range(len(header)):
			sheet1.write(0, icol, header[icol], boldstyle)

		#get inclusive overlap
		irow = 0
		for iregionfilename in range(len(regionfiles)):
			regionfilename = regionfiles[iregionfilename]
			pop = takepops[iregionfilename]
			#for regionfilename in regionfiles:
			genes = BedTool(args.gene_bedfile) 
			regions = BedTool(regionfilename)
			intersect = regions.intersect(genes,  wa = True, wb = True)

			#narrow down
			geneDict = {}
			for item in intersect:
				selregion_chr = item[0]
				selregion_start, selregion_end = item[1], item[2]
				key = (selregion_chr, selregion_start, selregion_end)
				generegion_id = item[6]
				if key not in geneDict.keys():
					geneDict[key] = [generegion_id]
				else:
					geneDict[key].append(generegion_id)

			#for irow in range(len(regions)):
			for region in regions:
				totalselregions +=1
				irow +=1
				#region = regions[irow]
				chrom, start, end = region[0], region[1], region[2]
				key = (chrom, start, end)
				regionlen = int(end) - int(start)
				kb_regionlen=round(regionlen/1000)

				if key in geneDict.keys():
					genelist = geneDict[key]
					genes = set(genelist)
					genestring = ""
					for gene in genes:
						genestring += gene + ", "
					genestring = genestring[:-2]
					sheet1.write(irow+1, 0, chrom, wrapstyle)
					sheet1.write(irow+1, 1, int(start), wrapstyle)
					sheet1.write(irow+1, 2, int(end), wrapstyle)
					sheet1.write(irow+1, 3, kb_regionlen, wrapstyle)
					sheet1.write(irow+1, 4, genestring, wrapstyle)
					sheet1.write(irow+1, 5, pop, wrapstyle)
				else:
					sheet1.write(irow+1, 0, chrom, wrapstyle)
					sheet1.write(irow+1, 1, int(start), wrapstyle)
					sheet1.write(irow+1, 2, int(end), wrapstyle)
					sheet1.write(irow+1, 3, kb_regionlen,wrapstyle)			
					sheet1.write(irow+1, 5, pop,wrapstyle)

	colWidths = [10, 15, 15, 10, 25, 10]
	for icol in range(len(colWidths)):
		sheet1.col(icol).width = colWidths[icol] * 256 #~x char wide	

	book.save(args.saveLog)
	book.save(TemporaryFile())
	print('wrote ' + str(totalselregions) + ' significant regions to: ' + args.saveLog)
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_power()
	args = runparser.parse_args()

	# if called with no arguments, print help
	if len(sys.argv)==1:
		runparser.parse_args(['--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline