#!/usr/bin/env python
## top-level script for combining scores into composite statistics as part of CMS 2.0.
## last updated: 06.19.2017 	vitti@broadinstitute.org #update docstrings

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from power.parse_func import get_neut_repfile_name, get_sel_repfile_name, get_emp_cms_file, get_sim_component_score_files, get_emp_component_score_files, load_empscores
from combine.input_func import write_perpop_ihh_from_xp, write_run_paramfile, write_pair_sourcefile, normalize
from combine.viz_func import hapSort_coreallele, hapSort, hapViz, readAnnotations, find_snp_index, pullRegion, load_from_hap
from dists.scores_func import calc_fst_deldaf, calc_delihh
from dists.freqbins_func import check_create_dir, execute 
from dists.likes_func import get_master_likefiles
import numpy as np
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

	##############################
	## RECALCULATE INPUT SCORES ##
	##############################
	freqscores_parser = subparsers.add_parser('freqscores', help="Calculate allele frequency-based scores (Fst and delDAF) for a pair of populations.") 
	if True:
		freqscores_parser.add_argument('inTped1', type=str, action="store", help="input tped 1")			
		freqscores_parser.add_argument('inTped2', type=str, action="store", help="input tped 2")
		freqscores_parser.add_argument('recomFile', type=str, action="store", help="input recombination file") #should work around this		 
		freqscores_parser.add_argument('outfile', type=str, action="store", help="file to write")
		freqscores_parser.add_argument('--modelpath', action='store', type=str, default='cms/cms/model/', help="path to model directory containing executables") #will become redundant with conda
	delihh_from_ihs_parser = subparsers.add_parser('delihh_from_ihs', help="Calculate delIHH values from iHS output files.")
	if True:
		delihh_from_ihs_parser.add_argument('readfile', type=str, action='store', help='input ihs file')
		delihh_from_ihs_parser.add_argument('writefile', type=str, action='store', help='delihh file to write')
	ihh_from_xp_parser = subparsers.add_parser('ihh_from_xp', help="extract per-pop iHH values from XP-EHH and write to individual files to facilitate arbitrary population comparisons ")
	if True:
		ihh_from_xp_parser.add_argument('inXpehh', type=str, help="input xpehh file")
		ihh_from_xp_parser.add_argument('outIhh', type=str, help="write to file")
		ihh_from_xp_parser.add_argument('takePop', default=1, type=int, help="write for first (1) or second (2) pop in XP-EHH file?")
	xp_from_ihh_parser = subparsers.add_parser('xp_from_ihh', help="Calculate XP-EHH based on two per-pop iHH files.")
	if True:
		xp_from_ihh_parser.add_argument('inIhh1', type=str, action='store', help="input ihh file 1")
		xp_from_ihh_parser.add_argument('inIhh2', type=str, action='store', help="input ihh file 2")
		xp_from_ihh_parser.add_argument('outfilename', type=str, action='store', help="write to file")

	######################################
	## CALCULATING COMPOSITE STATISTICS ##
	######################################
	composite_sims_parser = subparsers.add_parser('composite_sims', help='calculate composite scores for simulations')
	composite_sims_parser.add_argument('--regional_cms', action="store_true", default=False, help="calculate within-region CMS rather than genome-wide CMS")
	composite_emp_parser = subparsers.add_parser('composite_emp', help="calculate composite scores for empirical data")	
	composite_emp_parser.add_argument('--score_basedir', default="/n/regal/sabeti_lab/jvitti/clear-synth/1kg_scores/")
	composite_emp_parser.add_argument('--regional_cms_chrom', type=int, action="store", help="if included, calculate within-region CMS (rather than CMS_gw) for specified bounds at this chromosome")
	
	normsims_genomewide_parser = subparsers.add_parser('normsims_genomewide', help="normalize simulated composite scores to neutral")
	normemp_genomewide_parser = subparsers.add_parser('normemp_genomewide', help="normalize CMS scores to genome-wide") #norm emp REGIONS?
	normemp_genomewide_parser.add_argument('--score_basedir', default="/n/regal/sabeti_lab/jvitti/clear-synth/1kg_scores/")


	###############
	## VISUALIZE ##
	###############
	hapviz_parser = subparsers.add_parser('hapviz', help="Visualize haplotypes for region")
	if True:
		hapviz_parser.add_argument('inputfile', type=str, action="store", help="input tped")
		hapviz_parser.add_argument('--startpos',type=int, help="define physical bounds of region")
		hapviz_parser.add_argument('--endpos',type=int, help="define physical bounds of region")
		hapviz_parser.add_argument('out',type=str,default=None, help="save image as file")
		hapviz_parser.add_argument('--corepos',type=int,default=-1, help="partition haplotypes based on allele status at this position")
		hapviz_parser.add_argument('--title', type=str, default=None, help="title to give to plot")
		hapviz_parser.add_argument('--annotate', type=str, default=None, help="tab-delimited file where each line gives <chr.pos>\t<annotation>")			
		hapviz_parser.add_argument('--maf', type=str, default=None, help="filter on minor allele frequency (e.g. .01, .05)")			
		hapviz_parser.add_argument('--dpi', type=str, default=300, help="image resolution")			
	ml_region_parser = subparsers.add_parser('ml_region', help='machine learning algorithm (within-region)') #connect to AJ work- this could also go in above section
	ucsc_viz_parser = subparsers.add_parser('ucsc_viz', help="Generate trackfiles of CMS scores for visualization in the UCSC genome browser.")
	if True:
		ucsc_viz_parser.add_argument('infile_prefix', type=str, action="store", help="prefix of file containing scores to be reformatted (e.g. 'score_chr' for files named scores_chr#.txt)")
		ucsc_viz_parser.add_argument('outfile', type=str, action="store", help="file to write")
		ucsc_viz_parser.add_argument('--posIndex', type=int, action="store", default=1, help="index for column of datafile containing physical position (zero-indexed)")
		ucsc_viz_parser.add_argument('--scoreIndex', type=str, action="store", default=-2, help="index for column of datafile containing score (zero-indexed)")
		ucsc_viz_parser.add_argument('--strip_header', action="store_true", help="if input files include header line")

	################# 
	## SHARED ARGS ##
	################# 
	for commonparser in [xp_from_ihh_parser, composite_sims_parser, composite_emp_parser, normemp_genomewide_parser]:
		commonparser.add_argument('--cmsdir', help='TEMPORARY, will become redundant with conda packaging', action = 'store', default= "/idi/sabeti-scratch/jvitti/cms/cms/")
		commonparser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')

	for sim_parser in [composite_sims_parser, normsims_genomewide_parser]:
		sim_parser.add_argument('--nrep_sel', type= int, action='store', default='500')
		sim_parser.add_argument('--nrep_neut', type= int, action='store', default='1000')

	for composite_parser in [composite_sims_parser, composite_emp_parser]:
		composite_parser.add_argument('--likes_masterDir', type=str, default="/n/regal/sabeti_lab/jvitti/clear-synth/sims_reeval/likes_masters/", help="location of likelihood tables, defined")
		composite_parser.add_argument('--likes_nonSel', type=str, default="vsNeut", help='do we use completely neutral, or linked neutral SNPs for our non-causal distributions? by default, uses strict neutral (CMSgw)')
		composite_parser.add_argument('--likes_freqSuffix', type=str, default="allFreqs", help='for causal SNPs, include suffix to specify which selbins to include')
		composite_parser.add_argument('--cutoffline', type=str, default="250000\t1250000\t0\t1", help='specify bounds to include/exclude in calculations, along with MAF filter and likes decomposition')
		composite_parser.add_argument('--includeline', type=str, default="0\t0\t0\t0\t0\t0", help='specify (0:yes; 1:no) which scores to include: iHS,  ... ') #JV complete

	for cms_parser in [composite_sims_parser, composite_emp_parser, normsims_genomewide_parser, normemp_genomewide_parser]:
		cms_parser.add_argument('--writedir', type=str, default="", help="specify relative path") #ENFORCE CONSISTENCY - assumes e.g. model with sim scores live in this folder
		cms_parser.add_argument('--runSuffix', type=str, default=None, help='add a suffix to .cms file (corresponds to runparamfile)')

	for common_parser in [composite_sims_parser, composite_emp_parser, normsims_genomewide_parser]:
		common_parser.add_argument('--simpop', action='store', help='simulated population', default=1)

	for emp_parser in [composite_emp_parser, normemp_genomewide_parser]:
		emp_parser.add_argument('--emppop', action='store', help='empirical population', default="YRI")	

	for common_parser in [normemp_genomewide_parser, normsims_genomewide_parser, composite_sims_parser, composite_emp_parser]:
		common_parser.add_argument('--model', type=str, action="store", default="nulldefault")
		common_parser.add_argument('--checkOverwrite', action="store_true", default=False)

	#for commonparser in [composite_sims_parser, normsims_genomewide_parser, composite_emp_parser, normemp_genomewide_parser]:		=
	##	commonparser.add_argument('--suffix', type= str, action='store', default='')
	#	commonparser.add_argument('--simpop', action='store', help='simulated population', default=1)
	#	commonparser.add_argument('--emppop', action='store', help='empirical population', default="YRI")
	#	commonparser.add_argument('--model', type=str, default="nulldefault")
	#	commonparser.add_argument('--nrep', type=int, default=1000) #hmm remove for normemp_genomewide

	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
### Recalculate ancillary scores 
### from primary component scores
def execute_freqscores(args):
	''' python wrapper for program to calculate Fst and delDAF '''
	calc_fst_deldaf(args.inTped1, args.inTped2, args.recomFile, args.outfile, args.modelpath) #should obviate need for recom file. 
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
	if args.cmsdir is not None:
		cmd = args.cmsdir
	else:
		cmd = ""
	cmd += "combine/write_xpehh_from_ihh" #will become redundant with conda
	inputtped1 = args.inIhh1
	inputtped2 = args.inIhh2
	outfilename = args.outfilename
	argstring = inputtped1 + " " + inputtped2 + " " + outfilename
	cmdstring = cmd + " " + argstring
	if args.printOnly:
			print(cmdstring)
	else:
		subprocess.check_call( cmdstring.split() )	
	return

### Combine input component scores
### in user-defined CMS statistic
def execute_composite_sims(args):
	''' given simulated data and component scores (e.g. from likes_from_model.py) together with likelihood tables, generate CMS scores '''
	sel_freq_bins = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']

	if args.cmsdir is not None: 	#will be able to nix this construction 
		cmd = args.cmsdir 			#once conda packaging is complete
	else:							#but for now, keep things smooth.
		cmd = ""

	if args.regional_cms:
		cmd += "combine/combine_scores_local" 
		file_ending = ".cms.local.out"
	else:
		cmd += "combine/combine_scores_gw" #
		file_ending = ".cms.gw.out" 

	model = args.model
	selpop = args.simpop
	numPerBin_sel = args.nrep_sel
	numPerBin_neut = args.nrep_neut

	########################################
	## SPECIFY INPUT LIKELIHOOD FUNCTIONS ##
	########################################
	likes_masterDir = args.likes_masterDir
	likes_nonSel = args.likes_nonSel				
	likes_freqSuffix = args.likes_freqSuffix
	ihs_master, nsl_master, delihh_master, xpehh_master, fst_master, deldaf_master = get_master_likefiles(likes_masterDir, model, selpop, likes_nonSel, likes_freqSuffix)

	########################################################
	## RECORD INPUT PARAMETERS (scores, MAF filter, etc.) ##
	########################################################
	writedir = args.writedir
	cutoffline = args.cutoffline
	includeline = args.includeline 
	paramfilename = writedir + "run_params.txt"
	if args.runSuffix is not None:
		paramfilename += args.runSuffix
		suffix = args.runSuffix
	else:
		suffix = ""
	paramfilename = write_run_paramfile(paramfilename, ihs_master, nsl_master, delihh_master, xpehh_master, fst_master, deldaf_master, cutoffline, includeline)
	print("wrote CMS run parameters to: " + paramfilename)
	altpops = [1, 2, 3, 4]
	selpop = int(selpop)
	altpops.remove(selpop)

	##################################
	## CALCULATE CMS: ALL NEUT SIMS ##
	##################################
	scoremodeldir = writedir + model + "/neut/"
	compositedir = scoremodeldir + "composite/" 
	pairdir = scoremodeldir + "pairs/"
	check_create_dir(compositedir)
	check_create_dir(pairdir)
	for irep in range(1, numPerBin_neut +1):	
		altpairs = []
		for altpop in altpops:
			in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_sim_component_score_files(model, irep, selpop, altpop, selbin = "neut", filebase = writedir, normed = True)
			pairfilename = pairdir + "rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair" 
			if os.path.isfile(in_ihs_file) and os.path.isfile(in_nsl_file) and os.path.isfile(in_delihh_file) and os.path.isfile(in_xp_file) and os.path.isfile(in_fst_deldaf_file):
				write_pair_sourcefile(pairfilename, in_ihs_file, in_delihh_file, in_nsl_file, in_xp_file, in_fst_deldaf_file)
				altpairs.append(pairfilename)
		if len(altpairs) !=0:
			outfile = compositedir  + "rep" + str(irep) + "_" + str(selpop) + file_ending + suffix
			alreadyExists = False
			if args.checkOverwrite:
				if not os.path.isfile(outfile): #check for overwrite
					alreadyExists = False
				else:
					alreadyExists = True				
			if alreadyExists == False:
				argstring = outfile + " " + paramfilename + " "
				for pairfile in altpairs:
					argstring += pairfile + " "
				fullcmd = cmd + " " + argstring
				print(fullcmd)
				execute(fullcmd)
	
	#################################
	## CALCULATE CMS: ALL SEL SIMS ##
	#################################
	scoremodeldir = writedir + model + "/sel" + str(selpop) + "/"
	#check_create_dir(scoremodeldir)
	for sel_freq_bin in sel_freq_bins:
		this_bindir = scoremodeldir + "sel_" + str(sel_freq_bin) + "/"
		#check_create_dir(this_bindir)
		compositedir = this_bindir + "composite/" 
		pairdir = this_bindir + "pairs/"
		check_create_dir(compositedir)
		check_create_dir(pairdir)
		for irep in range(1, numPerBin_sel +1):
			altpairs = []
			for altpop in altpops:
				in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_sim_component_score_files(model, irep, selpop, altpop, selbin = sel_freq_bin, filebase = writedir, normed = True)
				pairfilename = pairdir + "rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair" 
				if os.path.isfile(in_ihs_file) and os.path.isfile(in_nsl_file) and os.path.isfile(in_delihh_file) and os.path.isfile(in_xp_file) and os.path.isfile(in_fst_deldaf_file):
					write_pair_sourcefile(pairfilename, in_ihs_file, in_delihh_file, in_nsl_file, in_xp_file, in_fst_deldaf_file)
					altpairs.append(pairfilename)
			if len(altpairs) !=0:
				outfile = compositedir + "rep" + str(irep) + "_" + str(selpop) + file_ending + suffix
				alreadyExists = False
				if args.checkOverwrite:
					if not os.path.isfile(outfile): #check for overwrite
						alreadyExists = False
					else:
						alreadyExists = True				
				if alreadyExists == False:
					argstring = outfile + " " + paramfilename + " "
					for pairfile in altpairs:
						argstring += pairfile +" "
					fullcmd = cmd + " " + argstring
					print(fullcmd)
					execute(fullcmd)
	print('calculated CMS scores for ' + str(numPerBin_neut) + ' neutral replicates and ' + str(numPerBin_sel) + " selection replicates per bin.")
	return
def execute_composite_emp(args):
	''' given component scores from empirical data (e.g. from scans.py) together with likelihood tables, generate CMS scores '''
	model_popsdict = {1:["YRI", "AFR", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
						2:["CEU", "EUR", "TSI", "FIN", "GBR", "IBS", "IRN"],
						3:["CHB", "EAS", "JPT", "CHS", "CDX", "KHV"],
						4:["BEB", "SAS", "GIH", "PJL", "STU", "ITU"],
						0:["MXL", "AMR", "PUR", "CLM", "PEL"]} #American populations excluded from model

	if args.cmsdir is not None: 	#will be able to nix this construction 
		cmd = args.cmsdir 			#once conda packaging is complete
	else:							#but for now, keep things smooth.
		cmd = ""

	if args.regional_cms_chrom is None:
		cmd += "combine/combine_scores_gw" #genome-wide
		chroms = range(1,23)
		file_ending = ".cms.gw.out"
	else:
		cmd += "combine/combine_scores_local" #within-region
		chroms = [args.regional_cms_chrom]
		file_ending = ".cms.local.out"

	model = args.model
	modelPop = args.simpop
	score_basedir = args.score_basedir

	########################################
	## SPECIFY INPUT LIKELIHOOD FUNCTIONS ##
	########################################
	likes_masterDir = args.likes_masterDir
	likes_nonSel = args.likes_nonSel				
	likes_freqSuffix = args.likes_freqSuffix
	ihs_master, nsl_master, delihh_master, xpehh_master, fst_master, deldaf_master = get_master_likefiles(likes_masterDir, model, modelPop, likes_nonSel, likes_freqSuffix)
	#build in a check here to enforce correct likes for within-region vs. genomew-wide?

	#########################
	## DESIGNATE OUTGROUPS ##
	#########################
	emp_selpop = args.emppop
	model_selpop = 0
	for modelpop in [1, 2, 3, 4, 0]:
		if emp_selpop in model_popsdict[modelpop]:
			model_selpop = modelpop
	if model_selpop == 0: 
		model_selpop = 4

	altmodelpops = [1, 2, 3, 4]
	altmodelpops.remove(model_selpop)
	altpops = []
	for altpop in altmodelpops:
		altpops.append(model_popsdict[altpop][0])

	########################################################
	## RECORD INPUT PARAMETERS (scores, MAF filter, etc.) ##
	########################################################
	writedir = args.writedir
	cutoffline = args.cutoffline
	includeline = args.includeline
	paramfilename = writedir + "run_params.txt"
	if args.runSuffix is not None:
		paramfilename += args.runSuffix
		suffix = args.runSuffix
	else:
		suffix = ""
	paramfilename += "_" + str(model_selpop)
	paramfilename = write_run_paramfile(paramfilename, ihs_master, nsl_master, delihh_master, xpehh_master, fst_master, deldaf_master, cutoffline, includeline)
	print("wrote CMS run parameters to: " + paramfilename)

	#############################################
	## CALCULATE CMS: ITERATE OVER CHROMOSOMES ##
	#############################################
	for chrom in chroms:
		altpairs = []
		for altpop in altpops:
			in_ihs_file, in_nsl_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_emp_component_score_files(chrom, emp_selpop, altpop=altpop, basedir = score_basedir) 
			pairfilename = score_basedir + "pairs/chr" + str(chrom) + "_" + str(emp_selpop) + "_" + str(altpop) + ".pair" 
			if os.path.isfile(in_ihs_file) and os.path.isfile(in_nsl_file) and os.path.isfile(in_delihh_file) and os.path.isfile(in_xp_file) and os.path.isfile(in_fst_deldaf_file):
				write_pair_sourcefile(pairfilename, in_ihs_file, in_delihh_file, in_nsl_file, in_xp_file, in_fst_deldaf_file)
				altpairs.append(pairfilename)
		if len(altpairs) !=0:
			outfile = score_basedir + "composite/"
			if args.regional_cms_chrom is not None:
				outfile += "regional/"
			outfile += "chr" + str(chrom) + "_" + str(emp_selpop) + file_ending + suffix
			alreadyExists = False
			if args.checkOverwrite:
				if not os.path.isfile(outfile): #check for overwrite
					alreadyExists = False
				else:
					alreadyExists = True				
			if alreadyExists == False:
				argstring = outfile + " " + paramfilename + " "
				for pairfile in altpairs:
					argstring += pairfile + " "
				fullcmd = cmd + " " + argstring
				print(fullcmd)
				execute(fullcmd)	
	print('calculated CMS scores for ' + str(len(chroms)) + ' chromosomes.')

	return
def execute_normsims_genomewide(args): 
	""" given output from composite_sims, normalize all replicates to neutral parameters """ 
	sel_freq_bins = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
	model = args.model
	selpop = args.simpop
	numPerBin_sel = args.nrep_sel
	numPerBin_neut = args.nrep_neut
	writedir = args.writedir
	suffix = args.runSuffix
	
	values = []
	##############################
	## LOAD STATS FROM NEUT SIMS #
	##############################
	for irep in range(1, numPerBin_neut +1):	
		outfile  = get_neut_repfile_name(model, irep, selpop, suffix = suffix, normed=False, basedir=writedir)
		if os.path.isfile(outfile):
			openfile = open(outfile, 'r')
			header = openfile.readline() 
			for line in openfile:
				entries = line.split()
				rawscore = np.log(float(entries[-1]))
				values.append(rawscore)
			openfile.close()
		else:
			print('missing: ' + outfile)

	print('loaded ' + str(len(values)) + ' values from neutral sims...')

	#check for nans
	values = np.array(values)
	values = values[~np.isnan(values)]
	values = list(values)

	#check for infs
	values = np.array(values)
	values = values[~np.isinf(values)]
	values = list(values)

	mean = np.mean(values)
	var = np.var(values)
	sd = np.sqrt(var)

	print("max: " + str(max(values)))
	print("min: " + str(min(values)))
	print("mean: " + str(np.mean(values)))
	print("var: " + str(np.var(values)))

	############################
	## NORMALIZE NEUTRAL SIMS ##
	############################
	
	for irep in range(1, numPerBin_neut +1):	
		outfile  = get_neut_repfile_name(model, irep, selpop, suffix = suffix, normed=False, basedir=writedir)
		if os.path.isfile(outfile):
			normedfile = outfile + ".norm"
			if True:
			#if not os.path.isfile(normedfile): #CHANGE FOR --checkOverwrite
				openfile = open(outfile, 'r')
				writefile = open(normedfile, 'w')
				header = openfile.readline()
				writefile.write(header)
				for line in openfile:
					entries = line.split()
					rawscore = np.log(float(entries[-1]))
					normalized = normalize(rawscore, mean, sd)
					writeline = line.strip('\n') + "\t" + str(normalized)+ "\n"
					writefile.write(writeline)
				openfile.close()
				writefile.close()
	print("wrote to eg: " + normedfile)	
	
	########################
	## NORMALIZE SEL SIMS ##
	########################
	for sel_freq_bin in sel_freq_bins:
		for irep in range(1, numPerBin_sel +1):
			rawfile = get_sel_repfile_name(model, irep, selpop, sel_freq_bin, suffix=suffix, normed = False, basedir=writedir)
			#print(rawfile)
			if os.path.isfile(rawfile):
				normedfile = rawfile + ".norm"
				if True:
				#if not os.path.isfile(normedfile):
					openfile = open(rawfile, 'r')
					writefile = open(normedfile, 'w')
					header = openfile.readline()
					writefile.write(header)
					for line in openfile:
						entries = line.split()
						rawscore = np.log(float(entries[-1]))
						normalized = normalize(rawscore, mean, sd)
						writeline = line.strip('\n') + "\t" + str(normalized) + "\n"
						writefile.write(writeline)
					openfile.close()
					writefile.close()
	print("wrote to eg: " + normedfile)	
	return
def execute_normemp_genomewide(args):
	""" given output from composite_emp, normalize CMS scores genome-wide """  #could also introduce a feature to normalize to explicitly neutral regions. 
	selpop = args.emppop
	model = args.model
	#suffix = args.suffix
	#basedir = args.writedir
	score_basedir = args.score_basedir

	if args.runSuffix is not None:
		suffix = args.runSuffix
	else:
		suffix = ""

	print('must check correct load_emp function')
	scores2 = load_empscores(model, selpop, normed=False, suffix=suffix, basedir=score_basedir) #HANDLE THIS MORE EXPLICITLY
	scores = [np.log(item) for item in scores2]

	#check for nans
	scores = np.array(scores)
	scores = scores[~np.isnan(scores)]
	scores = list(scores)

	#check for infs
	scores = np.array(scores)
	scores = scores[~np.isinf(scores)]
	scores = list(scores)
	
	print('loaded ' + str(len(scores)) + " scores")
	print("max: " + str(max(scores)))
	print("min: " + str(min(scores)))
	print("mean: " + str(np.mean(scores)))
	print("var: " + str(np.var(scores)))

	mean = np.mean(scores)
	var = np.var(scores)
	sd = np.sqrt(var)

	##############
	## NORMALIZE #
	############## 
	chroms = range(1,23)
	for thischrom in chroms:
		unnormedfile = get_emp_cms_file(selpop, thischrom, normed=False, basedir=score_basedir, suffix=suffix,) #model #selpop, chrom, normed=False, suffix=suffix, basedir = score_basedir)
		assert os.path.isfile(unnormedfile)
		normedfile = unnormedfile + ".norm"

		readfile = open(unnormedfile, 'r')
		writefile = open(normedfile, 'w')
		for line in readfile:
			line = line.strip('\n')
			entries=line.split()
			rawscore = np.log(float(entries[-1]))
			normedscore = normalize(rawscore, mean, sd)
			writeline = line + "\t" + str(normedscore) + '\n'
			writefile.write(writeline)
		readfile.close()
		writefile.close
		print('wrote to '  + normedfile)
	return

### Visualize and hone in 
### on variants within regions
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
