#!/usr/bin/env python
## top-level script for demographic modeling as part of CMS 2.0. 
## last updated: 02.26.18 vitti@broadinstitute.org
##	10.19.18: fitting ASW with custom tree topology.
##	12.17.18 fitting PEL

import matplotlib as mp 
mp.use('agg') 
from model.bootstrap_func import flattenList, checkFileExists, readFreqsFile, readLDFile, readFstFile, estimateFstByBootstrap, estimateFstByBootstrap_bysnp, estimateFreqSpectrum, estimatePi, estimater2decay, estimatedprimedecay
from model.params_func import get_ranges, generate_params, get_dict_from_paramfile, write_paramfile
#from params_func_ASW import get_ranges, generate_params
from model.error_func import calc_error, read_error_dimensionsfile
from model.search_func import read_dimensionsfile, sample_point, get_real_value, get_scaled_value
from model.plot_func import plot_comparison
from model.draw_tree import draw_tree
from scipy import optimize
import numpy as np
import subprocess
import argparse
import random
import sys
from itertools import product

def update_paramdict(paramDict, startPointFileName):
	#quick hack
	openfile = open(startPointFileName, 'r')
	for line in openfile:
		entries = line.split('\t')
		if "(" in entries[0]:
			key = eval(entries[0])
		else:
			key = entries[0]
		index = int(entries[1])
		value = float(entries[2])
		paramDict[key][index] = value
	openfile.close()
	return paramDict

#############################
## DEFINE ARGUMENT PARSER ###
#############################
def full_parser_cms_modeller():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for exploratory fitting of demographic models to population genetic data.")
	subparsers = parser.add_subparsers(help="sub-commands")

	######################
	## CALCULATE TARGET ##
	######################
	target_stats_parser = subparsers.add_parser('target_stats', help='Perform per-site(/per-site-pair) calculations of population summary statistics for model target values.')
	target_stats_parser.add_argument('inputTpeds', action='store', type=list, help='comma-delimited list of unzipped input tped files (only one file per pop being modelled; must run chroms separately or concatenate)')
	target_stats_parser.add_argument('recomFile', action='store', type=str, help='file defining recombination map for input') 
	target_stats_parser.add_argument('regions', action='store', type=str, help='tab-separated file with putative neutral regions') #OPTIONAL?
	target_stats_parser.add_argument('--freqs', action='store_true', help='calculate summary statistics from within-population allele frequencies') 
	target_stats_parser.add_argument('--ld', action='store_true', help='calculate summary statistics from within-population linkage disequilibrium') 
	target_stats_parser.add_argument('--fst', action='store_true', help='calculate summary statistics from population comparison using allele frequencies') 
	target_stats_parser.add_argument('out', action='store', type=str, help='outfile prefix') 
	target_stats_parser.add_argument('--modelpath', action='store', type=str, default='cms/model/', help="path to model directory containing executables")

	bootstrap_parser = subparsers.add_parser('bootstrap', help='Perform bootstrap estimates of population summary statistics from per-site(/per-site-pair) calculations in order to finalize model target values.')
	bootstrap_parser.add_argument('nBootstrapReps', action='store', type=int, help='number of bootstraps to perform in order to estimate standard error of the dataset (should converge for reasonably small n)')
	bootstrap_parser.add_argument('--in_freqs', action='store', help='comma-delimited list of infiles with per-site calculations for population. One file per population -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files') 
	bootstrap_parser.add_argument('--nFreqHistBins', action='store',type=int, default=6, help="number of bins for site frequency spectrum and p(der|freq)")
	bootstrap_parser.add_argument('--in_ld', action='store', help='comma-delimited list of infiles with per-site-pair calculations for population. One file per population -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files') 
	bootstrap_parser.add_argument('--mafcutoffdprime', action='store', type=float, default=.2, help="for D' calculations, only use sites with MAF > mafcutoffdprime") 	
	bootstrap_parser.add_argument('--nphysdisthist', action='store', type=int, default=14, help="nbins for r2 LD calculations") 	


	bootstrap_parser.add_argument('--in_fst', action='store', help='comma-delimited list of infiles with per-site calculations for population pair. One file per population-pair -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files') 	
	bootstrap_parser.add_argument('--ngendisthist', action='store', type=int, default=17, help="nbins for D' LD calculations") 	
	bootstrap_parser.add_argument('out', action='store', type=str, help='outfile prefix') 
	
	##########################
	### COSI - SHARED ARGS  ##
	##########################
	point_parser = subparsers.add_parser('point', help='Run simulates of a point in parameter-space.')
	grid_parser = subparsers.add_parser('grid', help='Perform grid search: for specified parameters and intervals, define points in parameter-space to sample and compare.')
	optimize_parser = subparsers.add_parser('optimize', help='Perform optimization algorithm (scipy.optimize) to fit model parameters robustly.')

	point_parser.add_argument('--inputParamFile', type=str, action='store', help='file with model specifications for input')
		
	for cosi_parser in [point_parser, grid_parser, optimize_parser]:
		#cosi_parser.add_argument('inputParamFile', type=str, action='store', help='file with model specifications for input')
		cosi_parser.add_argument('nCoalescentReps', type=int, help='number of coalescent replicates to run per point in parameter-space')
		cosi_parser.add_argument('outputDir', type=str, action='store', help='location in which to write cosi output')
		cosi_parser.add_argument('--cosiBuild', action='store', default="coalescent", help='which version of cosi to run?')  
		cosi_parser.add_argument('--dropSings', action='store', type=float, help='randomly thin global singletons from output dataset (i.e., to model ascertainment bias)')
		cosi_parser.add_argument('--genmapRandomRegions', action='store_true', help='cosi option to sub-sample genetic map randomly from input')
		cosi_parser.add_argument('--stopAfterMinutes', action='store', help='cosi option to terminate simulations')
		cosi_parser.add_argument('--calcError', type=str, action='store', help='file specifying dimensions of error function to use. if unspecified, defaults to all. first line = stats, second line = pops')

	######################
	## VISUALIZE MODEL  ##
	######################
	point_parser.add_argument('--targetvalsFile', action='store', type=str, help='file containing target values for model')	
	point_parser.add_argument('--plotStats', action='store_true', default=False, help='visualize goodness-of-fit to model targets')
	point_parser.add_argument('--plotbase', action="store", default="/web/personal/vitti/test_dem")
	#########################
	## FIT MODEL TO TARGET ##
	#########################
	grid_parser.add_argument('grid_inputdimensionsfile', type=str, action='store', help='file with specifications of grid search. each parameter to vary is indicated: KEY\tINDEX\t[VALUES]') #must be defined for each search 	
	#grid_parser.add_argument('--parallel', type=str, action='store', default="uger", help='if specified, launch points of grid search as tasks on a scheduler')
	optimize_parser.add_argument('optimize_inputdimensionsfile', type=str, action='store', help='file with specifications of optimization. each parameter to vary is indicated: KEY\tINDEX')
	optimize_parser.add_argument('--stepSize', action='store', type=float, help='scaled step size (i.e. whole range = 1)')
	optimize_parser.add_argument('--method', action='store', type=str, default='SLSQP', help='algorithm to pass to scipy.optimize')
	for common_parser in [optimize_parser, grid_parser]:
		common_parser.add_argument('--savePar', action='store', type=str, default="opt_out.par")
	for common_parser in [target_stats_parser, point_parser]:#[bootstrap_parser, grid_parser, optimize_parser]:
		common_parser.add_argument('--printOnly', action='store_true', help='print rather than execute pipeline commands')
	for common_parser in [point_parser, optimize_parser]:
		common_parser.add_argument('--startPointFileName', action='store', type=str, default=None, help='start the gradient at a point in parameter-space with these specififications: KEY\tINDEX\tVALUE')
	
	return parser

############################
## DEFINE EXEC FUNCTIONS ###
############################
def execute_target_stats(args):
	'''calls bootstrap_*_popstats_regions to get per-snp/per-snp-pair values'''
	#pathcmd = "export PATH=" + args.modelpath + ":$PATH"
	#subprocess.check_call(pathcmd.split())
	modelpath = args.modelpath
	if modelpath[-1] != "/":
		modelpath += "/"
	inputtpedstring = ''.join(args.inputTpeds)
	inputtpeds = inputtpedstring.split(',')
	npops = len(inputtpeds)
	print("calculating summary statistics for " +  str(npops) + " populations...")
	allCmds = []
	for ipop in range(npops):
		inputtped = inputtpeds[ipop]
		if args.freqs:
			freqCmd = [modelpath + 'bootstrap_freq_popstats_regions', inputtped, args.recomFile, args.regions, args.out + "_freqs_" + str(ipop)]
			allCmds.append(freqCmd)
		if args.ld:
			ldCmd = [modelpath + 'bootstrap_ld_popstats_regions', inputtped, args.recomFile, args.regions, args.out + "_ld_" + str(ipop)]
			allCmds.append(ldCmd)
		if args.fst:
			for jpop in range(ipop+1, npops):
				inputtped2 = inputtpeds[jpop]
				fstCmd = [modelpath + 'bootstrap_fst_popstats_regions', inputtped, inputtped2, args.recomFile, args.regions, args.out + "_fst_" + str(ipop) + "_" + str(jpop)]
				allCmds.append(fstCmd)
	for command in allCmds:
		command = [str(x) for x in command]
		if args.printOnly:
			commandstring = ""
			for item in command:
				commandstring += item + " "
			print(commandstring)
		else:
			subprocess.check_call( command )
	return
def execute_bootstrap(args):
	'''pulls all per-snp/per-snp-pair values to get genome-wide bootstrap estimates.'''
	nbootstraprep = args.nBootstrapReps
	print("running " + str(nbootstraprep) + " bootstrap estimates of summary statistics...")
	targetstats_filename = args.out + "_bootstrap_n" + str(nbootstraprep) + ".txt"
	writefile = open(targetstats_filename, 'w')

	#################
	### FREQ STATS ##
	#################
	if args.in_freqs is not None: 
		nhist = args.nFreqHistBins
		inputestimatefilenames = ''.join(args.in_freqs)
		inputfilenames = inputestimatefilenames.split(',')
		npops = len(inputfilenames)
		for ipop in range(npops):
			allRegionDER, allRegionANC, allRegionPI, allseqlens = [], [], [], []
			nsnps, totalregions, totallen = 0, 0, 0
			inputfilename = inputfilenames[ipop]
			print("reading allele frequency statistics from: " + inputfilename)
			writefile.write(str(ipop) + '\n')
			if checkFileExists(inputfilename):
				allpi, allnderiv, allnanc, nregions, seqlens = readFreqsFile(inputfilename)
				allRegionPI.extend(allpi)
				allRegionDER.extend(allnderiv)
				allRegionANC.extend(allnanc)
				allseqlens.extend(seqlens)
				totalregions += nregions
				totallen += sum(seqlens)
				for i in range(len(allpi)):
					nsnps += len(allpi[i])
			print("TOTAL: logged frequency values for " + str(nsnps) + " SNPS across " + str(totalregions) + ".")
			
			####################################
			#### PI: MEAN & BOOTSTRAP STDERR ###
			####################################
			pi_mean = estimatePi(allRegionPI, allseqlens)
			writefile.write(str(pi_mean)+'\t')
			estimates = []
			for j in range(nbootstraprep):
				rep_pis, rep_seqlens = [], []
				for k in range(totalregions):
					index = random.randint(0, totalregions-1)
					rep_pis.append(allRegionPI[index])
					rep_seqlens.append(allseqlens[index])
				rep_pi_mean = estimatePi(rep_pis, rep_seqlens)
				estimates.append(rep_pi_mean)
			pi_se = np.std(estimates)
			writefile.write(str(pi_se) + '\n')

			#########################################
			### SFS, ANC: MEAN ACROSS ALL REGIONS ###
			#########################################
			mafhist, anchist = estimateFreqSpectrum(allRegionDER, allRegionANC, nhist)
			npoly = sum(mafhist)
			sfs_mean = [float(x)/npoly for x in mafhist]
			anc_mean = [anchist[i]/float(mafhist[i]) for i in range(len(mafhist))]

			###################################################
			### SFS, ANC: STDERR ACROSS BOOTSTRAP ESTIMATES ###
			###################################################
			estimates_sfs, estimates_anc = [[] for i in range(nhist)], [[] for i in range(nhist)]
			for j in range(nbootstraprep):
				rep_all_nderiv, rep_all_nanc = [], []
				flatanc = flattenList(allRegionANC)
				flatder = flattenList(allRegionDER)
				for w in range(nsnps):
					index = random.randint(0, nsnps-1)
					rep_all_nderiv.append(flatder[index])
					rep_all_nanc.append(flatanc[index])
				repmafhist, repanchist = estimateFreqSpectrum(rep_all_nderiv, rep_all_nanc, nhist)
				npoly = sum(repmafhist)
				repsfs = [float(x)/npoly for x in repmafhist]
				for ibin in range(nhist):
					estimates_sfs[ibin].append(repsfs[ibin])
				repanc = [repanchist[i]/float(repmafhist[i]) for i in range(nhist)]
				for ibin in range(nhist):
					estimates_anc[ibin].append(repanc[ibin])

			sfs_se = [np.std(x) for x in estimates_sfs]
			anc_se = [np.std(x) for x in estimates_anc]
			writefile.write(str(sfs_mean) + '\n')
			writefile.write(str(sfs_se) + '\n')
			writefile.write(str(anc_mean) + '\n')
			writefile.write(str(anc_se) + '\n')
		
	#########
	### LD ##
	#########
	if args.in_ld is not None:
		nphysdisthist = args.nphysdisthist
		ngendisthist = args.ngendisthist
		inputestimatefilenames = ''.join(args.in_ld)
		inputfilenames = inputestimatefilenames.split(',')
		npops = len(inputfilenames)
		#print('npops ' + str(npops)) #debug
		for ipop in range(npops):
			inputfilename = inputfilenames[ipop]
			print("reading linkage disequilibrium statistics from: " + inputfilename)
			writefile.write(str(ipop) + '\n')
			N_r2regs, N_dprimeregs = 0, 0
			N_r2snps, N_dprimesnps = 0, 0
			allRegionDists, allRegionr2, allRegionGendists, allRegionDprime, nr2regions, ndprimeregions = readLDFile(inputfilename, dprimecutoff = args.mafcutoffdprime)
			N_r2regs += nr2regions
			N_r2snps += sum([len(x) for x in allRegionr2])
			N_dprimeregs += ndprimeregions
			N_dprimesnps += sum([len(x) for x in allRegionDprime])
			print("\tlogged r2 values for " + str(N_r2snps) + " SNP pairs across " + str(N_r2regs) + " regions.")
			print("\tlogged D' values for " + str(N_dprimesnps) + " SNP pairs across " + str(N_dprimeregs) + " regions.")


			###################################
			### r2: MEAN ACROSS ALL REGIONS ###
			###################################
			r2sums, physDistHist = estimater2decay(allRegionr2, allRegionDists)#, nphysdisthist)
			r2dist = [r2sums[u]/physDistHist[u] for u in range(len(r2sums))]
			writefile.write(str(r2dist) + "\n")

			############################################
			### r2: STDERR ACROSS BOOTSTRAP ESTIMATES ##
			############################################
			estimates_r2 = [[] for i in range(nphysdisthist)]
			while len(estimates_r2[0]) < nbootstraprep:
				rep_all_r2, rep_all_physdist = [], []
				flatr2 = flattenList(allRegionr2)
				flatregions = flattenList(allRegionDists)
				nsnppairs = len(flatr2)
				for w in range(nsnppairs):
					index_r2 = random.randint(0, nsnppairs-1)
					rep_all_r2.append(flatr2[index_r2])
					rep_all_physdist.append(flatregions[index_r2])

				#add pseudocount for empty bins
				repr2sum, repphysdisthist = estimater2decay(rep_all_r2, rep_all_physdist)#, nphysdisthist)
				for ibin in range(len(repphysdisthist)):
					if repphysdisthist[ibin] == 0:
						repphysdisthist[ibin] = 1
				r2estimate =[repr2sum[u]/repphysdisthist[u] for u in range(len(repr2sum))]
				for ibin in range(nphysdisthist):
					estimates_r2[ibin].append(r2estimate[ibin])

			r2_se = [np.std(x) for x in estimates_r2]
			writefile.write(str(r2_se) + "\n")

			####################################
			### D': MEAN ACROSS ALL REGIONS ###
			####################################
			compLDhist, genDistHist = estimatedprimedecay(allRegionDprime, allRegionGendists)#, ngendisthist)
			#add pseudocounts
			for ibin in range(len(genDistHist)):
				if genDistHist[ibin] == 0:
					genDistHist[ibin]+=1

			dprimedist = [float(compLDhist[x])/float(genDistHist[x]) for x in range(len(compLDhist))]
			writefile.write(str(dprimedist) + "\n")

			############################################
			### D': STDERR ACROSS BOOTSTRAP ESTIMATES ##
			############################################
			estimates_dprime = [[] for i in range(ngendisthist)]
			while len(estimates_dprime[0]) < nbootstraprep:
				rep_all_dprime, rep_all_gendist = [], []
				flatdprime = flattenList(allRegionDprime)
				flatgendist = flattenList(allRegionGendists)
				nsnppairs = len(flatdprime)

				for w in range(nsnppairs):
					index_dprime = random.randint(0, nsnppairs-1)
					rep_all_dprime.append(flatdprime[index_dprime])
					rep_all_gendist.append(flatgendist[index_dprime])

				repcompLDhist, repgenDistHist = estimatedprimedecay(rep_all_dprime, rep_all_gendist)#, ngendisthist)
				for ibin in range(len(repgenDistHist)):
					if repgenDistHist[ibin] == 0:
						repgenDistHist[ibin] = 1
				dprimeestimate = [float(repcompLDhist[x])/float(repgenDistHist[x]) for x in range(ngendisthist)]
				for ibin in range(ngendisthist):
					estimates_dprime[ibin].append(dprimeestimate[ibin])
			dprime_se = [np.std(x) for x in estimates_dprime]
			writefile.write(str(dprime_se) + "\n")

	##########
	### FST ##
	##########
	if args.in_fst is not None:
		inputestimatefilenames = ''.join(args.in_fst)
		inputfilenames = inputestimatefilenames.split(',')
		npopcomp = len(inputfilenames)
		for icomp in range(npopcomp):
			fstfilename	= inputfilenames[icomp]
			print("reading Fst values from: " + fstfilename)
			if checkFileExists(fstfilename):
				allfst, nregions = readFstFile(fstfilename)
			else:
				print('missing ' + fstfilename)
			target_mean, target_se = estimateFstByBootstrap_bysnp(allfst, nrep = nbootstraprep)
			writeline =  str(icomp) + "\t" + str(target_mean) + "\t" + str(target_se) + '\n'
			writefile.write(writeline)
			print("TOTAL: logged Fst values for " + str(len(allfst)) + " SNPs.\n")

	writefile.close()
	print("wrote to file: " + targetstats_filename)
	return
def execute_point(args):
	'''runs simulates of a point in parameter-space, comparing to specified target'''
	################
	## FILE PREP ###
	################
	if args.inputParamFile is not None:
		print("generating " + str(args.nCoalescentReps) + " simulations from model: " + args.inputParamFile)
		inputParamFile = args.inputParamFile
		paramDict = get_dict_from_paramfile(inputParamFile)

	else:
		inputParamFile = "test_args.par"
		#args.inputParamFile
		paramDict = generate_params()
		write_paramfile(inputParamFile, paramDict) #MAKE TWEAKABLE FROM CMDLINE

	statfilename = args.outputDir
	if args.outputDir[-1] != "/":
		statfilename += "/"
	statfilename += "n" + str(args.nCoalescentReps) + "stats.txt"

	draw_tree(paramDict, "/web/personal/vitti/test_tree.png")
	draw_tree(paramDict, "/web/personal/vitti/test_tree_recent.png", recent=True)
	draw_tree(paramDict, "/web/personal/vitti/test_tree_very_recent.png", veryrecent=True)

	###############
	## RUN SIMS ###
	###############
	runStatsCommand = args.cosiBuild + " -p " + inputParamFile + " -n " + str(args.nCoalescentReps) 
	if args.dropSings is not None:
		runStatsCommand += " --drop-singletons " + str(args.dropSings)
	if args.genmapRandomRegions:
		runStatsCommand += " --genmapRandomRegions"
	if args.stopAfterMinutes is not None:
		runStatsCommand += " --stop-after-minutes " + str(args.stopAfterMinutes)
	runStatsCommand += " --custom-stats "#"> " + statfilename

	print(runStatsCommand)
	if args.printOnly:
		print(runStatsCommand)
	else:
		with open(statfilename, 'w') as outfile:
			subprocess.run( runStatsCommand.split(), stdout=outfile)
		outfile.close()
		#print(x)
	#################
	## CALC ERROR ###
	#################
	if args.calcError is not None:
		print('calculating error from ' + statfilename)
		if args.calcError == '': #no error dimension file given
			error = calc_error(statfilename)
		else:
			stats, pops = read_error_dimensionsfile(args.calcError) 
			print('stats: ', stats)
			print('pops: ', pops)
			error = calc_error(statfilename, stats, pops)
	else:
		error = calc_error(statfilename)
	print(" error: " + str(error)) #record?
		
	################
	## VISUALIZE ###
	################		
	if args.plotStats:
		print('plotting stats from: ' + statfilename)
		plot_comparison(statfilename, args.nCoalescentReps, args.plotbase)

	return
def execute_grid(args):
	'''run points in parameter-space according to specified grid'''
	if args.calcError is not None:
		#print('calculating error from ' + statfilename)
		#if args.calcError == '': #no error dimension file given
		#	error = calc_error(statfilename)
		#else:
		stats, pops = read_error_dimensionsfile(args.calcError) 
	else:
		stats=['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
		pops = [1, 2, 3, 4]

	writeparamfilename = args.savePar
	print("loading dimensions of grid to search from: " + args.grid_inputdimensionsfile)
	gridname, keys, indices, values = read_dimensionsfile(args.grid_inputdimensionsfile, 'grid', probTakeSearchParam=100) #TAKE ALL FOR MANUALLY SPECIFIED GRID
	assert len(keys) == len(indices) 
	valuecombos = product(*values)
	combos =  [' '.join(str(y) for y in x) for x in valuecombos]
	errors = []
	for combo in combos:
		argstring = combo + "\n"
		theseValues = [float(item) for item in combo.split(' ')]#eval(combo) #list of values
		print('\n$$$$$$$$$$\n')
		print(stats, pops)
		print(keys, indices, theseValues)
		print(writeparamfilename)
		print('start call to SAMPLE POINT...')
		error = sample_point(args.nCoalescentReps, keys, indices, theseValues, writeparamfilename, False, stats, pops) # add update_paramdict in here?
		errors.append(error)


	print('\n******************')
	for icombo in range(len(combos)):
		print(str(combos[icombo]) + "\t" + str(errors[icombo]))
	return
def execute_optimize(args):
	'''run scipy.optimize module according to specified parameters'''
	if args.calcError is not None:
		#print('calculating error from ' + statfilename)
		#if args.calcError == '': #no error dimension file given
		#	error = calc_error(statfilename)
		#else:
		stats, pops = read_error_dimensionsfile(args.calcError) 
		print(stats, pops)
	else:
		stats=['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
		pops = [1, 2, 3, 4]

	print("loading dimensions to search from: " + args.optimize_inputdimensionsfile)
	runname, keys, indices = read_dimensionsfile(args.optimize_inputdimensionsfile, runType='optimize', probTakeSearchParam=20)

	rangeDict = get_ranges()
	paramDict = generate_params() #<- this is the starting point for gradient. make it customizeable?
	if args.startPointFileName is not None:
		paramDict = update_paramdict(paramDict, args.startPointFileName)

	x0 = []
	bounds = []
	values = []
	scaleds = []
	for i in range(len(keys)):
		key = keys[i]
		index = indices[i]
		value = paramDict[key][index]
		interval = rangeDict[key][index]
		low, high = float(interval[0]), float(interval[1])
		values.append(value)
		scaled = get_scaled_value(value, low, high)
		x0.append(scaled)
		bounds.append([0,1])
		scaleds.append(scaled)

	def sample_point_wrapper(values):
		'''function passed to scipy.optimize.'''
		#SET THESE UP TO GET FROM FILE!!!
		for i in range(len(values)):
			if values[i] > 1:
				values[i] = 1
			if values[i] < 0:
				values[i] = 0 #ENFORCE BOUNDS
		writeparamfilename = args.savePar#args.optimize_inputdimensionsfile + "_explore.par" #<--LETS TRY THIS FOR NOW TO TRACK BETTER
		return sample_point(args.nCoalescentReps, keys, indices, values, writeparamfilename, True, stats, pops)


	x0 = np.array(x0)
	stepdict = {'eps':float(args.stepSize)}
	result = optimize.minimize(sample_point_wrapper, x0, method=args.method, bounds=bounds, options=stepdict)
	print(result)

	print( "******************")
	#translate back to changes to model
	bestparams = []
	assert len(keys) == len(result.x)
	for i in range(len(keys)):
		key = keys[i]
		index = indices[i]
		interval = rangeDict[key][index]
		low, high = float(interval[0]), float(interval[1])
		realVal = get_real_value(result.x[i], low, high)
		bestparams.append(result.x[i])
		print("best " + str(key) + "|" + str(index) + "|" + str(realVal))
	return

##########
## MAIN ##
##########
if __name__ == '__main__':
	runparser = full_parser_cms_modeller()
	args = runparser.parse_args()

	# if called with no arguments, print help
	if len(sys.argv)==1:
		runparser.parse_args(['--help'])
	elif len(sys.argv)==2:
		runparser.parse_args([sys.argv[1], '--help'])

	subcommand = sys.argv[1]
	function_name = 'execute_' + subcommand + "(args)"
	eval(function_name) #points to functions defined above, which wrap other programs in the pipeline

