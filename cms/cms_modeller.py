## top-level script for demographic modeling as part of CMS 2.0. 
#(CURRENTLY: assumes C programs have been compiled in same directory -- JV must provide Makefile; assumes python in $PATH)
## last updated: 07.04.16 vitti@broadinstitute.org

prefixstring = "{CMS2.0}>>\t\t" #for stderr (make global?)
from model.bootstrap_func import *
from model.params import getTargetValues, generateParams, writeParamFile, getDictFromParamFile, getRanges, updateParams, getScaledValues, getScaledValue, getRealValue
from model.parse import givecommand, readlines, checkFailed, readParameterValues, readGlobalArgFile, readOutfile, readcustomstatfile
from model.error import RMSE, calcError, samplePoint

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
	target_stats_parser.add_argument('regions', action='store', help='putative neutral regions') #make this optional? and/or TSV
	target_stats_parser.add_argument('--freqs', action='store_true', help='calculate summary statistics from within-population allele frequencies') 
	target_stats_parser.add_argument('--ld', action='store_true', help='calculate summary statistics from within-population linkage disequilibrium') 
	target_stats_parser.add_argument('--fst', action='store_true', help='calculate summary statistics from population comparison using allele frequencies') 
	target_stats_parser.add_argument('out', action='store', help='outfile prefix') 
	#could add option for block bootstrap
	bootstrap_parser = subparsers.add_parser('bootstrap', help='perform bootstrap estimates of population summary statistics in order to finalize model target values')
	bootstrap_parser.add_argument('n', action='store', type=int, help='number of bootstraps to perform in order to estimate standard error of the dataset (should converge for reasonably small n)')
	bootstrap_parser.add_argument('--in_freqs', action='store', help='comma-delimited list of infiles with per-site calculations for population. One file per population -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 
	bootstrap_parser.add_argument('--in_ld', action='store', help='comma-delimited list of infiles with per-site-pair calculations for population. One file per population -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 
	bootstrap_parser.add_argument('--in_fst', action='store', help='comma-delimited list of infiles with per-site calculations for population pair. One file per population-pair -- for bootstrap estimates of genome-wide values, should first concatenate per-chrom files)') 	
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
	print prefixstring + "calculating summary statistics for " +  str(npops) + " populations..."
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
		print prefixstring + command
	return
def execute_bootstrap(args):
	'''pulls all per-snp/per-snp-pair values to get genome-wide bootstrap estimates. adapted from JV experimental: get_neutral_targetstats_from_bootstrap.py'''
	nbootstraprep = args.n
	print prefixstring + "running " + str(nbootstraprep) + " bootstrap estimates of summary statistics..."
	targetstats_filename = args.out + "_bootstrap_n" + str(nbootstraprep) + ".txt"
	writefile = open(targetstats_filename, 'w')

	#################
	### FREQ STATS ##
	#################
	if args.in_freqs is not None: 
		inputestimatefilenames = ''.join(args.in_freqs)
		inputfilenames = inputestimatefilenames.split(',')
		npops = len(inputfilenames)
		for ipop in range(npops):
			inputfilename = inputfilenames[i]
			print prefixstring + "reading allele frequency statistics from: " + inputfilename
			writefile.write(str(ipop) + '\n')
			if checkFileExists(inputfilename):
				allpi, allnderiv, allnanc, nregions, seqlens = readFreqsFile(inputfilename)
			print prefixstring + "TOTAL: logged frequency values for " + str(nsnps) + " SNPS across " + str(totalregions) + ".\n"
			
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
		inputestimatefilenames = ''.join(args.in_ld)
		inputfilenames = inputestimatefilenames.split(',')
		npops = len(inputfilenames)
		for ipop in range(npops):
			inputfilename = inputfilenames[i]
			print prefixstring + "reading linkage disequilibrium statistics from: " + inputfilename
			writefile.write(str(ipop) + '\n')
			alldists, allr2, allgendists, alldprime, nr2regions, ndprimeregions = readLDFile(ldfilename, dprimecutoff = mafcutoffdprime)
			print prefixstring + "TOTAL: logged r2 values for " + str(allr2) + " SNP pairs.\n\tlogged D' values for " + str(alldprime) + " SNP pairs.\n"

			###################################
			### r2: MEAN ACROSS ALL REGIONS ###
			###################################
			r2sums, physDistHist = estimater2decay(allRegionr2, allRegionDists)
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
				repr2sum, repphysdisthist = estimater2decay(rep_all_r2, rep_all_physdist)
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
			compLDhist, genDistHist = estimatedprimedecay(allRegionDprime, allRegionGendists)
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

				repcompLDhist, repgenDistHist = estimatedprimedecay(rep_all_dprime, rep_all_gendist)
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
		for icomp in range(len(npopcomp)):
			fstfilename	= inputfilenames[icomp]
			print prefixstring + "reading Fst values from: " + fstfilename	
			if checkFileExists(fstfilename):
				allfst, nregions = readFstFile(fstfilename)
			target_mean, target_se = estimateFstByBootstrap_bysnp(allRegionValues, nrep = nbootstraprep)
			writeline =  str(icomp) + "\t" + str(target_mean) + "\t" + str(target_se) + '\n'
			writefile.write(writeline)
			print prefixstring + "TOTAL: logged Fst values for " + str(len(allRegionValues)) + " SNPs for "".\n"

	writefile.close()
	print prefixstring + "wrote to file: " + targetstats_filename
	return
def execute_point(args):
	'''runs simulates of a point in parameter-space, comparing to specified target. adapted from JV experimental: grid_point.py'''
	print prefixstring + "generating " + str(n) + " simulations from model: " + args.modelfile
	#calls to samplePoint from error.py
	#print prefixstring + "MUST CONNECT PIPELINE from cms_modeller.py TO grid_point.py"
	return
def execute_grid(args):
	'''run points in parameter-space according to specified grid'''
	print prefixstring + "loading dimensions of grid to search from: " + args.grid_inputdimensionsfile

	print prefixstring + "MUST CONNECT PIPELINE from cms_modeller.py TO grid.py"
	return
def execute_optimize(args):
	'''run scipy.optimize module according to specified parameters'''
	print prefixstring + "loading dimensions to search from: " + args.optimize_inputdimensionsfile
	print prefixstring + "MUST CONNECT PIPELINE from cms_modeller.py TO optimize.py (; optimize_out.py...)"
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