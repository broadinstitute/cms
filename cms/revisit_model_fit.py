##	quantify goodness-of-fit to empirical for each demographic model 
##	12.02.2017

models = ['default_112115_825am', 'default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault', 'nulldefault_constantsize']
paramdir = "/idi/sabeti-scratch/jvitti/params/"
outputdir = "/idi/sabeti-scratch/jvitti/1kg_models/"
numCosiReps = 500
calcTarget, calcStats, calcError, plotStats = False, False, False, True
pops = ['YRI', 'CEU', 'CHB', 'BEB']
chroms = range(1,23)

import matplotlib as mp
mp.use('agg')

import os
import time
import subprocess
from jv_uger_util import *
import importlib.machinery
error_func = importlib.machinery.SourceFileLoader('error_func','/idi/sabeti-scratch/jvitti/cms/cms/model/error_func.py').load_module()
from error_func import calc_error

plot_func = importlib.machinery.SourceFileLoader('plot_func','/idi/sabeti-scratch/jvitti/cms/cms/model/plot_func.py').load_module()
from plot_func import plot_comparison


def main():
	###############################
	## RECALCULATE TARGET VALUES ##
	###############################
	ncmds_script, basecmd, dispatchNow = 10, "", True
	scriptbase = "/idi/sabeti-scratch/jvitti/scripts/revisit_model_target"
	if calcTarget:
		allcmds = []
		for pop in pops:
			for chrom in chroms:
				in_tped_filename = "/idi/sabeti-scratch/jvitti/1kg_tpeds/unmasked_aa_v2/full_1kg_chr" + str(chrom) + "_" + pop + ".tped.gz"
				recomfilename = "/idi/sabeti-scratch/ilya/gsvn/Data/Ilya_Data/genmaps/hm2/genetic_map_chr" + str(chrom) + "_b37.txt"#GRCh37_chr" + str(chrom) + ".txt"
				regionfilename = "/idi/sabeti-scratch/jvitti/1kg_models/NRE_regions/chr" + str(chrom) + "_putative_neutral_regions.bed"#test_regions.txt"
				outfilename = outputdir + 'target_stats_' + str(pop) + "_chr" + str(chrom) 
				#print('gunzip ' + in_tped_filename)
				#test_cmd = "/idi/sabeti-scratch/jvitti/cms/cms/model/bootstrap_freq_popstats_regions " + in_tped_filename.strip('.gz') + " " + recomfilename + " " + regionfilename + " " + outfilename
				#print(test_cmd)
				freq_cmd = "/idi/sabeti-scratch/jvitti/cms/cms/model/bootstrap_freq_popstats_regions " + in_tped_filename.strip('.gz') + " " + recomfilename + " " + regionfilename + " " + outfilename +"_freqs"
				ld_cmd = "/idi/sabeti-scratch/jvitti/cms/cms/model/bootstrap_ld_popstats_regions " + in_tped_filename.strip('.gz') + " " + recomfilename + " " + regionfilename + " " + outfilename+"_ld"
				#if not os.path.isfile(outfilename +"_freqs"):
				#	allcmds.append(freq_cmd)
				if not os.path.isfile(outfilename +"_ld"):
					allcmds.append(ld_cmd)
				altpops = pops[:]
				altpops.remove(pop)
				#for altpop in altpops:
				#	in_tped2_filename = "/idi/sabeti-scratch/jvitti/1kg_tpeds/unmasked_aa_v2/full_1kg_chr" + str(chrom) + "_" + altpop + ".tped.gz"			
				#	fst_cmd = "/idi/sabeti-scratch/jvitti/cms/cms/model/bootstrap_fst_popstats_regions " + in_tped_filename.strip('.gz') + " "+ in_tped2_filename.strip('.gz') + " " + recomfilename + " " + regionfilename + " " + outfilename+ "_" + altpop + "_fst"
				#	allcmds.append(fst_cmd)
		print('loaded a total of ' + str(len(allcmds)) + " commands.")
		iscript = 0
		for argchunk in chunks(allcmds,ncmds_script):
			iscript +=1
			scriptname = scriptbase + "_" + str(iscript)		
			slurmcmd = writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, UGERoptions=['-l h_vmem=10g', '-cwd', '-P sabeti_lab', '-l h_rt=48:00:00'])
			if dispatchNow:
				subprocess.check_output( slurmcmd.split() )
			time.sleep(1)		
	print("run bootstrap, then put together with get_targetstats_frombootstrap.py")

	#######################
	## CALC CUSTOM STATS ##
	#######################
	if calcStats:
		for model in models:
			paramfilename = paramdir  + model + ".par"
			statfilename = outputdir + model + "_stats_n" + str(numCosiReps) + ".txt"
			fullcmd =  "coalescent -p " + paramfilename + " --genmapRandomRegions --drop-singletons .25 "#--tped " + writefilename + " --output-gen-map"
			fullcmd += " --output-gen-map -n " + str(numCosiReps)
			fullcmd += " --custom-stats > " + statfilename
			#fullcmd += " --stop-after-minutes " + str(args.stopAfterMinutes)
			print(fullcmd)

	#################
	## CALC ERROR ###
	#################
	if calcError:
		for model in models:
			statfilename = outputdir + model + "_stats_n" + str(numCosiReps) + ".txt"
			error = calc_error(statfilename, verbose=True)
			print(error, statfilename)

	if plotStats:
		for model in models:
			statfilename = outputdir + model + "_stats_n" + str(numCosiReps) + ".txt"
		
			plot_comparison(statfilename, 2, "/web/personal/vitti/revisit_" + model)


main()

'''
#extant neutral(/sim) files galore...
if model == "gradient_101915_treebase_6_best":
tped_dir = "/idi/sabeti-scratch/jvitti/standing_5mb/" + model + "/tpeds/"
elif model == "default_112115_825am"
tped_dir = "/idi/sabeti-scratch/jvitti/standing_3mb/" + model + "/tpeds/"
else:
tped_dir = "/idi/sabeti-scratch/jvitti/standing_1mb/" + model + "/tpeds/"
'''

