## top-level script for passing varying input parameters to power.py (et al) in parallel via UGER.
## last updated 11.29.16

checkOverwrite = True

################
### TOGGLE RUN #
################
rundir = "/idi/sabeti-scratch/jvitti/cms/cms/"
runstring ="sel_poppair"
runNeutSims, runNeutRepscores, normNeut_repscores, normNeut_frombins = False, False, False, False
runPoppair = True
compositeSims, normSims = False, False
runLikes, runLikesComp = False, False
compositeEmp, normEmp = False, False
repViz, distViz = False, False
runCDF, runFPR, runTPR, runROC = False, False, False, False
runGwRegions, runRegionLog, runManhattan, runExtendedManhattan = False, False, False, False
cmdPrefix = "python /idi/sabeti-scratch/jvitti/power.py "

alt_cmdPrefix = "python /idi/sabeti-scratch/jvitti/cms/cms/"

baseDirs = {'sims':'/idi/sabeti-scratch/jvitti/clean/sims/', 'scores':'/idi/sabeti-scratch/jvitti/clean/scores/', 'plot':"/web/personal/vitti/", 'sel_basedir':'/idi/sabeti-scratch/jvitti/scores/'}

models = [ 'nulldefault', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm', 'default_112115_825am','nulldefault_constantsize',]
likes_suffices = ["neut"]#, "linked"]
regionlens = [10000, 25000, 50000, 100000]
thressholds = [25, 30, 35, 40, 45, 50]
cutoffs = [2, 2.5, 3, 3.5, 4, 4.5, 5]
modelpops = [1, 2, 3, 4]
numNeutReps, numSelReps = 1000, 500
empPops = ['CEU', 'YRI', 'CHB', 'BEB']
selbin_clusters = [["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]] #["0.10", "0.20", "0.30"], ["0.40", "0.50", "0.60",], ["0.70", "0.80", "0.90",],
cluster_names = [ "allfreq"]#"low", "mid", "hi",]
scoredirs =['delihh', 'fst_deldaf', 'ihs', 'nsl', 'xpehh']
##########
## MAIN ##
##########

from parse_func import writeJobsToTaskArray_fromArglist, check_create_dir
from parse_func import getSelFiles_asList, getNeutFiles_asList, getNeutCompFiles, getSelCompFiles #REPLACE THESE?
import sys
import os

def main():
	allOpts = [runNeutSims, runNeutRepscores, normNeut_repscores, normNeut_frombins, runPoppair, compositeSims, normSims, runLikes, runLikesComp, 
		compositeEmp, normEmp, repViz, distViz, runCDF, runFPR, runTPR, runROC, runGwRegions, runRegionLog, runManhattan, runExtendedManhattan]
	if allOpts.count(True) != 1:
		print("call one command to dispatch at a time.")
		print(allOpts)
		sys.exit(0)


	args = []
	discardOut = False
	mem = "1g"

	########	Manipulate simulated data 
	########	and provide composite scores.
	if runNeutSims:
		cmd = cmdPrefix + "run_neut_sims"
		basedir = baseDirs['sims']
		check_create_dir(basedir)
		for model in models:
			modeldir = basedir + model + "/"
			check_create_dir(modeldir)
			neutdir = modeldir + "neut/"
			check_create_dir(neutdir)
			argstring = "--model " + model + " --writedir " + neutdir + " --nrep " + str(numNeutReps)
			args.append(argstring)
	if runNeutRepscores:
		cmd = cmdPrefix + "run_neut_repscores"
		basedir = baseDirs['scores']
		check_create_dir(basedir)
		for model in models:
			modeldir = basedir + model
			check_create_dir(modeldir)
			neutmodeldir = modeldir + "/neut/"
			check_create_dir(neutmodeldir)

			tpeddir = baseDirs['sims'] + model + "/neut/"

			for score in scoredirs:
				modelscoredir = neutmodeldir + score + "/"
				check_create_dir(modelscoredir)

			for modelpop in modelpops:
				for irep in range(1, numNeutReps +1):
					argstring = "--model " + model + " --simpop " + str(modelpop) + " --irep " + str(irep) + " --writedir " + neutmodeldir + " --tpedfolder " + tpeddir + " --checkOverwrite"
					args.append(argstring)
	if normNeut_repscores:
		cmd = cmdPrefix + "run_norm_neut_repscores"
		basedir = baseDirs['scores']
		check_create_dir(basedir)
		for model in models:
			modeldir = basedir + model
			check_create_dir(modeldir)
			neutmodeldir = modeldir + "/neut/"
			check_create_dir(neutmodeldir)
			for modelpop in modelpops:
				for score in ['delihh']:#['ihs', 'delihh', 'nsl',]:
					argstring = "--model " + model + " --simpop " + str(modelpop) + " --writedir " + basedir + " --nrep " + str(numNeutReps) + " --score " + score
					args.append(argstring)
				#for score in ['fst', 'xpehh']:
				#	altpops = modelpops[:]
				#	altpops.remove(modelpop)
				#	for altpop in altpops:
				#		argstring = "--model " + model + " --simpop " + str(modelpop) + " --writedir " + basedir + " --nrep " + str(numNeutReps) + " --score " + score + " --altpop " + str(altpop)
				#		args.append(argstring)
	if normNeut_frombins:
		cmd = cmdPrefix + "norm_from_binfile"
		basedir = baseDirs['scores']
		check_create_dir(basedir)
		for model in models:
			modeldir = basedir + model
			check_create_dir(modeldir)
			neutmodeldir = modeldir + "/neut/"
			check_create_dir(neutmodeldir)
			for modelpop in modelpops:
				for score in ['delihh']:#['ihs', 'delihh', 'nsl',]:
					argstring = "--model " + model + " --simpop " + str(modelpop) + " --writedir " + basedir + " --nrep " + str(numNeutReps) + " --score " + score
					args.append(argstring)
				#for score in ['fst', 'xpehh']:
				#	altpops = modelpops[:]
				#	altpops.remove(modelpop)
				#	for altpop in altpops:
				#		argstring = "--model " + model + " --simpop " + str(modelpop) + " --writedir " + basedir + " --nrep " + str(numNeutReps) + " --score " + score + " --altpop " + str(altpop)
				#		args.append(argstring)
	if runPoppair:
		basedir = baseDirs['scores']
		cmd = cmdPrefix + "run_poppair"
		for model in models:
			neutmodeldir = basedir + model + "/neut/"
			pairdir = neutmodeldir + "pairs/"
			check_create_dir(pairdir)
			for modelpop in modelpops:
				altpops = modelpops[:]
				altpops.remove(modelpop)
				for altpop in altpops:
					sel_basedir = baseDirs['sel_basedir']
					argstring = "--model " + model + " --simpop " + str(modelpop) + " --altpop " + str(altpop) + " --writedir " + sel_basedir + " --nrep " + str(numSelReps)
				#	argstring = "--model " + model + " --simpop " + str(modelpop) + " --altpop " + str(altpop) + " --writedir " + basedir + " --nrep " + str(numNeutReps) #NEUT
					args.append(argstring)
	if compositeSims:
		cmd = cmdPrefix + "composite_sims" 
		for model in models:
			writedir = baseDirs['scores'] + model + "/neut/composite/"
			check_create_dir(writedir)

			for modelpop in modelpops:
				argstring = "--model " + model + " --simpop " + str(modelpop) + " --writedir " + writedir
				args.append(argstring)
	if normSims:
		cmd = cmdPrefix + "normsims"
		for model in models:
			for modelpop in modelpops:
				argstring = "--model " + model + " --simpop " + str(modelpop) + " --nrep 1000"
				args.append(argstring)

	########	Define likelihood functions
	########	from simulated data. (points to likes_from_model.py)
	if runLikes:
		cmd = alt_cmdPrefix + "likes_from_model.py likes_from_scores"
		mem = "3g"
		scores = ['ihs','nsl', 'delihh',]
		for model in models:
			model_write_dir = writedir + model + "/"
			check_create_dir(model_write_dir)
			for score in scores:
				
				print("finish patch")
				sys.exit()
				scorelistdir = ""
				writedir = ""

				model_write_likes_dir = model_write_dir + score + "/"
				check_create_dir(model_write_likescorelistdirs_dir)
				for pop in pops:
					neutlistfilename = scorelistdir + model +"/" + score + "/repfiles_neut_" + str(pop) + ".list"	
					neutfilename = getNeutFiles_asList(model, score, pop, neutlistfilename, overwrite=overwriteList, numPerBin=numPerBin)
					for icluster in range(len(cluster_names)):
						selbins = selbin_clusters[icluster]
						binclustername = cluster_names[icluster]
						listfilename = scorelistdir +model +"/" + score + "/repfiles_sel" + str(pop) + "_" + binclustername + ".list"	
						selfilename = getSelFiles_asList(model, score, pop, listfilename, selbins, binclustername, overwrite=False)

						outprefix = model_write_likes_dir + "likes_sel" + str(pop) + "_" + binclustername 		

						outfilenames = [outprefix +"_causal.txt", outprefix + "_neut.txt"]
						if not (os.path.isfile(outfilenames[0]) and os.path.isfile(outfilenames[1])):
						#if True:
							argstring = "--" + score + " " + neutfilename + " " + selfilename + " " + str(selPos) +  " " + outprefix + " --nLikesBins " + str(numLikesBins)
							arguments.append(argstring)
	if runLikesComp:
		cmd = alt_cmdPrefix + "likes_from_model.py likes_from_scores"
		scores = ['deldaf', 'fst',  'xpehh']
		mem = "10g"
		check_create_dir(writedir)
		for model in models:
			model_likes_dir = scorelistdir + model + "/"
			model_write_dir = writedir + model + "/"
			check_create_dir(model_likes_dir)	
			check_create_dir(model_write_dir)
			for score in scores:
				model_score_likes_dir = model_likes_dir + score + "/"
				check_create_dir(model_score_likes_dir)
				write_score_dir = model_write_dir + score + "/"
				check_create_dir(write_score_dir)
					
				neutfilename = getNeutCompFiles(model, pop, model_score_likes_dir, runPrefix, overwrite=overwriteList)

				for icluster in range(len(cluster_names)):
					selbins = selbin_clusters[icluster]
					binclustername = cluster_names[icluster]

					selfilename = getSelCompFiles(model, pop, model_score_likes_dir, runPrefix, selbins, binclustername, overwrite=overwriteList)
				
					outprefix = write_score_dir  + "likes_sel" + str(pop) + "_choose_" + binclustername

					#CHECK IF COMPLETED
					#if True:
					if not os.path.isfile(outprefix + "_causal.txt") or not os.path.isfile(outprefix + "_neut.txt"):
						argstring = "--" + score + " " + neutfilename + " " + selfilename + " " + str(selPos) +  " " + outprefix + " --nLikesBins " + str(numLikesBins)
						arguments.append(argstring)

	########	Manipulate empirical data
	########	and provide composite scores
	if compositeEmp:
		cmd = cmdPrefix + "composite_emp"
		for model in model:
			for emppop in empPops:
				argstring = "--model " + model + " --emppop " + str(emppop)
				args.append(argstring)
	if normEmp:
		cmd = cmdPrefix + "normemp"
		for emppop in empPops:
			for model in models:
				outfilename = "/"
				argstring = "--model " + model + " --emppop " + emppop 

				args.append(argstring)

	########	Visualize composite score 
	########	output for a given CMS run.
	if repViz:

		vizReps = range(1,20)#
		plotBase = baseDirs['plot']  + "repviz/"
		cmd = cmdPrefix + "repviz"
		for model in models:
			modeldir = plotBase + model
			check_create_dir(modeldir)
			for modelpop in modelpops:
				for irep in vizReps:
					savefilename = modeldir + "/rep" + str(irep) + "_" + str(modelpop) + ".png"
					argstring = "--model " + model + " --simpop " + str(modelpop) + " --savefilename " + savefilename + " --vizRep " + str(irep)
					args.append(argstring)
	if distViz:
		plotBase = baseDirs['plot']
		cmd = cmdPrefix + "distviz --simsDist"
		for model in models:
			for simpop in modelpops:
				savefilename = plotBase + model +"_" + str(simpop) + ".png"
				argstring = "--simpop " + str(simpop) + " --model " + model + " --savefilename " + savefilename
				args.append(argstring)

	########	Quantify and visualize power
	########	across significance cutoffs.
	if runCDF:
		cmd = cmdPrefix + "cdf"
		for model in models:
			argstring = "--model " + model
	if runFPR:
		savebasedir = "/idi/sabeti-scratch/jvitti/clean/fpr/"
		cmd = cmdPrefix +'fpr'
		for model in models:
			for likessuffix in likes_suffices:
				for regionlen in regionlens:
					for thresshold in thressholds:
						for cutoff in cutoffs:
							for modelpop in modelpops:
								savefilename = savebasedir + runstring + "_" + model + "_" + str(modelpop) + "_" + likessuffix + "_" + str(regionlen) +"_" + str(thresshold) +"_" + str(cutoff) + ".txt"
								argstring = str(regionlen) + " " + str(thresshold) + " " + str(cutoff) + " --likessuffix " + likessuffix + " --model " + model + " --saveLog " + savefilename + " --simpop " + str(modelpop) + " --nrep 1000"
								#if True:
								if not os.path.isfile(savefilename):
									args.append(argstring)
	if runTPR:
		savebasedir = "/idi/sabeti-scratch/jvitti/cms2_power/tpr_4c/"
		cmd = cmdPrefix + 'tpr'
		for model in models:
			for likessuffix in likes_suffices:
				for regionlen in regionlens:
					for thresshold in thressholds:
						for cutoff in cutoffs:
							for modelpop in modelpops:
								savefilename = savebasedir + runstring + "_" + model + "_" + str(modelpop) + "_" + likessuffix + "_" + str(regionlen) +"_" + str(thresshold) +"_" + str(cutoff) + ".txt"
								argstring = str(regionlen) + " " + str(thresshold) + " " + str(cutoff) + " --likessuffix " + likessuffix + " --model " + model + " --saveLog " + savefilename + " --simpop " + str(modelpop)
								if not os.path.isfile(savefilename):
									args.append(argstring)

	########	Apply significance cutoff 
	########	to empirical results.
	if runGwRegions:
		cmd = cmdPrefix + "gw_regions"
		savebasedir = "/idi/sabeti-scratch/jvitti/cms2_power/regions/"
		for model in models:
			for empPop in empPops:
				for likessuffix in likes_suffices:
					for regionlen in regionlens:
						for thresshold in thressholds:
							for cutoff in cutoffs:
								savefilename = savebasedir + runstring + "_" + empPop +"_" + model +"_" + likessuffix + "_" + str(regionlen) +"_" + str(thresshold) +"_" + str(cutoff) + ".txt"
								argstring = str(regionlen) + " " + str(thresshold) + " " + str(cutoff) + " --emppop " + empPop +  " --likessuffix " + likessuffix + " --model " + model + " --saveLog " + savefilename 
								args.append(argstring)
		#pass
	if runManhattan:
		savebasedir  ='/web/personal/vitti/manhattan/'
		cmd = cmdPrefix + "manhattan"
		for model in models:
			for empPop in empPops:
				savefilename = savebasedir + empPop + "_" + model + "_" + runstring + '.png' 
				argstring = '--emppop ' + empPop + " --model " + model + " --zscore --savefilename " + savefilename
				args.append(argstring)
		pass
	if runExtendedManhattan:
		savebasedir  ='/web/personal/vitti/manhattan/'
		cmd = cmdPrefix + "extended_manhattan"
		for model in models:
			for empPop in empPops:
				for likessuffix in likes_suffices:
					for regionlen in regionlens:
						for thresshold in thressholds:
							for cutoff in cutoffs:
								regionfilename = "/idi/sabeti-scratch/jvitti/cms2_power/regions/regions_111716_" + empPop +"_" + model +"_" + likessuffix + "_" + str(regionlen) +"_" + str(thresshold) +"_" + str(cutoff) + ".txt"
								savefilename = savebasedir + empPop + "_" + model + "_" + str(regionlen) +"_" + str(thresshold) +"_" + str(cutoff)  + '.png' 
								argstring = '--emppop ' + empPop + " --model " + model + " --savefilename " + savefilename + " --regionsfile " + regionfilename
								args.append(argstring)
	
	writeJobsToTaskArray_fromArglist(cmd, args, rundir + runstring, discardOut=discardOut, memory=mem, UGERoptions = ['-P sabeti_lab', '-wd ' + rundir])	

main()

