## {{DEM MODEL --> COMPONENT SCORE PROBABILITY DISTS}} -- script ii.
## last updated 05.31.16 	vitti@broadinstitute.org

############################################
##TOGGLE OPTIONS FOR THIS RUN OF PIPELINE ##
###### (i.e., across many runs of script) ##
scriptName = "dryrun_053116_debug" 			#subfolder for this pipeline run
numPerScenario, numNeut = 300, 500 			#for each selfreq bin in each pop in each demographic model
neutScenarios = ['default_default_101715_12pm']#, 'nulldefault', 'nulldefault_constantsize','gradient_101915_treebase_6_best', 'default_112115_825am', 'DEBUG_default_default_101715_12pm']
selPops, nThreads = [1, 2, 3, 4], 35

##########################################
##TOGGLE OPTIONS FOR THIS RUN OF SCRIPT ##
################# (i.e., w/in pipeline) ##

makeDirs = False
neutIhs, neutNormIhs, selIhs, selNormIhs = False, False, False, False
neutXp, neutNormXp, selXp, selNormXp = False, False, True, False
neutDelIhh,neutNormDelIhh, selDelIhh, selNormDelIhh  = False, False, False, False
neutFstDelDaf, selFstDelDaf, runPopPairs = False, False, [(4, 1), (4, 2), (4, 3), (3, 2), (3, 1), (2, 1)] 

##########
## MAIN ##
##########

from func_clean import *
import subprocess
import os

def main():

	############################
	### DEFINE BINS OF DISTS ##
	###########################
	if True:
		fullrange, binlen = [.05, .95], .1 #define selfreq bins
		numbins = int((fullrange[1] - fullrange[0])/binlen)
		bin_starts = [fullrange[0] + binlen * i for i in range(numbins+1)]
		bin_ends = [fullrange[0] + binlen * i for i in range(1, numbins+2)]
		bin_medians = [fullrange[0] + (binlen/2) + i*.1 for i in range(numbins+1)]
		bin_medians_str = [str(x)[:3] for x in bin_medians]

	######################################
	## MAKE DIRS FOR SCORES BY SCENARIO ##
	######################################
	tpedloc = 'sims/cosi_tped/'
	ihsloc, delihhloc, xploc, fstloc = "scores/cosi_ihs/", "scores/cosi_delihh", "scores/cosi_xp", 'scores/cosi_fst_deldaf'
	parentdirs = [ihsloc, delihhloc, xploc, fstloc]
	verboseDirs = False
	if makeDirs == True:
		for parentdir in parentdirs:
			checkMakeDir(parentdir, verbose = verboseDirs)
			for neutScenario in neutScenarios:
				scenariodir = parentdir + neutScenario
				checkMakeDir(scenariodir, verbose = verboseDirs)
				neutscenariodir = scenariodir + "/neut"
				checkMakeDir(neutscenariodir, verbose = verboseDirs)
				
				for selPop in selPops:
					selpopdir = parentdir +  neutScenario + "/sel" + str(selPop)
					checkMakeDir(selpopdir, verbose = verboseDirs)
					for selfreq in bin_medians_str:
						selfreqdir =  parentdir +  neutScenario + "/sel" + str(selPop) +  "/sel" + str(selfreq)
						checkMakeDir(selfreqdir, verbose = verboseDirs)

						terminalfolder = selfreqdir + "/" + scriptName
						checkMakeDir(terminalfolder, verbose = verboseDirs)

	##############
	## RUN IHS ###
	##############
	if neutIhs == True:
		for iscenario in range(len(neutScenarios)):
			tpeddir = tpedloc + neutScenarios[iscenario] + "/neut/"
			ihsdir = ihsloc  + neutScenarios[iscenario] + "/neut/"
			inputs, outputs = [], []
			for iTped in range(numPerScenario):
				for selPop in selPops:
					inputfilename = tpeddir + "rep" + str(iTped) +"_pop" + str(selPop) +".tped"
					#print inputfilename
					#assert os.path.isfile(inputfilename)
					if os.path.isfile(inputfilename):
						#print inputfilename
						outputfilename = ihsdir + "rep" + str(iTped) +"_pop" + str(selPop)
						inputs.append(inputfilename)
						outputs.append(outputfilename)
			ihs_scenarioname = "neut_ihs_" + neutScenarios[iscenario] + "_" + scriptName 
			launch_calc_ihs_batch(ihs_scenarioname, inputs, outputs, numThreads = nThreads)
	if neutNormIhs == True:
		for iscenario in range(len(neutScenarios)):
			ihsdir = ihsloc  + neutScenarios[iscenario] + "/neut/"
			print "normalizing neutral scores from " + ihsdir 
			normalizeiHS(neutScenarios[iscenario], numPerScenario, ihsloc =  ihsdir, pops = [1, 2, 3, 4])
	if selIhs == True:
		for iscenario in range(len(neutScenarios)):
			for selPop in selPops:
				for selfreq in bin_medians_str:
					tpeddir = tpedloc +  neutScenarios[iscenario] + "/sel" + str(selPop) +  "/sel" + str(selfreq) + "/" + scriptName + "/"
					ihsdir =  ihsloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/" 
					inputs, outputs = [], []
					for iTped in range(numPerScenario):
						inputfilename = tpeddir + "rep" + str(iTped) +"_pop" + str(selPop) +".tped"
						#print inputfilename
						#assert os.path.isfile(inputfilename)
						#if os.path.isfile(inputfilename):
							#print inputfilename

						if os.path.isfile(inputfilename):
							outputfilename = ihsdir + "rep" + str(iTped) +"_pop" + str(selPop)
							inputs.append(inputfilename)
							outputs.append(outputfilename)
					ihs_scenarioname = "ihs_" + neutScenarios[iscenario] + "_sel" + str(selPop) + "_" + str(selfreq) + "_" + scriptName 
					launch_calc_ihs_batch(ihs_scenarioname, inputs, outputs, numThreads = nThreads)
	if selNormIhs == True:
		neutihsparamsloc = "neut_norm_params/ihs/"
		entries = os.listdir(neutihsparamsloc)

		for iscenario in range(len(neutScenarios)):
			outfileprefix = "norm_ihs_" + neutScenarios[iscenario] + ".sh."
			model_neutparamfiles = [item for item in entries if outfileprefix in item]
			#print model_neutparamfiles
			assert len(model_neutparamfiles) == 4
			model_neutparamfiles.sort()
			for i in range(len(selPops)):
				selPop = selPops[i]
				#get neutral parameters to normalize
				norm_param_filename = neutihsparamsloc + model_neutparamfiles[i]
				print "loading neutral parameters from: " + norm_param_filename
				bins, nums, means, variances = readNeutNormFile(norm_param_filename, scoretype='ihs')
				print "n neut snp: " + str(sum(nums))

				for selfreq in bin_medians_str:
					ihsdir = ihsloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					
					#check to make sure we have the right number of replicates
					ihsfiles = [item for item in os.listdir(ihsdir) if ".out" in item]
					#print ihsfiles
					if len(ihsfiles) != numPerScenario:
						print "found " + str(len(ihsfiles)) + " ihs files"
					normalize_Ihs_manually(ihsdir, numPerScenario, selPop, means, variances, bins)			

	#################
	## RUN delIHH ###
	#################
	if neutDelIhh == True:
		for iscenario in range(len(neutScenarios)):
			ihsdir = ihsloc + neutScenarios[iscenario] + "/neut/"
			delihhdir = delihhloc + neutScenarios[iscenario] + "/neut/"
			for selPop in selPops:
				inputfilename = ihsdir + "pop" + str(selPop) + ".ihs.concat"
				if not os.path.isfile(inputfilename):
					print inputfilename
				outpufilename = delihhdir + "pop" + str(selPop) + ".delihh.concat"
				write_delIHH_file(inputfilename, outpufilename)
				#launch_calcdelihh(scenario, ihsloc, delihhloc)
	if neutNormDelIhh == True:
		for iscenario in range(len(neutScenarios)):
			delihhdir = delihhloc + neutScenarios[iscenario] + "/neut/"
			normalize_delIhh(neutScenarios[iscenario], delihhdir)
	if selDelIhh == True:
		for iscenario in range(len(neutScenarios)):
			for selPop in selPops:
				for selfreq in bin_medians_str:
					ihsdir = ihsloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					inputfilename = ihsdir + "n300.ihs.concat"
					if not os.path.isfile(inputfilename):
						print inputfilename
					delihhdir = delihhloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					outpufilename = delihhdir + "n300.delihh.concat"
					write_delIHH_file(inputfilename, outpufilename)
					#launch_calcdelihh(scenario, ihsloc, delihhloc)
	if selNormDelIhh == True:
		neutdelihhparamsloc = "/scores/neut_norm_params/delihh/"
		entries = os.listdir(neutdelihhparamsloc)
		for iscenario in range(len(neutScenarios)):
			outfileprefix = "norm_delihh_" + neutScenarios[iscenario] + ".sh."
			model_neutparamfiles = [item for item in entries if outfileprefix in item]
			#print model_neutparamfiles
			assert len(model_neutparamfiles) == 4
			model_neutparamfiles.sort()
			for i in range(len(selPops)):
				selPop = selPops[i]
				norm_param_filename = neutdelihhparamsloc + model_neutparamfiles[i]
				print "loading neutral parameters from: " + norm_param_filename
				bins, nums, means, variances = readNeutNormFile(norm_param_filename, scoretype='ihs')
				print "n neut snp: " + str(sum(nums))
				for selfreq in bin_medians_str:
					delihhdir = "/scores/cosi_delihh/" + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					delihhfilename = delihhdir + "n300.delihh.concat"
					assert os.path.isfile(delihhfilename)
					normalize_delihh_manually(delihhfilename, means, variances, bins)			

	#############
	## RUN XP ###
	#############
	if neutXp == True:
		for iscenario in range(len(neutScenarios)):
			tpeddir = tpedloc + neutScenarios[iscenario] + "/neut/"

			#tpeddir = tpedloc + neutScenarios[iscenario] + "/neut/" 
			xpdir = xploc + neutScenarios[iscenario] + "/neut/"
			inputs_1, inputs_2, outputs = [], [], []	
			for iTped in range(numPerScenario):
				for selPop in selPops:
					inputfilename1 = tpeddir + "rep" + str(iTped) +"_pop" + str(selPop) +".tped"
					if os.path.isfile(inputfilename1):
					#assert os.path.isfile(inputfilename1)
						for altpop in selPops:
			 				if altpop != selPop:
								inputfilename2 = tpeddir + "rep" + str(iTped) + "_pop" + str(altpop) + ".tped"
								outputfilename = xpdir + "rep" + str(iTped) +"_pop" + str(selPop) +"_pop" +str(altpop)
								inputs_1.append(inputfilename1)
								inputs_2.append(inputfilename2)
								outputs.append(outputfilename)
			xp_scenarioname = "xp_" + neutScenarios[iscenario] + "_neut"
			launch_calc_xp_batch(xp_scenarioname, inputs_1, inputs_2, outputs, numThreads = nThreads)
	if neutNormXp == True:
		popPairs = []
		for pop in selPops:
			for altpop in selPops:
				if altpop != pop:
					popPairs.append((pop, altpop))
		for iscenario in range(len(neutScenarios)):
			xpdir = xploc+ neutScenarios[iscenario] + "/neut/"
			normalizeXpehh(neutScenarios[iscenario], numPerScenario, Xploc = xpdir, popPairs = popPairs)
	if selXp == True:
		for iscenario in range(len(neutScenarios)):
			for selPop in selPops:
				for selfreq in bin_medians_str:
					tpeddir = tpedloc +  neutScenarios[iscenario] + "/sel" + str(selPop) +  "/sel" + str(selfreq) + "/" + scriptName + "/"

					#tpeddir = tpedloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"+ scriptName +"n" + str(numPerScenario) + "/"
					xpdir = xploc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					inputs_1, inputs_2, outputs = [], [], []
					for iTped in range(numPerScenario):
						inputfilename1 = tpeddir + "rep" + str(iTped) +"_pop" + str(selPop) +".tped"
						if os.path.isfile(inputfilename1):
						#assert os.path.isfile(inputfilename1)
							for altpop in selPops:
		 						if altpop != selPop:
									inputfilename2 = tpeddir + "rep" + str(iTped) + "_pop" + str(altpop) + ".tped"
									outputfilename = xpdir + "rep" + str(iTped) +"_pop" + str(selPop) +"_pop" +str(altpop)
									inputs_1.append(inputfilename1)
									inputs_2.append(inputfilename2)
									outputs.append(outputfilename)
					xp_scenarioname = "xp_" + neutScenarios[iscenario] + "_sel" + str(selPop) + "_" + str(selfreq) + "_" + scriptName + "n" + str(numPerScenario)
					launch_calc_xp_batch(xp_scenarioname, inputs_1, inputs_2, outputs, numThreads = nThreads)
	if selNormXp == True:
		neutxpparamsloc = "/scores/neut_norm_params/xp/"
		entries = os.listdir(neutxpparamsloc)
		for iscenario in range(len(neutScenarios)):
			outfileprefix = "norm_xp_" + neutScenarios[iscenario] + ".sh."
			model_neutparamfiles = [item for item in entries if outfileprefix in item]
			#print model_neutparamfiles
			assert len(model_neutparamfiles) == 12
			model_neutparamfiles.sort()
			iFile = 0
			for i in range(len(selPops)):
				selPop = selPops[i]
				for altPop in selPops:
					if altPop != selPop:
						print "\n\t selpop: " + str(selPop) + " ; altpop: " + str(altPop) + "\n"
						#get neutral parameters to normalize

						norm_param_filename = neutxpparamsloc + model_neutparamfiles[iFile]
						print "loading neutral parameters from: " + norm_param_filename
						num, mean, var = readNeutNormFile(norm_param_filename, scoretype='xp')
						print "n neut snp: " + str(num)

						for selfreq in bin_medians_str:
							xpdir = xploc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"	
							normalize_Xp_manually(xpdir, numPerScenario, mean, var, selPop, altPop)
						iFile +=1 

	#####################	
	## RUN FST/DELDAF ###
	#####################
	if neutFstDelDaf== True:
		for iscenario in range(len(neutScenarios)):
			tpeddir = tpedloc + neutScenarios[iscenario] + "/neut/" 
			fstdir = fstloc + neutScenarios[iscenario] + "/neut/"
			launch_fst_deldaf(neutScenarios[iscenario], 300, tpeddir, fstdir, popPairs = runPopPairs)
	if selFstDelDaf == True:
		for iscenario in range(len(neutScenarios)):
			tpedbase = tpedloc + neutScenarios[iscenario] + "/" 
			fstbase = fstloc + neutScenarios[iscenario] + "/"
			for i in range(len(selPops)):
				selPop = selPops[i]
				tpeddir = tpedbase + "sel" + str(selPop) + "/"
				fstdir =  fstbase + "sel" + str(selPop) + "/"
				for selfreq in bin_medians_str:
					fulldir_tped = tpeddir + "sel" + str(selfreq) + "/sweep_full_032516n300/" #for now
					fulldir_fst = fstdir + "sel" + str(selfreq) + "/"
					fullscenarioname = neutScenarios[iscenario] + "_sel" + str(selPop) + "_" + str(selfreq)
					launch_fst_deldaf(fullscenarioname, 300, fulldir_tped, fulldir_fst, popPairs = runPopPairs)

main()