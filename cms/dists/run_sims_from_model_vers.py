## {{DEM MODEL --> COMPONENT SCORE PROBABILITY DISTS}} -- script i.
## last updated 06.01.16 	vitti@broadinstitute.org

############################################
##TOGGLE OPTIONS FOR THIS RUN OF PIPELINE ##
###### (i.e., across many runs of script) ##

scriptName = "dryrun_053116_debug" 			#subfolder for this pipeline run
numPerScenario, numNeut = 300, 500 			#for each selfreq bin in each pop in each demographic model
paramsloc, neutScenarios = "sims/params/", ['default_default_101715_12pm']#, 'DEBUG_default_default_101715_12pm']#, 'nulldefault', 'nulldefault_constantsize','gradient_101915_treebase_6_best', 'default_112115_825am']
cosi_build = "/idi/sabeti-scratch/ilya/cosi/builds/bnew52opt/coalescent_may25" #"/idi/sabeti-scratch/jvitti/synth/likes/debug_052316/cosi-2.0/coalescent"
singRate, sweepPos = .25, .5
selPops, selCoeffs = [1,2,3,4], {1:"E(20)", 2:"E(20)", 3:"E(20)", 4:"E(20)"} 

##########################################
##TOGGLE OPTIONS FOR THIS RUN OF SCRIPT ##
################# (i.e., w/in pipeline) ##

makeDirs = False   											 		
generateTrajs, nTrajPerScenario, presentDayFreq = False, 750, ".05-.95"		
checkSelFreqs = False			
runCosiMS, scaleSims, simScaleFactor = False, False, 10
cleanUpMS = False
countMS, trimMS, maxToKeep = False, False, 1000 #optional
convertTped = False											
runNeut, cleanNeut, convertNeutTped = False, False, True 									

##########
## MAIN ##
##########

from func_clean import *
from random import choice
import os
import subprocess

def main():
	cwd = os.getcwd()

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

	####################################
	## MAKE DIRS FOR SIMS BY SCENARIO ##
	####################################
	checkMakeDir('sims/')
	trajloc, msloc, tpedloc =  'sims/cosi_seltrajs/', "sims/cosi_ms/", 'sims/cosi_tped/' #will create
	parentdirs = [trajloc,  msloc, tpedloc]
	verboseDirs = False #for debug
	if makeDirs == True:
		for parentdir in parentdirs:
			checkMakeDir(parentdir, verbose = verboseDirs)
			for neutScenario in neutScenarios:
				scenariodir = parentdir + neutScenario
				checkMakeDir(scenariodir, verbose = verboseDirs)
				neutdir = scenariodir + "/neut"
				checkMakeDir(neutdir, verbose=verboseDirs)
				for selPop in selPops:
					selpopdir = parentdir +  neutScenario + "/sel" + str(selPop)
					checkMakeDir(selpopdir, verbose = verboseDirs)
					for selfreq in bin_medians_str:
						if parentdir != msloc:
							selfreqdir =  parentdir +  neutScenario + "/sel" + str(selPop) +  "/sel" + str(selfreq)
							checkMakeDir(selfreqdir, verbose = verboseDirs)
							if parentdir == tpedloc:
								terminalfolder = selfreqdir + "/" + scriptName
								checkMakeDir(terminalfolder, verbose = verboseDirs)

	#####################################
	## GENERATE SELECTION TRAJECTORIES ##
	#####################################
	if generateTrajs == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			startGens = {1:'U(0, 3000)', 2:'U(0, ' + str(paramDict[('split', 2)][0]-1) +')', 3:'U(0, ' + str(paramDict[('split', 3)][0]-1) +')', 4:'U(0, ' + str(paramDict[('split', 4)][0]-1) +')'} 
			for selPop in selPops:
				startGen, selCoeff = startGens[selPop], selCoeffs[selPop]
				startGenString = rewriteDistributionString(startGen)
				selCoeffString = rewriteDistributionString(selCoeff)
				sweepscenario = "sweep_" + str(selPop) 
				sweepParamFilename = paramsloc + neutScenarios[iscenario] +"_" + sweepscenario + ".par"
				label = '\"' + sweepscenario  + '\"'
				sweepline = "pop_event sweep_mult " + label + " " + str(selPop) + " " + str(startGen) + " " + str(selCoeff) + " " + str(sweepPos) + " " + str(presentDayFreq)
				sweepParamFilename = writeSweepToParamFile(neutParamFilename, sweepline, sweepscenario)
				argstring =  " -p " + sweepParamFilename + " --drop-singletons " + str(singRate) + " --genmapRandomRegions --output-gen-map --traj-only" #-j? 
				arguments = [argstring for w in range(nTrajPerScenario)]
				trajectoryname = trajloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/" + str(scriptName) + "_task${SGE_TASK_ID}_trajectory.tsv"
				commandstring = "env COSI_NEWSIM=1 COSI_SAVE_TRAJ=" + trajectoryname + " "+ cosi_build 
				writeJobsToTaskArray_fromArglist(commandstring, arguments, "generate_traj_" +  neutScenarios[iscenario] + "_" +sweepscenario)

	#####################
	#### CHECK OUTPUT ###
	#####################
	if checkSelFreqs == True:
		tallies = {}  #tallies[(selPop, freq)] = [list of trajectoryfiles]
		for iscenario in range(len(neutScenarios)):
			for selPop in selPops:
				fulltrajloc = trajloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/"
				alltrajectoryfiles = [fulltrajloc + item for item in os.listdir(fulltrajloc)]
				for selfreq in bin_medians_str:
					key = (str(selPop), str(selfreq))
					tallies[key] = []
				allfreqs = []
				for trajectoryfilename in alltrajectoryfiles:
					if os.path.isfile(trajectoryfilename):
						finalselfreq = getFinalSelFreq(trajectoryfilename, selpop = selPop)
						allfreqs.append(selfreq)
						for j in range(numbins+1):
							if finalselfreq >= bin_starts[j] and finalselfreq <= bin_ends[j]:
								binname = bin_medians_str[j]
								key = (str(selPop), binname)
								tallies[key].append(trajectoryfilename)
								break
				print "\tpulled final frequencies for " + str(len(allfreqs)) + " replicates in " + str(selPop)
			binkeys = tallies.keys()
			binkeys.sort()
			for key in binkeys:
				writefilename = trajloc + neutScenarios[iscenario] + "_" +str(key[0]) +"_" + str(key[1]) + "_trajs.txt"
				writefile = open(writefilename, 'w')
				for item in tallies[key]:
					writefile.write(item+'\n')
				writefile.close()
				print "wrote " + str(len(tallies[key])) + " trajectories to file: " + writefilename

	########################
	## RUN SEL SIMULATES ##
	########################
	if runCosiMS == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			startGens = {1:'U(0, 3000)', 2:'U(0, ' + str(paramDict[('split', 2)][0]-1) +')', 3:'U(0, ' + str(paramDict[('split', 3)][0]-1) +')', 4:'U(0, ' + str(paramDict[('split', 4)][0]-1) +')'} 
			for selPop in selPops:
				startGen, selCoeff = startGens[selPop], selCoeffs[selPop]
				startGenString = rewriteDistributionString(startGen)
				selCoeffString = rewriteDistributionString(selCoeff)
				sweepscenario = "sweep_" + str(selPop) #+ "_" + str(startGenString) + "_" + str(selCoeffString) 
				sweepParamFilename = paramsloc + neutScenarios[iscenario] +"_" + sweepscenario + ".par"
				assert os.path.isfile(sweepParamFilename)

				for selfreq in bin_medians_str:
					trajfilename = trajloc + neutScenarios[iscenario] + "_" +str(selPop) +"_" + str(selfreq) + "_trajs.txt"
					arguments = []
					key = (str(selPop), str(selfreq))
					openfile = open(trajfilename, 'r')
					for line in openfile:
						trajfile = line.strip('\n')
						fullargstring = "env COSI_LOAD_TRAJ="+ trajfile + " " + cosi_build + " -p " + sweepParamFilename +  " --drop-singletons " + str(singRate) + " --genmapRandomRegions --output-gen-map -m -P 15"
						arguments.append(fullargstring)
					commandstring = "env COSI_NEWSIM=1"

					print "running up to " + str(numPerScenario) + " reps per bin..."	
					if len(arguments) >	numPerScenario:
						arguments = arguments[:numPerScenario]

					if scaleSims == True:
						all_args = []
						for i in range(simScaleFactor):
							all_args.extend(arguments)
						arguments = all_args
					writeJobsToTaskArray_fromArglist(commandstring, arguments, "run_ms_sweeps_" + scriptName + neutScenarios[iscenario] + "_sel" + str(key[0]) + "_" + str(key[1]))

	######################
	## CLEAN UP MS SIMS ##
	######################
	if cleanUpMS == True:
		print "cleaning ms files generated from cosi run on UGER"
		allfiles = os.listdir(cwd)
		allfiles_filtered = [item for item in allfiles if "run_ms_sweeps" in item]
		print "loaded " + str(len(allfiles_filtered)) + " out/errfiles..."
		rm_args = []
		mv_args = []
		minSimLines = 50 #ensure outfiles is at least this long; depends on npop + nsamp but in general file should be a few hundred lines

		for iscenario in range(len(neutScenarios)):
			for selPop in selPops:
				for selfreq in bin_medians_str:
					key = (str(selPop), str(selfreq))
					selscenariostring = neutScenarios[iscenario] + "_sel" + str(key[0]) + "_" + str(key[1])
					thisScenarioFiles = [item for item in allfiles_filtered if selscenariostring in item and "taskargs" not in item]
					for item in thisScenarioFiles:
						if ".sh.o" in item:
							warn = True
						else:
							warn = False
						if checkFileLen(item, minNumLines = minSimLines, warn = warn) == item:
							movedir =  msloc + neutScenarios[iscenario]  + "/sel" + str(selPop) + "/sel" + str(selfreq)
							mv_argstring = item + " " + movedir
							mv_args.append(mv_argstring)
						else:
							rm_args.append(item)
		writeJobsToTaskArray_fromArglist("mv", mv_args, "mv_" + scriptName)# + neutScenarios[iscenario] + "_sel" + str(key[0]) + "_" + str(key[1]))
		writeJobsToTaskArray_fromArglist("rm", rm_args, "rm_" + scriptName)# + neutScenarios[iscenario] + "_sel" + str(key[0]) + "_" + str(key[1]))				

	####################
	## COUNT MS SIMS ##
	####################
	if countMS == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			for selPop in selPops:
				for selfreq in bin_medians_str: 
					key = (str(selPop), str(selfreq))
					msdir = msloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					print "loading from " + msdir
					ms_files = []
					
					#quick check for num files in dir
					entries = os.listdir(msdir)				
					print "\tloaded " + str(len(entries)) + " cosi ms files for " + str(key) + " under model " + str(neutScenarios[iscenario])
					if len(entries) < numPerScenario:
						print "\t\tERROR: <" + str(numPerScenario)
						#break

					if trimMS == True:
						if len(entries) > maxToKeep:
							keep = []
							while len(keep) < maxToKeep:
								random_outfile = choice(entries)
								if random_outfile not in keep:
									keep.append(random_outfile)
							discard = []
							for item in entries:
								if item not in keep:
									discard.append(item)
							print "deleting " + str(len(discard)) + " files randomly..."
							arguments = []
							for item in discard:
								arguments.append(msdir + item)
							writeJobsToTaskArray_fromArglist("rm", arguments, "trim_" + neutScenarios[iscenario] + "_sel" + str(selPop) + "_" + str(selfreq) + "_" + scriptName + "keep" + str(maxToKeep))

	##############################
	## CONVERT MS FILES TO TPED ##
	##############################
	if convertTped == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			for selPop in selPops:
				for selfreq in bin_medians_str: 
					key = (str(selPop), str(selfreq))
					msdir = msloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"
					print "loading from " + msdir
					ms_files = []
					
					entries = os.listdir(msdir)
					for item in entries:
						task_outfilename = msdir + item
						ms_files.append(task_outfilename)
					
					print "\tloaded " + str(len(ms_files)) + " cosi ms files for " + str(key) + " under model " + str(neutScenarios[iscenario])
					
					# SUBSETTING MS TO TPED
					print "\trandomly selecting " + str(numPerScenario) + " replicates for conversion..."
					convert_sims = []

					numPerScenario_mod = numPerScenario
					if len(ms_files) < numPerScenario_mod:
						numPerScenario_mod = len(ms_files)
						print "\t\tFewer than " + str(numPerScenario) + " ms simulations available. Limiting selection to " + str(numPerScenario_mod)

					while len(convert_sims) < numPerScenario_mod:
						random_outfile = choice(ms_files)
						if random_outfile not in convert_sims:
							convert_sims.append(random_outfile)
					#convert_sims = ms_files[:]

					tpeddir = tpedloc + neutScenarios[iscenario] + "/sel" + str(selPop) + "/sel" + str(selfreq) + "/"+ scriptName #+ "n" + str(numPerScenario)
					print "writing tpeds to: " + tpeddir

					argstrings = []
					for iTped in range(len(convert_sims)):
						sourcefilename = convert_sims[iTped]
						writefileprefix = tpeddir + "/rep" + str(iTped)
						argstring = sourcefilename + " " + writefileprefix + " " + str(regionlen) + " True"
						argstrings.append(argstring)
					commandstring = "python ms_to_tped.py"
					writeJobsToTaskArray_fromArglist(commandstring, argstrings, "convert_" + neutScenarios[iscenario] + "_sel" + str(selPop) + "_" + str(selfreq) + "_" + scriptName)

	########################
	## RUN NEUT SIMULATES ##
	########################
	if runNeut == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			fullargstring = "-p " + neutParamFilename +  " --drop-singletons " + str(singRate) + " --genmapRandomRegions --output-gen-map -m -P 15"
			arguments = [fullargstring for x in range(numNeut)]
			commandstring = cosi_build #"env COSI_NEWSIM=1 " + 
			writeJobsToTaskArray_fromArglist(commandstring, arguments, "run_ms_neutsims_" + scriptName + neutScenarios[iscenario])

	######################
	## CLEAN UP NEUT SIMS ##
	######################
	if cleanNeut == True:
		print "cleaning ms files generated from cosi run on UGER"
		allfiles = os.listdir(cwd)
		allfiles_filtered = [item for item in allfiles if "run_ms_neutsims_" in item]
		print "loaded " + str(len(allfiles_filtered)) + " out/errfiles..."
		rm_args = []
		mv_args = []
		minSimLines = 50 #ensure outfiles is at least this long; depends on npop + nsamp but in general file should be a few hundred lines
		for iscenario in range(len(neutScenarios)):
			selscenariostring = neutScenarios[iscenario]
			thisScenarioFiles = [item for item in allfiles_filtered if selscenariostring in item and "taskargs" not in item]
			for item in thisScenarioFiles:
				if ".sh.o" in item:
					warn = True
				else:
					warn = False
				if checkFileLen(item, minNumLines = minSimLines, warn = warn) == item:
					movedir = msloc + neutScenarios[iscenario]  + "/neut"
					mv_argstring = item + " " + movedir
					mv_args.append(mv_argstring)
				else:
					rm_args.append(item)
		writeJobsToTaskArray_fromArglist("mv", mv_args, "mv_" + scriptName + neutScenarios[iscenario] + "_neut")
		writeJobsToTaskArray_fromArglist("rm", rm_args, "rm_" + scriptName + neutScenarios[iscenario] + "_neut")				


	##############################
	## CONVERT NEUT SIMS TO TPED ##
	##############################
	if convertNeutTped == True:
		for iscenario in range(len(neutScenarios)):
			neutParamFilename = paramsloc + neutScenarios[iscenario] + ".par"
			print "**********\nLoading neutral parameters for scenario from " + neutParamFilename
			paramDict = getDictFromParamFile(neutParamFilename)
			regionlen = paramDict['chromlength']
			msdir = msloc + neutScenarios[iscenario] + "/neut/"
			print "loading from " + msdir
			ms_files = []
			entries = os.listdir(msdir)
			for item in entries:
				task_outfilename = msdir + item
				ms_files.append(task_outfilename)
			print "\tloaded " + str(len(ms_files)) + " cosi ms files for model " + str(neutScenarios[iscenario])

			tpeddir = tpedloc + neutScenarios[iscenario] + "/neut"
			print "writing tpeds to: " + tpeddir

			argstrings = []
			for iTped in range(len(ms_files)):
				sourcefilename = ms_files[iTped]
				writefileprefix = tpeddir + "/rep" + str(iTped)
				argstring = sourcefilename + " " + writefileprefix + " " + str(regionlen) + " False"
				argstrings.append(argstring)
			commandstring = "python ms_to_tped.py"
			writeJobsToTaskArray_fromArglist(commandstring, argstrings, "convert_" + neutScenarios[iscenario] + "_neut" + scriptName)

main()