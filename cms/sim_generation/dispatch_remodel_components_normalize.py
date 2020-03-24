##	run new, heterogenous demography model set; generate component scores
##	02.21.2019 	06.28: run2

ncmds_script, dispatchNow, subopts = 300, True, " "
basecmd = "python /home/unix/vitti/remodel_components_normalize.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/remodel/", "rem_norm"

import time
import os
from random import shuffle
import subprocess

models = [ 'gradient15', 'defdef15', 'def15']
regimes = ['lct', 'slc24a5', 'edar']
reps = list(range(0,1001))
basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run4/"  
#models, regimes, pops, reps = ['nd'], ['neut'], [1,2,3, 4], range(0,1001) #'gravel', 
#models, regimes, pops, reps = ['gradient15', 'defdef15', 'def15', 'gravel'], ['neut'], [1,2,3, 4], range(1,1001) #'gravel', 
#models, regimes, pops, reps = ['gradient15', 'defdef15', 'def15', 'gravel'], ['neut'], [1,2,3, 4], range(1,1001) #'gravel', 
#models, regimes, pops, reps = ['defdef15', 'nd', 'ndcs', 'gradient15', 'def15'], ['neut'], [1,2,3, 4], range(0,1001) #'gravel', 
#models, regimes, pops, reps = ['gravel',], ['neut'], [1,2,3, 4], range(0,1001) #'gravel', 
#basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run3/"

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-l h_rt=48:00:00'], activate_venv = True, dispatch = False, subopts = ''): #re-add arg for dependencies
	""" writes a script for all of these tasks to be performed sequentially """
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("\n")
	scriptfile.write('export PATH=//home/unix/vitti/miniconda3/bin:$PATH\n')
	if activate_venv:
		scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/py2\n')
		#scriptfile.write('source /broad/software/scripts/useuse\nreuse BEDTools\nreuse VCFtools\n')
	for iarg in range(len(arguments)):
		this_cmd = commandstring + " " + arguments[iarg]
		scriptfile.write(this_cmd + "\n")
	#writestring = commandstring + " $SAMPLE\n"
	#scriptfile.write(writestring)
	scriptfile.close()
	subcommand = "qsub " + subopts + scriptname
	print(subcommand)
	if dispatch:
		subprocess.check_output( subcommand.split() )
	return subcommand
def checkNeeded(argstring):
	needed = False
	entries = argstring.split('\t')
	model, regime, pop, irep = entries
	pop = int(pop)
	selpop = pop
	irep = int(irep)

	replicate_idstring = "rep" + str(irep) + "_pop" + str(pop)
	
	basedir = basewritedir + model + "_" + regime  + "/"#+ "_sel" + str(selpop) + "/"


	#altpops = [1, 2, 3]
	altpops = [1, 2, 3, 4]
	altpops.remove(selpop)
	alt1, alt2, alt3 = altpops
	#alt1, alt2 = altpops


	ihh12_filename = basedir + "ihh12/" + replicate_idstring + ".ihh12.out.norm"
	if not os.path.isfile(ihh12_filename):
		ihh12_filename += ".gz"



	ihs_filename = basedir + "ihs/" + replicate_idstring + ".ihs.out.norm"
	if not os.path.isfile(ihs_filename):
		ihs_filename += ".gz"

	delihh_filename =  basedir + "delihh/" + replicate_idstring + ".norm"
	if not os.path.isfile(delihh_filename):
		delihh_filename += ".gz"

	nsl_filename = basedir + "nsl/" + replicate_idstring + ".nsl.out.norm"
	if not os.path.isfile(nsl_filename):
		nsl_filename += ".gz"

	xpehh_filename1 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt1) + ".xpehh.out.norm"
	if not os.path.isfile(xpehh_filename1):
		xpehh_filename1 += ".gz"

	
	xpehh_filename2 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt2) + ".xpehh.out.norm"
	if not os.path.isfile(xpehh_filename2):
		xpehh_filename2 += ".gz"

	xpehh_filename3 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt3) + ".xpehh.out.norm"
	if not os.path.isfile(xpehh_filename2):
		xpehh_filename3 += ".gz"


	#checkfilenames = 
	#[ihs_filename, delihh_filename, nsl_filename, xpehh_filename1,
	checkfilenames = [delihh_filename, nsl_filename, ihs_filename, xpehh_filename1, xpehh_filename2,xpehh_filename3,]
						#]#, freq_filename1, freq_filename2]
						# freq_filename1, 
						#freq_filename2, freq_filename3]

	for checkfilename in checkfilenames:
		if not (os.path.isfile(checkfilename)):# or os.path.getsize(checkfilename) == 0):
			needed = True
			#print(checkfilename)
			break
		elif os.path.getsize(checkfilename) == 0:
			needed = True
			#print(checkfilename)
			break
	return needed

def main():
	arguments = []
	for model in models:
		for regime in regimes:
			if regime == "edar":
				pops = [3]
			else:
				pops = [2]
			#for pop in pops:
			for pop in pops:
				if not (model == "gravel" and pop == 4):
					this_rundir = basewritedir + model +"_" + regime +"/"# "_sel" + str(pop) + "/"
					for irep in reps:
						tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
						if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
							argstring = "\t".join([model, regime, str(pop), str(irep)])
							if checkNeeded(argstring):		
							#if True:
								arguments.append(argstring)
							#else:
							#	print('passing')


	##############
	## DISPATCH ##
	##############
	#shuffle(arguments)
	print(('loaded a total of ' + str(len(arguments)) + " commands."))
	iscript = 0
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()
