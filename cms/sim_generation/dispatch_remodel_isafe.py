##	/idi/sabeti-scratch/jvitti/remodel/run/
##	02.20.2019	#last used 11.17.2019

ncmds_script, subopts, scriptprefix, dispatchNow = 40, "-l h_vmem=9g ", "isafe_remodel", False
scriptdir, basecmd = "/idi/sabeti-scratch/jvitti/scripts/", "python /home/unix/vitti/run_remodel_isafe.py"

import os
import sys
import subprocess
from random import shuffle
import time

models = [ 'gradient15', 'defdef15', 'def15']#'nd', 'ndcs', 'gravel']#'gravel'] #--drop-singletons " + str(singrate) <--- absent from gravel as run! presnt in gradient15!
regimes = ['lct', 'edar', 'slc24a5']#['neut']
pops = [1,2,3,4]
reps = list(range(0,1001))
#models, regimes, pops, reps = ['def15',], ['neut'], [1,2,3,4], range(0,1001)#range(1001, 5000)#
#models, regimes, pops, reps = ['defdef15', 'gradient15', 'def15', 'gravel', 'nd', 'ndcs'], ['neut'], [1,2,3,4], range(0,1001)#range(1001, 5000)#
#models, regimes, pops, reps = ['gravel'], ['neut'], [1,2,3,], range(0,1001)#range(1001, 5000)#
basedir = "/idi/sabeti-scratch/jvitti/remodel/run4/"

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-l h_rt=48:00:00'], activate_cms = True, dispatch = False, subopts = ''): #re-add arg for dependencies
	""" writes a script for all of these tasks to be performed sequentially """

	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("\n")
	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	if activate_cms:
		scriptfile.write('source activate /idi/sabeti-scratch/jvitti/miniconda3/envs/cms-env3\n')
	else:
		scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/iSafe\n')
	scriptfile.write('source /broad/software/scripts/useuse\nreuse BEDTools\nreuse BCFtools\nreuse use .samtools-1.8\nreuse .tabix-0.2.6\nreuse .bcftools-1.8\n')
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

def main():
	arguments = []
	for model in models:
		for regime in regimes:
			if regime == "edar":
				pops = [3]
			else:
				pops = [2]
			for pop in pops:
				if not (pop == 4 and model == "gravel"):
					this_rundir = basedir + model +"_" + regime + "/"#"_sel" + str(pop) + "/"
					for extension in ['isafe', 'isafe_data']:
						create_dircmd = "mkdir -p " + this_rundir + extension + "/"
						subprocess.check_output( create_dircmd.split() )
					for irep in reps:
						tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
						if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
							checkfilename = this_rundir + "isafe/OG_rep" + str(irep) + "_0_pop" + str(pop) + ".iSAFE.out"
							if not os.path.isfile(checkfilename):
								argstring = "\t".join([model, regime, str(pop), str(irep)])
								arguments.append(argstring)
					

	print(('loaded a total of ' + str(len(arguments)) + " commands."))
	iscript = 0
	shuffle(arguments)
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, activate_cms = False, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()
