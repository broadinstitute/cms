##	/idi/sabeti-scratch/jvitti/remodel/run/
##	02.20.2019 	#last used 06.13.2019 #last used neut 11.16.19

ncmds_script, dispatchNow = 5, True
subopts, scriptprefix = "-l h_vmem=9g ", "sims_pc"
scriptdir = "/idi/sabeti-scratch/jvitti/scripts/"

import os
import sys
import subprocess
from random import shuffle
import time
#'def15']#,
models = [ 'gradient15', 'defdef15', 'def15']#'nd', 'ndcs', 'gravel']#'gravel'] #--drop-singletons " + str(singrate) <--- absent from gravel as run! presnt in gradient15!
regimes = ['lct', 'edar', 'slc24a5']#['neut']
basedir = "/idi/sabeti-scratch/jvitti/remodel/run4/" #NEUT #run/"
reps = range(0, 1001)
basecmd = "python /home/unix/vitti/rerun_cosi.py"

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

def main():
	arguments = []
	for model in models:
		for regime in regimes:
				this_rundir = basedir + model +"_" + regime  + "/"
				for item in ['sampled_vars', 'tpeds']:
					mkdir_cmd = "mkdir -p " + this_rundir + '/' + item + "/"
					subprocess.check_output( mkdir_cmd.split() )
				paramfilename = basedir + "params/" + model + "_" + regime + ".par" 
				if not os.path.isfile(paramfilename) or os.path.getsize(paramfilename) == 0:
					print('missing ' + paramfilename)
				for irep in reps:
					trajfilename = this_rundir + "sampled_vars/rep" + str(irep) + ".traj"
					save_sample_filename = this_rundir + "sampled_vars/rep" + str(irep) + ".sampled"
					checkfilename = this_rundir + "tpeds/rep" + str(irep)
					argstring = "env COSI_NEWSIM=1 COSI_MAXATTEMPTS=1000000 COSI_SAVE_TRAJ=" + trajfilename + " COSI_SAVE_SAMPLED=" + save_sample_filename + " coalescent -p " + paramfilename + " --genmapRandomRegions --drop-singletons .25 --output-gen-map "
					argstring += " --tped  " + checkfilename  #MUST BE LAST ARG
					if not os.path.isfile(checkfilename + "_0_1.tped"):# and not os.path.getsize(checkfilename + "_0_1.tped" == 0):#os.path.isfile(checkfilename + "_0_1.tped.gz"):
						arguments.append(argstring)

	
	print('loaded a total of ' + str(len(arguments)) + " commands.")
	iscript = 0
	shuffle(arguments)
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, activate_venv = True, dispatch=dispatchNow, subopts = subopts)

		time.sleep(1)


main()
