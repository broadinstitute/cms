##	/n/scratchlfs/sabeti_lab/vitti/remodel_samplesize/
##	06.26.2019

ncmds_script, dispatchNow = 5, False
scriptdir, scriptprefix = "/n/home08/jvitti/remodel_rc/scripts/", "missing_sims"

import os
import sys
import subprocess
from random import shuffle
import time

basedir = "/n/scratchlfs/sabeti_lab/vitti/remodel_samplesize/"
reps, ns = range(1, 1001), [25, 50, 100, 150, 200, 250, 500]#25, 50, 100, 150] #ADD 200, 250, 500 #,,,____ SOFT AND HARD
basecmd = "python /n/home08/jvitti/remodel_rc/rerun_cosi.py"

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def write_quick_slurm_script(commandstrings, scriptname):
	openfile = open(scriptname, 'w')
	openfile.write('#!/bin/bash\n')
	openfile.write('#SBATCH -p general\n')
	openfile.write('#SBATCH -c 1\n')
	openfile.write('#SBATCH -N 1\n')
	openfile.write('#SBATCH -t 6-12:00\n')
	openfile.write('#SBATCH --mem 50000\n') #BASTA? SATIS? 
	openfile.write('#SBATCH -o  '+ scriptname+ "_o_%j\n")
	openfile.write('#SBATCH -e ' + scriptname + '_e_%j\n')
	openfile.write('source activate cms-env3\n')
	for commandstring in commandstrings:
		openfile.write(commandstring + "\n")
	openfile.close()
	dispatchcmd = "sbatch " + scriptname
	return dispatchcmd

def main():
	commands = []
	for n in ns:
		for batch in ['', '_soft']:
			this_rundir = "/n/scratchlfs/sabeti_lab/vitti/remodel_samplesize" + batch + "/n" + str(n) + "/"#basedir + model +"_" + regime + "_sel" + str(pop) + "/"
			paramfilename = "/n/home08/jvitti/remodel_samplesize/params/default_n" + str(n) + batch + ".par"
			if not os.path.isfile(paramfilename) or os.path.getsize(paramfilename) == 0:
				print('missing ' + paramfilename)
			for irep in reps:
				trajfilename = this_rundir + "sampled_vars/rep" + str(irep) + ".traj"
				save_sample_filename = this_rundir + "sampled_vars/rep" + str(irep) + ".sampled"
				checkfilename = this_rundir + "tpeds/rep" + str(irep)
				argstring = "env COSI_NEWSIM=1 COSI_MAXATTEMPTS=1000000 COSI_SAVE_TRAJ=" + trajfilename + " COSI_SAVE_SAMPLED=" + save_sample_filename + " coalescent -p " + paramfilename + " --genmapRandomRegions --drop-singletons .25 --output-gen-map "
				argstring += " --tped  " + checkfilename  #MUST BE LAST ARG
				if not os.path.isfile(checkfilename + "_0_1.tped") and not os.path.isfile(checkfilename + "_0_1.tped.gz"):
					cmdstring = basecmd + " " +argstring
					commands.append(cmdstring)
					print(cmdstring)


	print('loaded a total of ' + str(len(commands)) + " commands.")
	iscript = 0
	shuffle(commands)
	for commandchunk in chunks(commands,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		slurmcmd = write_quick_slurm_script(commandchunk, scriptname)
		print(slurmcmd)
		if dispatchNow:
			subprocess.check_output( slurmcmd.split() )
		time.sleep(1)


main()
