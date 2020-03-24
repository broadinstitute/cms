##	run new, heterogenous demography model set; generate component scores
##	02.21.2019

ncmds_script, dispatchNow, subopts = 100, False, "-l h_vmem=5g -l os=RedHat6 "#-pe smp 8 -binding linear:8 -R y "
basecmd = "python /idi/sabeti-scratch/jvitti/remodel_freqs_4pop.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/remodel/", "gravel_coverage"

import time
import os
from random import shuffle
import subprocess

models, regimes, pops, reps = ['gradient15'], ['hard', 'soft'], [1,2,3,4], list(range(1,5001))#range(1001, 5000)#
basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run/"

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
	scriptfile.write('#$ -o ' + fileprefix + ".o\n")
	scriptfile.write('#$ -e ' + fileprefix + ".e\n")	

	scriptfile.write('export PATH=/idi/sabeti-scratch/jvitti/miniconda3/bin:$PATH\n')
	if activate_venv:
		scriptfile.write('source activate /idi/sabeti-scratch/jvitti/miniconda3/envs/cms-env3\n')
	for iarg in range(len(arguments)):
		this_cmd = commandstring + " " + arguments[iarg]
		scriptfile.write(this_cmd + "\n")
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

	replicate_idstring = "rep" + str(irep)
	
	basedir = basewritedir + model + "_" + regime + "_sel" + str(selpop) + "/"

	altpops = [1, 2, 3, 4]
	altpops.remove(selpop)
	alt1, alt2, alt3 = altpops

	ihs_filename = basedir + "ihs/" + replicate_idstring + ".ihs.out"
	delihh_filename = basedir + "delihh/" + replicate_idstring
	nsl_filename = basedir + "nsl/" + replicate_idstring + ".nsl.out"
	xpehh_filename1 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt1) + ".xpehh.out"
	xpehh_filename2 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt2) + ".xpehh.out"
	#xpehh_filename3 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt3) + ".xpehh.out"
	freq_filename1 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt1)
	freq_filename2 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt2)
	freq_filename3 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt3)

	checkfilenames = [freq_filename1, freq_filename2, freq_filename3]

	for checkfilename in checkfilenames:
		if not os.path.isfile(checkfilename) or os.path.getsize(checkfilename) == 0:
			needed = True
			break
	return needed

def main():
	arguments = []
	for model in models:
		for regime in regimes:
			for pop in pops:
				this_rundir = basewritedir + model +"_" + regime + "_sel" + str(pop) + "/"
				for irep in reps:
					tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
					if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
						argstring = "\t".join([model, regime, str(pop), str(irep)])
						if checkNeeded(argstring):		
							arguments.append(argstring)
							full_cmd = basecmd + " " + argstring
							print(full_cmd)
							subprocess.check_output(full_cmd.split())


	##############
	## DISPATCH ##
	##############
	shuffle(arguments)
	print(('loaded a total of ' + str(len(arguments)) + " commands."))
	iscript = 0
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()