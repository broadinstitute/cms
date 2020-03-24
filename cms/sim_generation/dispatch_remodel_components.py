##	run new, heterogenous demography model set; generate component scores
##	02.21.2019

ncmds_script, dispatchNow, subopts = 15, True, " -l h_vmem=1g  -pe smp 8 -binding linear:8 -R y " ##-l os=RedHat6 "# "-l h_vmem=15g -l os=RedHat6 
basecmd = "python /home/unix/vitti/remodel_components.py"
#remodel_freqs.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/remodel/", "jv_remodel4"

import time
import os
from random import shuffle
import subprocess

models = [ 'gradient15', 'defdef15', 'def15']#'nd', 'ndcs', 'gravel']#'gravel'] #--drop-singletons " + str(singrate) <--- absent from gravel as run! presnt in gradient15!
regimes = ['lct', 'slc24a5']
#'edar']#'lct', 'edar', 'slc24a5']#['neut']
pops = [1,2,3,4]
reps = list(range(0,1001))
#models, regimes, pops, reps = ['def15', 'defdef15', 'gradient15', 'nd', 'ndcs', 'gravel'], ['neut'], [1,2,3, 4], range(1,1001)
#models, regimes, pops, reps = ['gradient15',], ['neut'], [1,2,3, 4], range(0,1001)
basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run4/" #'nd', 'ndcs', 

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-l h_rt=48:00:00'], activate_venv = False, dispatch = False, subopts = ''): #re-add arg for dependencies
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
	
	basedir = basewritedir + model + "_" + regime + "/"#+ "_sel" + str(selpop) + "/"

	altpops = [1, 2, 3, 4]
	altpops.remove(selpop)
	alt1, alt2, alt3 = altpops

	ihs_filename = basedir + "ihs/" + replicate_idstring + ".ihs.out"
	delihh_filename = basedir + "delihh/" + replicate_idstring
	nsl_filename = basedir + "nsl/" + replicate_idstring + ".nsl.out"
	xpehh_filename1 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt1) + ".xpehh.out"
	xpehh_filename2 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt2) + ".xpehh.out"
	xpehh_filename3 = basedir + "xpehh/" + replicate_idstring + "_vs" + str(alt3) + ".xpehh.out"
	freq_filename1 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt1)

	if not os.path.isfile(freq_filename1):
		freq_filename1 += ".gz"
	freq_filename2 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt2)
	if not os.path.isfile(freq_filename2):
		freq_filename2 += ".gz"

	freq_filename3 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt3)
	if not os.path.isfile(freq_filename3):
		freq_filename3 += ".gz"

	#checkfilenames = [delihh_filename]
	#[ihs_filename, delihh_filename, nsl_filename, xpehh_filename1,
	#					xpehh_filename2, freq_filename1, freq_filename2]
	#checkfilenames = [freq_filename1, freq_filename2, freq_filename3]
	checkfilenames = [xpehh_filename3, xpehh_filename1, xpehh_filename2, delihh_filename, ihs_filename, nsl_filename]
						#freq_filename2, freq_filename3]

	for checkfilename in checkfilenames:
		if not os.path.isfile(checkfilename) or os.path.getsize(checkfilename) == 0:
			needed = True
			print(checkfilename)
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
			for pop in pops:
				this_rundir = basewritedir + model +"_" + regime + "/"
				
				for necessary_dir in ['ihh12', 'ihs', 'delihh', 'freqs', 'xpehh', 'nsl']:
					mkdir_cmd = "mkdir -p " + this_rundir + necessary_dir + "/"
					subprocess.check_output( mkdir_cmd.split() )

				#this_rundir = basewritedir + model +"_" + regime + "_sel" + str(pop) + "/"
				for irep in reps:
					tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
					if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
						#print(tpedfilename)
						argstring = "\t".join([model, regime, str(pop), str(irep)])
						if checkNeeded(argstring):		
						#if True:
							arguments.append(argstring)
							#full_cmd = basecmd + " " + argstring
							#print(full_cmd)
							#subprocess.check_output(full_cmd.split())


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

