##	run new, heterogenous demography model set; generate component scores
##	02.21.2019

ncmds_script, dispatchNow, subopts = 100, True, "-l h_vmem=5g "
basecmd = "python /home/unix/vitti/remodel_freqs_4pop.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/remodel/", "jv_freqs"

import time
import os
from random import shuffle
import subprocess
models = [ 'gradient15', 'defdef15', 'def15']#'nd', 'ndcs', 'gravel']#'gravel'] #--drop-singletons " + str(singrate) <--- absent from gravel as run! presnt in gradient15!
regimes = ['lct', 'edar', 'slc24a5']#['neut']
pops = [1,2,3,4]
reps = list(range(0,1001))
#models, regimes, pops, reps = ['gravel'], ['neut'], [1,2,3], range(0,1001)
#models, regimes, pops, reps = ['defdef15', 'def15', 'nd', 'ndcs', 'gradient15'], ['neut'], [1,2,3, 4], range(0,1001)
basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run4/"

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
	scriptfile.write('#$ -o ' + fileprefix + ".o\n")
	scriptfile.write('#$ -e ' + fileprefix + ".e\n")	

	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	if activate_venv:
		scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/py2\n')
	scriptfile.write('source /broad/software/scripts/useuse\n')
	scriptfile.write('reuse -q .google-cloud-sdk\n')
	for iarg in range(len(arguments)):
		this_cmd = commandstring + " " + arguments[iarg]
		scriptfile.write(this_cmd + "\n")
	#scriptfile.write('gsutil -m cp -r /idi/sabeti-scratch/jvitti/remodel/ gs://vhoorl/jvitti/remodel/')
	#scriptfile.write('gsutil -m rsync /idi/sabeti-scratch/jvitti/remodel gs://vhoorl/jvitti/remodel/\n')
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
	
	basedir = basewritedir + model + "_" + regime + '/'#"_sel" + str(selpop) + "/"
	altpops = [1, 2,3, 4]
	altpops.remove(selpop)
	if not "gravel" in argstring:
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


		checkfilenames = [freq_filename1, freq_filename2, freq_filename3]
		#[ihs_filename, delihh_filename, nsl_filename, xpehh_filename1,
		#					xpehh_filename2, xpehh_filename3, ] #freq_filename1, freq_filename2]
							#
							#freq_filename1, 
							#freq_filename2, freq_filename3]
	else:
		altpops.remove(4)
		alt1, alt2 = altpops		
		freq_filename1 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt1)
		if not os.path.isfile(freq_filename1):
			freq_filename1 += ".gz"
		freq_filename2 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt2)
		if not os.path.isfile(freq_filename2):
			freq_filename2 += ".gz"
		checkfilenames = [freq_filename1, freq_filename2]
			

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
#				if True:
				if not (model == "gravel" and pop == 4):
					this_rundir = basewritedir + model +"_" + regime + "/" #+ "_sel" + str(pop) + "/"
					mkdir_cmd = "mkdir -p " + this_rundir + "freqs/"
					for irep in reps:
						tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
						if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
							argstring = "\t".join([model, regime, str(pop), str(irep)])
							#print(tpedfilename)
							if checkNeeded(argstring):		
								arguments.append(argstring)

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
