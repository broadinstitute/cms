##	run new, heterogenous demography model set; generate component scores
##	02.21.2019	## RUN HAPLOTYPE SCORES FOR POP AS OUTGROUP FOR NORMALIZE

ncmds_script, dispatchNow, subopts = 1, True, "-l h_vmem=1g -pe smp 8 -binding linear:8 -R y "
basecmd = "/idi/sabeti-scratch/jvitti/selscan/bin/linux/selscan"
#basecmd = "python /home/unix/vitti/remodel_components_neut.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/remodel/", "remodel_neut_ex"

import time
import os
from random import shuffle
import subprocess

#'defdef15', 'gradient15'
models, regimes, pops, reps = ['nd'], ['neut'], [1,2,3,4], range(1,1000)#range(1001, 5000)#
#models, regimes, pops, reps = ['nd', 'ndcs', 'def15', 'defdef15', 'gradient15', 'gravel' ], ['neut'], [1,2,3,4], range(0,1001)#range(1001, 5000)#
#models, regimes, pops, reps = [ 'gradient15'], ['neut'], [1,2,3,4], range(0,1001)#range(1001, 5000)#

basewritedir = "/idi/sabeti-scratch/jvitti/remodel/run3/"

def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-l h_rt=48:00:00'], activate_venv = True, dispatch = False, subopts = ''): #re-add arg for dependencies
	""" writes a script for all of these tasks to be performed sequentially """
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("\n")
	scriptfile.write('export PATH=//home/unix/vitti/miniconda3/bin:$PATH\n')
	#if activate_venv:
	#	scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/py2\n')
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
def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def checkNeeded(argstring):
	needed = False
	entries = argstring.split('\t')
	model, regime, pop, irep, use_pop = entries
	pop = int(pop)
	selpop = pop
	irep = int(irep)
	use_pop = int(use_pop)

	replicate_idstring = "rep" + str(irep)
	
	basedir = basewritedir + model + "_" + regime + "_sel" + str(use_pop) + "/"

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
	freq_filename2 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt2)
	#freq_filename3 = basedir + "freqs/"  + replicate_idstring + "_vs" + str(alt3)

	#checkfilenames = [ihs_filename, delihh_filename, nsl_filename, xpehh_filename1,
	#					xpehh_filename2, freq_filename1, freq_filename2,
	#					xpehh_filename3]#, freq_filename1, 
						#freq_filename2, freq_filename3]
	#checkfilenames = [xpehh_filename3, xpehh_filename2, xpehh_filename1]
	checkfilenames = [delihh_filename]
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
				if not (model == "gravel" and pop == 4):
						
					use_pop = pop
					this_rundir = basewritedir + model +"_" + regime + "/"
					for irep in reps:
						tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
						if not os.path.isfile(tpedfilename):
							tpedfilename += ".gz"
						#else:
						#	print('missing ' + tpedfilename)
						
						#argstring = model + " " + regime + " " + str(pop)  + " " +  str(irep)  + " " + str(pop)
						#arguments.append(argstring)
						
						
						for score in ['ihh12', 'ihs', 'nsl']:
							outfilename = this_rundir + score + "/rep" + str(irep) + "_pop" + str(pop) 
							if os.path.isfile(tpedfilename) and (not os.path.isfile(outfilename + "." + score + ".out") or os.path.getsize(outfilename + "." + score + ".out") == 0):
								argstring = " --" + score + " --tped " + tpedfilename + " --out " + outfilename + " --threads 8"
								arguments.append(argstring)
						

						for score in ['xpehh']:
							altpops = pops[:]
							altpops.remove(pop)
							for altpop in altpops:
								if not (model == "gravel" and altpop == 4):
									outfilename = this_rundir + score + "/rep" + str(irep) + "_pop" + str(pop) + "_vs" + str(altpop)
									alt_tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(altpop) + ".tped"
									if not os.path.isfile(alt_tpedfilename):
										alt_tpedfilename += ".gz"
									if os.path.isfile(tpedfilename) and (not os.path.isfile(outfilename + "." + score + ".out") or os.path.getsize(outfilename + "." + score + ".out") == 0):
										argstring = " --" + score + " --tped " + tpedfilename + " --out " + outfilename + " --threads 8"
										argstring += " --tped-ref " + alt_tpedfilename
										arguments.append(argstring)

	##############
	## DISPATCH ##
	##############
	shuffle(arguments)
	print('loaded a total of ' + str(len(arguments)) + " commands.")
	iscript = 0
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()
