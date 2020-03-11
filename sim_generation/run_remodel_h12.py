##	01.16.2019

ncmds_script, dispatchNow, subopts = 10, True, " -l os=RedHat6 -l h_vmem=3g -M vitti@broadinstitute.org -m eas "#"-l h_vmem=6g "
h12_script = "/idi/sabeti-scratch/jvitti/h12/SelectionHapStats/scripts/H12_H2H1.py"
scriptdir, scriptprefix = "/idi/sabeti-scratch/jvitti/scripts/", "PEL_h12"

import time
import os
from random import shuffle
#import numpy as np
import subprocess

freqBins = [0.05, 0.15, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
#np.arange(.05,1,.05)  #rounding error
nPerBin_native, nPerBin_foreign, selpops = 25, 100, [1, 2, 3, 4]
basewritedir = "/idi/sabeti-scratch/jvitti/PEL_model/run2/"#standing_sims_batch4/"
numIndivs = 170

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-l h_rt=48:00:00',], activate_venv = True, dispatch = False, subopts = ''): #re-add arg for dependencies
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
	scriptfile.write('source /broad/software/scripts/useuse\nuse .python-2.7.9-sqlite3-rtrees\nuse .numpy-1.8.0-python-2.7.8-sqlite3-rtrees\n')
	
	for iarg in range(len(arguments)):
		this_cmd = commandstring + " " + arguments[iarg]
		scriptfile.write(this_cmd + "\n")
	scriptfile.close()
	subcommand = "qsub " + subopts + scriptname
	print(subcommand)
	if dispatch:
		subprocess.check_output( subcommand.split() )
	return subcommand
def get_num_from_formatted(infilename):
	openfile = open(infilename, 'r')
	firstline = openfile.readline()
	openfile.close()
	numIndivs = firstline.count(',')
	return numIndivs

def main():
	arguments = []

	##############################
	## ITERATE OVER SELDAF BINS ## we want the distribution to be ~uniform.
	##############################
	for freqBin in freqBins:
		#####################
		## DOMESTIC SWEEPS ##  we want this to be ~60% of sims.
		#####################
		for selpop in [4]:#selpops:
			ancestralpop = selpop
			basedir = basewritedir + "sel" + str(selpop) + "/"
			for irep in range(1, nPerBin_native+1):
				replicate_idstring = "pop" + str(selpop) + "_native_seldaf" + str(freqBin) + "_rep" + str(irep) 
				inputfilename = basedir + "h12_data/" + replicate_idstring + "_0_" + str(selpop)
				if os.path.isfile(inputfilename) and os.path.getsize(inputfilename) > 0:
					h12_cmd = "python2.7 " + h12_script 
					#numIndivs = 172#get_num_from_formatted(inputfilename)
					#print(numIndivs)
					outfilename = basedir + 'h12/' + replicate_idstring
					h12_argstring = inputfilename + " " + str(numIndivs) + " -j 1 > " + outfilename
					h12_fullcmd = h12_cmd + " " +h12_argstring
					if not os.path.isfile(outfilename) or os.path.getsize(outfilename) == 0:
						arguments.append(h12_argstring)
						
		
		####################
		## FOREIGN SWEEPS ##  we want this to be ~40% of sims.
		####################	6 pop/ancpop combos
		for selpop in [4]:#selpops:
			basedir = basewritedir + "sel" + str(selpop) + "/"
			ancestralpops = [item for item in selpops if item < selpop]
			tped_dir = basedir + "tpeds/"	
			for ancestralpop in ancestralpops:
				for irep in range(1, nPerBin_foreign+1):
					replicate_idstring = "pop" + str(selpop) + "_foreign" + str(ancestralpop) + "_seldaf" + str(freqBin) + "_rep" + str(irep) 
					inputfilename = basedir + "h12_data/" + replicate_idstring + "_0_" + str(selpop)
					if os.path.isfile(inputfilename) and os.path.getsize(inputfilename) > 0:
						h12_cmd = "python2.7 " + h12_script 
						#numIndivs = get_num_from_formatted(inputfilename)
						outfilename = basedir + 'h12/' + replicate_idstring
						h12_argstring = inputfilename + " " + str(numIndivs) + " -j 1 > " + outfilename
						h12_fullcmd = h12_cmd + " " +h12_argstring
						if not os.path.isfile(outfilename) or os.path.getsize(outfilename) == 0:
							arguments.append(h12_argstring)
							


	##############
	## DISPATCH ##
	##############
	shuffle(arguments)
	print('loaded a total of ' + str(len(arguments)) + " commands.")
	iscript = 0
	for argchunk in chunks(arguments,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(h12_cmd, argchunk, scriptname, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()
