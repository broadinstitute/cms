##	02.27.2019		## run2! 08.06.2019 #remodel_extra

dispatchNow, subopts, ncmds_script = True, "", 50
scriptdir, scriptprefix  = "/idi/sabeti-scratch/jvitti/scripts/", "collate_remodel2"
this_basedir  = "/idi/sabeti-scratch/jvitti/remodel_extra/"
basecmd = "/idi/sabeti-scratch/jvitti/cms/cms/combine/collate_scores_v2_maf"

import os
import time
import subprocess
import importlib.machinery
parse_func = importlib.machinery.SourceFileLoader('parse_func','/idi/sabeti-scratch/jvitti/cms/cms/power/parse_func.py').load_module()
from parse_func import write_pair_sourcefile
jv_uger_util = importlib.machinery.SourceFileLoader('jv_uger_util','/home/unix/vitti/jv_uger_util.py').load_module()
from jv_uger_util import chunks, check_create_dir, writeJobsToBatch_fromArglist

#models, regimes, pops, reps = ['ndcs', 'nd',],  ['hard', 'soft'], [1,2,3, 4], range(1,1001)#range(1001, 5000)#
#['gradient15', 'defdef15',],  ['hard', 'soft'], [1,2,3, 4], range(1,1001)#range(1001, 5000)#

#'gradient15', 'defdef15', 'from_RC/nd', 'from_RC/ndcs'

subdirs = ['negsel/']
"""
ns = [25, 50, 100, 150, 200, 250, 500]
for n in ns:
	subdirs.extend(['samplesize/n' + str(n) + "/", 'samplesize_soft/n' + str(n) + "/", ])
"""
def get_remodel_components(basedir, replicate_idstring, altpop, thispop):
	#if "negsel" not in basedir:

	ihs_filename = basedir + "ihs/" + replicate_idstring + "_" + str(thispop) +".ihs.out" + ".norm"
	if not os.path.isfile(ihs_filename):
		ihs_filename += ".gz"

	delihh_filename =  basedir + "delihh/" + replicate_idstring + "_1.txt.norm"
	if not os.path.isfile(delihh_filename):
		delihh_filename += ".gz"

	nsl_filename = basedir + "nsl/" + replicate_idstring + "_" + str(thispop) + ".nsl.out" + ".norm"
	if not os.path.isfile(nsl_filename):
		nsl_filename += ".gz"

	xpehh_filename = basedir + "xpehh/" + replicate_idstring + "_" + str(altpop) + "_" + str(thispop)+ ".xpehh.out" + ".norm"
	if not os.path.isfile(xpehh_filename):
		xpehh_filename += ".gz"

	freqs_filename = basedir + "freqs/" + replicate_idstring + "_" + str(altpop) + "_" + str(thispop) #rep" + str(repNum) + "_sel" + str(selpop) + "_" + str(altpop) + "_orig" + str(ancpop) + "_" + str(thispop)
	if not os.path.isfile(freqs_filename):
		freqs_filename += ".gz"

	h12_filename = basedir + "h12/" + replicate_idstring + "_" + str(thispop)#rep" + str(repNum) + "_sel" + str(selpop) + "_orig" + str(ancpop) + "_" + str(thispop)
	if not os.path.isfile(h12_filename):
		h12_filename += ".gz"

	isafe_filename = basedir +  "isafe/OG_" + replicate_idstring[:-5] + "_0_pop" + str(thispop) + ".iSAFE.out"
	if not os.path.isfile(isafe_filename):
		isafe_filename += ".gz"
	

	component_files = [ihs_filename, delihh_filename, nsl_filename, xpehh_filename, freqs_filename, h12_filename, isafe_filename]
	#for item in component_files:
	#	if not os.path.isfile(item):
	#		item += ".gz" #not sure if this works 

	exists = [os.path.isfile(item) and os.path.getsize(item) > 500 for item in component_files]

	"""
	if os.path.isfile(freqs_filename):
		if os.path.getsize(freqs_filename) < 500: #QUICK HACK TO REMOVE MISMADE
			rm_cmd = "rm " + freqs_filename
			print(rm_cmd)
			subprocess.check_output(rm_cmd.split())
	"""
	
	if False in exists:
		missingstring = ""
		for ifile in range(len(exists)):
			if not exists[ifile]:
				missingstring += " " + component_files[ifile]
		print(missingstring)
		return []
	return component_files

def main():
	allargs = []
	#for model in models:
	#	for regime in regimes:
	#		for pop in pops:
	#			this_rundir = this_basedir + model +"_" + regime + "_sel" + str(pop) + "/"
	for subdir in subdirs:
		if True:
			for pop in [1]:
				reps = range(1,1001)
				this_rundir = this_basedir + subdir
				for irep in reps:
					tpedfilename = this_rundir + "tpeds/rep" + str(irep) + "_0_" + str(pop) + ".tped"
					if os.path.isfile(tpedfilename) or os.path.isfile(tpedfilename + ".gz"):
	
						outdir = this_rundir + "collated/"
						pairdir = this_rundir + "pairs/"
						#replicate_idstring = "rep" + str(irep) #+ "_" + str(pop)
						replicate_idstring = "rep" + str(irep) + "_sel1"#_1"
						tomake_filename = outdir + replicate_idstring + ".collated"
						if not os.path.isfile(tomake_filename) and not os.path.isfile(tomake_filename + ".gz"):

							pairfiles = []
							altpops = [1,2,3,4]#pops[:]
							altpops.remove(pop)
							for altpop in altpops:
								infostring = "_".join([str(x) for x in [pop, irep, altpop]])
								pairfilename = pairdir + infostring + ".pair"
								if not os.path.isfile(pairfilename) or os.path.getsize(pairfilename) == 0:
									components = get_remodel_components(this_rundir, replicate_idstring, altpop, pop)
									if len(components) > 0:
										ihs_filename, delihh_filename, nsl_filename, xpehh_filename, freqs_filename, h12_filename, isafe_filename = components
										pairfilename = pairdir + infostring + ".pair"
										combined_filename = nsl_filename + "\n" + h12_filename + "\n" + isafe_filename #QUICK HACK
										write_pair_sourcefile(pairfilename, ihs_filename, delihh_filename, combined_filename, xpehh_filename, freqs_filename)
								if os.path.isfile(pairfilename):
									pairfiles.append(pairfilename)
							if len(pairfiles) > 0:
								fullarg = tomake_filename
								for pairfile in pairfiles:
										fullarg += " " + pairfile
								fullcmd = basecmd + " " + fullarg
								allargs.append(fullarg)
								#print(fullcmd)
								#if not (regime == "hard" and irep == 111): #quick hack
								#subprocess.check_output(fullcmd.split())


	print('loaded ' + str(len(allargs)) + " cmds")
	####################
	## DISPATCH JOBS ###
	####################
	iscript = 0
	for argchunk in chunks(allargs,ncmds_script):
		iscript +=1
		scriptname = scriptdir + scriptprefix + "_" + str(iscript)
		writeJobsToBatch_fromArglist(basecmd, argchunk, scriptname, dispatch=dispatchNow, subopts = subopts)
		time.sleep(1)

main()
