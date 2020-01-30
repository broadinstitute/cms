## /idi/sabeti-scratch/jvitti last used/updated 11.2.16

runPrefix = "composite_sims_111716"
basedir = "/idi/sabeti-scratch/jvitti/scores_composite_5/"
pairbasedir = "/idi/sabeti-scratch/jvitti/scores_composite/"
likesdir = '/idi/sabeti-scratch/jvitti/likes_111516_b/'

import os
import subprocess

likessuffix = "neut"
pwd = "/idi/sabeti-scratch/jvitti/cms/cms/"
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/composite.py outgroups"
models = ['default_112115_825am']#,'default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault_constantsize', 'nulldefault']
sel_freq_bins = ['0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90']
numPerBin = 1000
pops = [1, 2, 3, 4] 

def get_likesfiles(model, selpop, likessufix, likesdir, allfreqs = True):
	if not allfreqs:
		hi_likesfile = likesdir+ model + "/master/likes_" + str(selpop) + "_hi_vs_" + likessuffix + ".txt"
		mid_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_mid_vs_" + likessuffix + ".txt"
		low_likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_low_vs_" + likessuffix + ".txt"
		return hi_likesfile, mid_likesfile, low_likesfile 
	else:
		likesfile = likesdir + model + "/master/likes_" + str(selpop) + "_allfreq_vs_" + likessuffix + ".txt"
		return likesfile
	return
def get_inscorefiles(model, selpop, sel_freq_bin, irep, pairbasedir):
	inscorefiles = []
	altpops = pops[:]
	altpops.remove(selpop)
	for altpop in altpops:
		pairfilename = pairbasedir + model + "/sel" + str(selpop) + "/sel_" + str(sel_freq_bin) + "/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		if os.path.isfile(pairfilename):
			inscorefiles.append(pairfilename)
	return inscorefiles
def get_neutinscorefiles(model, selpop, irep, pairbasedir):
	inscorefiles = []
	altpops = pops[:]
	altpops.remove(selpop)
	for altpop in altpops:
		pairfilename = pairbasedir + model + "/neut/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		if os.path.isfile(pairfilename):
			inscorefiles.append(pairfilename)
	return inscorefiles
def check_create_dir(directory):
	if not os.path.isdir(directory):
		subprocess.check_output(['mkdir', directory])
def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-P sabeti_lab'], memory="5g", discardOut=False):
	"""flexible function for generating scripts and argfiles to dispatch parallel jobs on Univa Grid Engine.
	arglist: a list of strings. each string represents one argstring, for one task on the array."""
	#write file with arguments for all tasks
	t = 0
	taskargfilename = filewd + fileprefix + "_taskargs.txt"
	taskargfile = open(taskargfilename, "w")
	for argstring in arguments:
		taskargfile.write(argstring + "\n")
		t += 1
	taskargfile.close() 
	#write script for one job that points to task-argument-file
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("#$ -t 1-" + str(t) + "\n")
	scriptfile.write("\n")

	scriptfile.write('export PATH=/idi/sabeti-scratch/jvitti/cms/cms:$PATH\n')
	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/cms-env3\n')
	#scriptfile.write("source activate cms-env3\n")
	
	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	if discardOut:
		subcommand = "qsub -e /tmp/test_error -o /tmp/test_out -l h_vmem=" + str(memory) + " " + scriptname
	else:
		subcommand = "qsub -l h_vmem=" + str(memory) + " " + scriptname
	print(subcommand)
	return

def main():
	arguments = []
	for model in models:
		modelpop = basedir + model +'/'
		check_create_dir(modelpop)
		for selpop in pops:
			#arguments = []
			hi_likesfile = get_likesfiles(model, selpop, likessuffix, likesdir, allfreqs=True)
			mid_likesfile, low_likesfile = hi_likesfile, hi_likesfile

			writeseldir = modelpop + "/sel" + str(selpop)
			writeneutdir = modelpop + "/neut"			
			check_create_dir(writeseldir)
			check_create_dir(writeneutdir)


			#ALL SEL SIMS
			"""
			for sel_freq_bin in sel_freq_bins:
				writeselfreqbin = writeseldir + "/sel_" + str(sel_freq_bin)
				check_create_dir(writeselfreqbin)
				for irep in range(1, numPerBin +1):
					inscorefilelist = get_inscorefiles(model, selpop, sel_freq_bin, irep, pairbasedir) 
					if len(inscorefilelist) > 0:
						scorefilelist = ""
						for filename in inscorefilelist:
							scorefilelist += filename + ","
						scorefilelist = scorefilelist[:-1]
						outfile = basedir + model + "/sel" + str(selpop) + "/sel_" + str(sel_freq_bin) + "/rep" + str(irep) + "_" + str(selpop) +"_vsneut.cms.out"
						#if True:
						if os.path.isfile(outfile) and os.path.getsize(outfile) > 0:
							pass
						else:
							argstring = scorefilelist + " " + hi_likesfile + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + str(selpop) + " " + outfile 
							arguments.append(argstring)		
			"""
			#ALL NEUT
			#for irep in range(1, numPerBin +1):
			for irep in range(1, 100):#1):#range(500, 1001):	
				neutinscorefilelist = get_neutinscorefiles(model, selpop, irep, pairbasedir) 
				if len(neutinscorefilelist) > 0:
					scorefilelist = ""
					for filename in neutinscorefilelist:
						scorefilelist += filename + ","
					scorefilelist = scorefilelist[:-1]
					outfile = basedir  + model + "/neut/rep" + str(irep) + "_" + str(selpop) +".cms.out"
					
					if True:
					#if os.path.isfile(outfile) and os.path.getsize(outfile) > 0:
					#	pass
					#else:
						argstring = scorefilelist + " " + hi_likesfile + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + str(selpop) + " " + outfile
						arguments.append(argstring)

			#writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix + "_" + str(model) + "_sel" + str(selpop), discardOut=True)
	writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix, discardOut=False)

main()



