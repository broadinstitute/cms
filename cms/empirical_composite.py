## /idi/sabeti-scratch/jvitti last used/updated 11.05.16

import os

runPrefix = "composite_emp_111616"
basedir = "/idi/sabeti-scratch/jvitti/scores_composite4_b/"
likesdir = '/idi/sabeti-scratch/jvitti/likes_111516_b/'

likessuffix = "neut" #linked
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/composite.py outgroups"
pwd = "/idi/sabeti-scratch/jvitti/cms/cms/"
models = ['nulldefault', 'default_112115_825am','default_default_101715_12pm','nulldefault_constantsize', 'gradient_101915_treebase_6_best',]
pops = ['YRI', 'CEU', 'CHB', 'BEB']
modelPops = [1, 2, 3, 4] 
chroms = range(1,23)

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
def get_pairfiles(chrom, selpop):
	altpops = pops[:]
	altpops.remove(selpop)
	inscorefilelist = []
	for altpop in altpops:
		filename = '/idi/sabeti-scratch/jvitti/synth/cms_pairs/' + str(selpop) + "_" + str(altpop) + "_chr" + str(chrom) + "_062215_2pm_nomultialleles_strictMask.cms"
		if os.path.isfile(filename):
			inscorefilelist.append(filename)
	return inscorefilelist
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
		for ipop in range(len(pops)):
			selpop, modelPop = pops[ipop], modelPops[ipop]
			hi_likesfile = get_likesfiles(model, modelPop, likessuffix, likesdir, allfreqs=True)
			mid_likesfile, low_likesfile = hi_likesfile, hi_likesfile

			for chrom in chroms:
				inscorefilelist = get_pairfiles(chrom, selpop) 
				if len(inscorefilelist) > 0:
					scorefilelist = ""
					for filename in inscorefilelist:
						scorefilelist += filename + ","
					scorefilelist = scorefilelist[:-1]
					outfile = basedir + selpop  +"_" + str(model) + "_" + likessuffix + ".chr" + str(chrom) + ".txt"
					if not os.path.isfile(outfile):
						argstring = scorefilelist + " " + hi_likesfile + " " + str(modelPop) + " --likesfile_low " + low_likesfile + " --likesfile_mid " + mid_likesfile + " " + outfile
						arguments.append(argstring)		

	writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix)

main()



