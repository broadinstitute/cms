##	vitti@broadinstitute.org
##	last updated 10.06.2017

import subprocess

def chunks(l, n):
	for i in range(0, len(l), n):
		yield  l[i:i+n]
def check_create_dir(dirname):
	mkdircmd = "mkdir -p " + dirname
	subprocess.check_output( mkdircmd.split() )
	return
def writeJobsToBatch_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-P sabeti_lab', '-l h_rt=48:00:00'], 
activate_venv = True, dispatch = False, subopts = ''): #re-add arg for dependencies
	""" writes a script for all of these tasks to be performed sequentially """

	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("\n")
	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	if activate_venv:
		scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/cms-env3\n')
	scriptfile.write('source /broad/software/scripts/useuse\nreuse BEDTools\nreuse VCFtools\n')
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


def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd', '-P sabeti_lab'], activate_venv = True, 
dispatch = False, subopts = ''): #re-add arg for dependencies
	"""flexible function for generating scripts and argfiles to dispatch parallel jobs on Univa Grid Engine.
	arglist: a list of strings. each string represents one argstring, for one task on the array."""

	##########################################
	#write file with arguments for all tasks #
	##########################################
	t = 0
	taskargfilename = filewd + fileprefix + "_taskargs.txt"
	taskargfile = open(taskargfilename, "w")
	for argstring in arguments:
		taskargfile.write(argstring + "\n")
		t += 1
	taskargfile.close() 

	#############################################################
	#write script for one job that points to task-argument-file #
	#############################################################
	scriptname = filewd + fileprefix + '.sh'
	scriptfile = open(scriptname, 'w')
	for UGERoption in UGERoptions:
		scriptfile.write("#$ " + UGERoption + "\n")
	scriptfile.write("#$ -t 1-" + str(t) + "\n")
	scriptfile.write("\n")
	scriptfile.write('export PATH=/home/unix/vitti/miniconda3/bin:$PATH\n')
	if activate_venv:
		scriptfile.write('source activate /home/unix/vitti/miniconda3/envs/cms-env3\n')
	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	subcommand = "qsub " + subopts + scriptname
	print(subcommand)
	if dispatch:
		subprocess.check_output( subcommand.split() )
	return subcommand

