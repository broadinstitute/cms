## functions to facilitate parallel computing
## last updated: 09.05.16 vitti@broadinstitute.org

def uger_array(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd']):
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
	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	subcommand = "qsub " + scriptname
	print(subcommand)
	return

def slurm_array(writefilename = "slurm_array.sh", arrayCmd = "./my_process $SLURM_ARRAY_TASK_ID", numTasks=1, dispatch=False, limit=5, mem=10, venv_command ="source activate cms-env3"):
	writefile = open(writefilename, 'w')
	writefile.write('#!/bin/bash\n')
	writefile.write('#SBATCH --array=1-' + str(numTasks) + "\n")
	if limit is not None: #in minutes
		writefile.write('#SBATCH --time=' + str(limit) + "\n")
	writefile.write('#SBATCH --mem=' + str(mem) + "\n")	
	writefile.write(venv_command + '\n')
	writefile.write(arrayCmd)
	writefile.write('\n')
	writefile.close()
	subcommand = "sbatch --array=1-" + str(numTasks) + " --time=" + str(limit) + " --mem=" + str(mem)
	subcommand += " " + writefilename
	print(subcommand)
	return
