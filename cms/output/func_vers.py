# /idi/sabeti-scratch/jvitti/synth
# 01.13.16

def writeJobsToTaskArray_fromArglist(commandstring, arguments, fileprefix, filewd = "", UGERoptions = ['-cwd'], dependencies = ['.python-2.7.8-sqlite3-rtrees', '.zlib-1.2.6', '.gsl-1.14', '.matplotlib-1.1.1rc-python-2.6.5', '.numpy-1.8.0-python-2.7.1-sqlite3-rtrees',
	'.samtools-1.2', '.tabix-0.2.6', '.biopython-1.64-python-2.7.1-sqlite3-rtrees']):
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

	scriptfile.write("source /broad/software/scripts/useuse\n")
	for dependency in dependencies:
		scriptfile.write("reuse -q " + dependency  + "\n")

	scriptfile.write("MYFILE=" + str(taskargfilename) + "\n")
	scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
	writestring = commandstring + " $SAMPLE\n"
	scriptfile.write(writestring)
	scriptfile.close()
	subcommand = "qsub " + scriptname
	print subcommand
	return

def givecommand(commandstr):
	commands = commandstr.split(' ')
	subprocess.check_output(commands)