## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING IN MAIN 
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

##############################################
## TOGGLE THIS RUN OF PIPELINE: DEFINE GRID ##
##############################################

cosi_build = "cosi_coalescent-2.0/coalescent" 
gridloc, gridname = "grids/", "fit1_120715_5pm"
keys = [('Ne', 1), ('Ne', 1), ('Ne', 1), ('Ne', 1), ('Ne_times', 1), ('Ne_times', 1), ('Ne_times', 1),]
indices = [0, 1, 2, 3, 1, 2, 3]
values = [[7500000], [40000], [7500, 8000, 8500], [15000, 18400, 25000], [50, 60, 70, 80, 90, 100], [350, 390], [400, 500, 600,]]
pops = [1, 2]
stats =['pi', 'sfs', 'anc', 'fst']
nCosiReps = 400 #how many simulations to run per point in parameter-space?
paramsAsDists = False#False #encode parameters as distributions, or points?
randomStart = False#False #start gradient at current optcat imum, or a random spot w.r.t params being searched?

###############################
## TOGGLE THIS RUN OF SCRIPT ##
############(i.e. w/in pipeline)

mode = "reportError" #"runGrid"

if mode == "runGrid":
	run, wait, reportError = True, False, False
elif mode == "reportError":
	run, wait, reportErr = False, False, True

from parse import givecommand
from itertools import product, izip
import sys
import os
import time

##########
## MAIN ##
##########

def main():
	#for error function 
	errorstring = "" #encode as a string to facilitate parsing
	for stat in stats:
		errorstring += stat + "_"
	for pop in pops:
		errorstring += str(pop) + "_"

	#record parameter with parameter value
	assert len(keys) == len(indices) and len(indices) == len(values)
	combos =  [' '.join(str(y) for y in x) for x in product(*values)]

	runDir = gridloc + gridname + "/"
	globalargfilename = runDir + gridname + "_globalargs.txt"
	argfilename = runDir + gridname + "_grid_args.txt"
	scriptname = runDir + gridname + "_grid.sh"
	outfilename = runDir + gridname + ".o"
	errorfilename = runDir + gridname + ".e"

	if run == True:
		mkdircommand = "mkdir " + runDir
		givecommand(mkdircommand)

		argfile = open(argfilename, 'w')
		numTests = 0
		for combo in combos:
			argstring = combo + "\n"
			argfile.write(argstring)
			numTests +=1
		argfile.close()

		globalargfile = open(globalargfilename, 'w')
		globalargfile.write(gridname + "\n")
		globalargfile.write(errorstring + "\n")
		globalargfile.write(str(keys) + "\n")
		globalargfile.write(str(indices) + "\n")
		globalargfile.write(str(nCosiReps) + "\n")
		globalargfile.write(str(paramsAsDists) + "\n")
		globalargfile.write(str(randomStart) + '\n')	
		globalargfile.close()

		scriptfile = open(scriptname, 'w')
		scriptfile.write("#$ -cwd\n")
		scriptfile.write("#$ -t 1-" + str(numTests) + "\n")
		#scriptfile.write("#$ -tc " + str(maxSimulTests) + "\n")
		if broad == True:
			scriptfile.write("source /broad/software/scripts/useuse\n")
			scriptfile.write("reuse -q .matplotlib-1.1.1rc-python-2.6.5\n")
			scriptfile.write("reuse -q .numpy-1.9.1-python-2.7.1-sqlite3-rtrees\n")
			scriptfile.write("reuse -q .scipy-0.15.1-python-2.7.1-sqlite3-rtrees\n")
			scriptfile.write("\n")
		scriptfile.write("MYFILE=" + str(argfilename) + "\n")
		scriptfile.write("SAMPLE=$(awk \"NR==$SGE_TASK_ID\" $MYFILE)\n")
		commandstring = "python grid_point.py " + globalargfilename + " $SAMPLE\n"
		scriptfile.write(commandstring)
		scriptfile.close()
		subcommand = "qsub " + scriptname 
		subprocess.call(subcommand, shell=True)

	errs, scenarrs = [], []
	if reportErr == True:
		logfilename = gridname + ".log"
		logfile = open(logfilename, 'w')
		for combo in combos:
			errorfilename = runDir + gridname 
			combo = combo.split(' ')
			for i in range(len(combo)):
				errorfilename += "_" + str(float(combo[i]))
			errorfilename += ".err"
			#print errorfilename
			
			if wait == True:
				while not os.path.exists(errorfilename):
					time.sleep(1)

			if os.path.isfile(errorfilename):
				errorfile = open(errorfilename, 'r')
				error = float(errorfile.readline())
				errs.append(error)
				scenarrs.append(combo)
				errorfile.close()
			else:
				print "missing: " + errorfilename
		sorted_errors = sorted(izip(errs, scenarrs)) 
		for error in sorted_errors:
			print error
			logfile.write(str(error))

main()