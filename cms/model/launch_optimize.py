## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

##############################################
## TOGGLE THIS RUN OF PIPELINE: DEFINE GRID ##
##############################################
if True:
	print "I should just put this in one spot and import it henceforth"
	gradientname = "fit_1_120715_3pm"
	gradloc = "gradients/"
	runDir = gradloc + gradientname + "/"

	#define dimensions of the objective function used for gradient descent
	pops = [1]
	stats = ['sfs', 'anc', 'pi']

	#define dimensions of search grid from which to launch gradient descent
	keys = [('Ne', 1), ('Ne', 1), ('Ne', 1), ('Ne', 1), ('Ne_times', 1), ('Ne_times', 1), ('Ne_times', 1),]
	indices = [0, 1, 2, 3, 1, 2, 3]
	assert len(keys) == len(indices)

	nCosiReps = 500 #how many simulations to run per point in parameter-space?
	paramsAsDists = False#False #encode parameters as distributions, or points?
	randomStart = False#False #start gradient at current optimum, or a random spot w.r.t params being searched?
	steprange = (.001, .5) 
	methods = ['L-BFGS-B', 'SLSQP'] 

###############################
## TOGGLE THIS RUN OF SCRIPT ##
############(i.e. w/in pipeline)

startrep = 0
numReplicates = 75
mkdir, launch, = True, True

from parse import givecommand
from itertools import product, izip
from random import choice, uniform
import subprocess
import time
import sys
import os

##########
## MAIN ##
##########

def main():

	###################
	## launch search ##
	###################

	if mkdir == True:
		mkdircommand = "mkdir " + runDir
		givecommand(mkdircommand)

	#encode as a string to facilitate parsing
	errorstring = "" 
	for stat in stats:
		errorstring += stat + "_"
	for pop in pops:
		errorstring += str(pop) + "_"

	globalargfilename = runDir + gradientname + "_globalargs.txt"
	globalargfile = open(globalargfilename, 'w')
	globalargfile.write(gradientname + "\n")
	globalargfile.write(errorstring + "\n")
	globalargfile.write(str(keys) + "\n")
	globalargfile.write(str(indices) + "\n")
	globalargfile.write(str(nCosiReps) + "\n")
	globalargfile.write(str(paramsAsDists) + "\n")
	globalargfile.write(str(randomStart) + '\n')
	globalargfile.close()

	for j in range(startrep, numReplicates):

		stepsize = uniform(steprange[0], steprange[1])
		method = choice(methods)

		scriptfilename = runDir + gradientname + "_" + str(j) + ".sh"
		errorfilename = runDir + gradientname + "_" + str(j) + ".e"
		outfilename = runDir + gradientname + "_" + str(j) + ".o"

		scriptfile= open(scriptfilename, 'w')
		scriptfile.write("#$ -cwd\n")
		scriptfile.write("python optimize.py " + str(stepsize) + " " + method + " " + globalargfilename)
		scriptfile.close()

		if launch == True:
			runscriptcommand = "qsub -o " + outfilename + " -e " + errorfilename + " " + scriptfilename
			givecommand(runscriptcommand)
		else:
			print runscriptcommand

main()