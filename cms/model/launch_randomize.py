## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

##############################################
## TOGGLE THIS RUN OF PIPELINE: DEFINE GRID ##
##############################################
if True:
	print "I should just put this in one spot and import it henceforth"
	singrate = .25
	cosi_build = "cosi_coalescent-2.0/coalescent" 
	nReps, datestring = 300, "120915_245pm"
	stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
	pops = [1, 2, 3, 4]
	numReplicates = 10
	skip = [] #in case don't want to redo all the replicates for this grid; 
	startRep, bigNumReps = 0, 500 #how many cosi simulations to run for corroborating/visualizing gradient output
			
###############################
## TOGGLE THIS RUN OF SCRIPT ##
############(i.e. w/in pipeline)
gradientprefix = "randomize_120915_430pm"
numSets = 20

from parse import givecommand
from params import *
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
	for i in range(numSets):
		gradientname = gradientprefix + "_" + str(i)
		numParams = 100
		startrep = 0
		numReplicates = 5
		mkdir, launch, = True, True

		pops = [1,2,3,4]
		stats = ['pi', 'sfs', 'anc', 'fst', 'r2', 'dprime']#, 'fst'

		nCosiReps = 500 #how many simulations to run per point in parameter-space?
		paramsAsDists = False#False #encode parameters as distributions, or points?
		randomStart = False#False #start gradient at current optimum, or a random spot w.r.t params being searched?

		steprange = (.01, .3) 
		methods = ['L-BFGS-B', 'SLSQP'] 

		queue = "long"
		maxSimulTests = 50
		dependencies = ['.python-2.7.1-sqlite3-rtrees',  '.numpy-1.8.0-python-2.7.1-sqlite3-rtrees', '.matplotlib-1.3.1-python-2.7.1-sqlite3-rtrees', '.scipy-0.15.1-python-2.7.1-sqlite3-rtrees']

		###################
		## launch search ##
		###################

		runDir = "gradients/" + gradientname + "/"
		#runDir = '/idi/sabeti-scratch/jvitti/dem_model/NRE_110915/gradients/' + gradientname + "/"

		if mkdir == True:
			mkdircommand = "mkdir " + runDir
			givecommand(mkdircommand)

		#choose params
		variableparams = generateParams()
		stable = ['chromlength', 'mutation_rate','num_indivs_per_sample','gene_conv_rel_rate', 'singrate','recomfilename','numPops', 'labels', ('Ne_change', 1), ('Ne_change', 2), ('Ne_change', 3), ('Ne_change', 4)]
		dependent = [
			[('Ne', 1), 0],
			[('Ne', 2), 0],
			[('Ne', 3), 0],
			[('Ne', 4), 0],
			[('bn_times', 2), 0],
			[('bn_times', 3), 0],
			[('bn_times', 4), 0],
			[('bn_times', 4), 1]
		]

		take_keys, take_indices = [], []
		take_index = 0 #initialize
		taken = []
		while len(take_keys) < numParams:
			#add a parameter
			key = choice(variableparams.keys())
			if key not in stable:
				numIndices = len(variableparams[key])
				take_index = choice(range(0,numIndices))
				key_ind = [key, take_index]
				if key_ind not in dependent and key_ind not in taken:
					#print "key: " 
					#print key
					#print "take index: " + str(take_index)
					take_keys.append(key)
					take_indices.append(take_index)
					taken.append(key_ind)
				#else:
					#print "passing on " 
					#print key_ind #for debug

		assert len(take_keys) == len(take_indices)
		keys, indices = take_keys, take_indices

		print keys
		print indices

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
			if broad == True:
				scriptfile.write("#$ -q " + queue + "\n")
				scriptfile.write("source /broad/software/scripts/useuse\n")
				for dependency in dependencies:
					scriptfile.write("reuse -q " + dependency  + "\n")
				scriptfile.write("python optimize.py " + str(stepsize) + " " + method + " " + globalargfilename)
			scriptfile.close()

			if launch == True:
				runscriptcommand = "qsub -o " + outfilename + " -e " + errorfilename + " " + scriptfilename
				givecommand(runscriptcommand)
			else:
				print runscriptcommand

main()
