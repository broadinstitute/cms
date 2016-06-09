## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## CALL PROGRAM AS: python optimize.py <stepsize> <descentmethod> <config file>
## last updated: 06.04.16 	vitti@broadinstitute.org

cosi_build = "cosi_coalescent-2.0/coalescent" 
nReps, datestring = 300, "120915_245pm"
stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
pops = [1, 2, 3, 4]

plotNow, plotPrefix = False, '/home/unix/vitti/public_html/dem_model/afr/trees/default_'

from parse import readGlobalArgfile, giveCommand
from error import calcError
from params import *
from scipy import optimize
import numpy as np
import subprocess
import string
import random
import sys
import os

if len(sys.argv) != 4:
	print "call program as: python optimize.py <stepsize> <descentmethod> <config file>"
	sys.exit()

stepsize = float(sys.argv[1])
descentmethod = sys.argv[2]
configfilename = sys.argv[3]
gradientname, stats, pops, keys, indices, nreps, paramsAsDists, randomStart = readGlobalArgfile(configfilename)

######################
## DEFINE FUNCTIONS ##
######################

def samplePoint(nRep, scenario, values, keys, indices, paramsAsDists = False, plot=False):
	"""values, keys, indices are as described in updateParams()"""
	print "running sample point in parameter space..."
	print values
	sys.stderr.write(str(values))
	sys.stderr.write("\n")
	ranges = getRanges()
	actualValues = []
	for i in range(len(keys)):
		key = keys[i]
		index = indices[i]
		value = values[i]
		interval = rangeDict[key][index]
		low, high = float(interval[0]), float(interval[1])
		actual = getRealValue(value, low, high)
		actualValues.append(actual)
	sys.stderr.write(str(actualValues))
	sys.stderr.write("\n")

	savestring = ''.join(random.choice(string.ascii_uppercase) for i in range(12))
	fileprefix = "gradients/" + scenario + "/" + scenario + "_" + savestring

	#for item in actualValues:
	#if len(str(item)) > 6:
		#item = str(item)[:6]
		#fileprefix += "_" + str(item)
	#if len(fileprefix) > 250:
		#fileprefix = fileprefix[:250] 
	paramfilename = fileprefix + ".par"
	statfilename = fileprefix +  ".stat"
	errorfilename = fileprefix + ".err"
	paramDict = generateParams()
	paramDict_updated = updateParams(paramDict, keys, indices, actualValues)
	if paramsAsDists == False:
		writeParamFile(paramfilename, paramDict_updated)
	else:
		distDict = paramsToDists(paramDict_updated)
		writeParamFile(paramfilename, distDict)

	singrate = paramDict_updated['singrate']
	#runstatscommand = cosi_build + " -p " + paramfilename + " -n " + str(nRep) + " --drop-singletons " + str(singrate) + "  --custom-stats > " + statfilename	
	runstatscommand = cosi_build + " -p " + paramfilename + " -n " + str(nRep) + " --drop-singletons " + str(singrate) + "  --custom-stats --genmapRandomRegions > " + statfilename

	print runstatscommand
	subprocess.call(runstatscommand, shell=True)
	print "calculating error with " + statfilename
	if os.path.isfile(statfilename) and os.path.getsize(statfilename) > 0:
		error = calcError(statfilename, stats, pops)
	else:
		error = 100
	print "error: " + str(error)
	sys.stderr.write(str(error) + "\n")
	errorfile = open(errorfilename, 'w')
	errorfile.write(str(error))
	errorfile.close()
	return error
def samplePoint_wrapper(values):
	return samplePoint(nreps, gradientname, values, keys, indices, paramsAsDists)

#############
### MAIN ###
#############

def main():

	rangeDict = getRanges()
	paramDict = generateParams()
	x0 = []
	bounds = []
	for i in range(len(keys)):
		key = keys[i]
		index = indices[i]
		print str(key)
		print str(index)
		value = paramDict[key][index]
		interval = rangeDict[key][index]
		low, high = float(interval[0]), float(interval[1])
		scaled = getScaledValue(value, low, high)
		if randomStart == False:
			x0.append(scaled)
		else:
			randomval = random.random()
			x0.append(randomval)
		#line = "Low: " + str(low) + " high: " + str(high) + " scaled: " + str(scaled) + " actual: " + str(value) + "\n"
		#sys.stderr.write(line)
		bounds.append([0,1])

	x0 = np.array(x0)

	stepdict = {'eps':float(stepsize)}

	print "START: "
	print descentmethod
	print str(stepsize)
	print stats
	print pops
	print keys
	print indices
	print x0
	print " ---------------- "

	result = optimize.minimize(samplePoint_wrapper, x0, method=descentmethod, bounds=bounds, options=stepdict)
	print result

	print "******************"
	#translate back to changes to model
	bestparams = []
	assert len(keys) == len(result.x)
	for i in range(len(keys)):
		key = keys[i]
		index = indices[i]
		interval = rangeDict[key][index]
		low, high = float(interval[0]), float(interval[1])
		realVal = getRealValue(result.x[i], low, high)
		bestparams.append(result.x[i])
		print "best " + str(key) + "|" + str(index) + "|" + str(realVal)
main()