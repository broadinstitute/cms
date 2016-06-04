## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

cosi_build = "cosi_coalescent-2.0/coalescent" 
gridloc = "grids/"

from parse import readGlobalArgfile
from error import calcError
from params import *
import subprocess
import sys

#####################
## DEFINE FUNCTION ##
####################

def samplePoint(nRep, scenario, values, keys, indices):
	"""values, keys, indices are as described in updateParams()"""
	print "running sample point in parameter space..."
	paramfilename = "points/" + scenario + ".par"
	paramDict = generateParams()
	paramDict_updated = updateParams(paramDict, keys, indices, values)
	writeParamFile(paramfilename, paramDict)
	statfilename = "points/" + scenario + ".stat"
	errorfilename = "points/" + scenario + ".err"
	runstatscommand = cosi_build + " -p " + paramfilename + " -n " + str(nRep) + " --drop-singletons " + str(paramDict['singrate']) + " --custom-stats --genmapRandomRegions > " + statfilename
	print runstatscommand
	subprocess.call(runstatscommand, shell=True)
	print "calculating error with " + statfilename
	error = calcError(statfilename, stats, pops)
	print "error: " + str(error)
	errorfile = open(errorfilename, 'w')
	errorfile.write(str(error))
	errorfile.close()
	return error

##########
## MAIN ##
##########

def main():

	globalargfilename = sys.argv[1]
	gridname, stats, pops, keys, indices, nreps, paramsAsDists, randomStart = readGlobalArgfile(globalargfilename)
	values = [float(x) for x in sys.argv[2:]]

	scenario = gridname
	for value in values:
		scenario += "_" + str(value)

	error = samplePoint(nreps, scenario, values, keys, indices)
	errorfile = gridloc + gridname + "/" + scenario + ".err"
	openfile = open(errorfile, 'w')
	openfile.write(str(error))
	openfile.close()

main()