## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

scenarioname = "test_121015_1pm" 
nRep = 500
pointnames = ['test.par'] 
runPoint, printMergedChanges = True, True

runDir = "/idi/sabeti-scratch/jvitti/dem_model/AWS/gradients/"
cosi_build = "/idi/sabeti-scratch/ilya/cosi/builds/bnew/coalescent_joe"
singrate = .25

from plot import plot_comparison
from parse import readOutfile
from error import calcError
from draw_tree import *
from params import *
import subprocess
import sys
import os

#####################
## DEFINE FUNCTION ##
#####################

def samplePointAndPlot(nRep, scenario, values, keys, indices):
	"""values, keys, indices are as described in updateParams()"""
	print "running sample point in parameter space..."
	print keys
	print indices
	print values
	fileprefix = "points/" + scenario 
	paramfilename = fileprefix + ".par"
	statfilename = fileprefix +  ".stat"
	errorfilename = fileprefix + ".err"
	paramDict = generateParams()
	paramDict_updated = updateParams(paramDict, keys, indices, values)
	writeParamFile(paramfilename, paramDict_updated)

	runstatscommand = cosi_build + " -p " + paramfilename + " -n " + str(nRep) + " --drop-singletons " + str(paramDict['singrate']) + "  --custom-stats --genmapRandomRegions > " + statfilename
	print runstatscommand
	subprocess.call(runstatscommand, shell=True)
	print "calculating error with " + statfilename
	error = calcError(statfilename)
	print "error: " + str(error)
	sys.stderr.write(str(error) + "\n")
	errorfile = open(errorfilename, 'w')
	errorfile.write(str(error))
	errorfile.close()

	print "plotting..."
	plot_comparison(statfilename, scenario, nRep)
	print "Drawing tree.."
	draw_tree(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + ".png")
	draw_tree(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + "_recent.png", recent=True)
	draw_tree(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + "_veryrecent.png", veryrecent=True)
	print "Drawing migration matrix..."
	draw_mig_matrix(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/'  + scenario +  "_mig.png")
	return error

##############
###  MAIN ####
##############

def main():
	## Get info from outfiles
	allMods = [] #same param given 2x will overwrite e.o. downstream
	for pointname in pointnames:
		gradientname = ""
		for item in pointname.split('_')[:-1]:
			gradientname += item + "_"
		outfilename = runDir + gradientname[:-1] + "/" + pointname +  ".o"
		if os.path.isfile(outfilename) and sum(1 for line in open(outfilename, 'r')) > 0:
			pointParamss, errors, numIters, errorDict, keys, indices, descentmethod, stepsize, beststat, bestparam, firststat, firstparam, modifications, pops, stats = readOutfile(outfilename)
			allMods.extend(modifications)
		else:
			print "eh? not a point: " + outfilename

	values, keys, indices = [], [], []
	for item in allMods:
		entries = item.split()
		if 'presentSizes' in entries[0] or 'admix_time' in entries[0] or 'admix_rate' in entries[0]:
			key = entries[0]
		else:
			key = eval(entries[0])
		index, value = eval(entries[1]), eval(entries[2])
		keys.append(key)
		indices.append(index)
		values.append(value)

	# Run point with modifications prescribed by specified outfiles
	if runPoint == True:
		error = samplePointAndPlot(nRep, scenarioname, values, keys, indices)

	# Print updated param arrays in an easy-to-copy format to allow changing params_exp.py
	if printMergedChanges == True:
		paramDict = generateParams()
		paramDict_updated = updateParams(paramDict, keys, indices, values)
		printed = []
		for key in keys:
			if key not in printed:
				print key
				print paramDict_updated[key]
				printed.append(key)

main()