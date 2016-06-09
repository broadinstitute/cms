## {{POP GEN DATA --> DEM MODEL}}
## this script uses the specified binary of cosi 2.0 (Shlyakhter et al., 2014: https://github.com/broadinstitute/cosi2) to run the default point in
## parameter-space, using model topology/defaults as specified in params.py and error function as specified in error.py
## last updated: 06.04.16 	vitti@broadinstitute.org

cosi_build = "cosi_coalescent-2.0/coalescent" 
nReps, datestring = 300, "120915_245pm"
stats = ['pi', 'sfs', 'anc', 'r2', 'dprime', 'fst']
pops = [1, 2, 3, 4]

runprefix = "points/default_" #rename this variable
plotNow, plotPrefix = False, '/home/unix/vitti/public_html/dem_model/afr/trees/default_'

from params import *
from error import calcError
import subprocess
import sys 

#####################
## DEFINE FUNCTION ##
#####################

def default(nRep, datestring, stats = stats, pops = pops, plot=False):
	print "running default point in parameter space..."
	paramfilename = runprefix + str(datestring) + ".par"
	paramDict = generateParams()
	writeParamFile(paramfilename, paramDict)
	statfilename = runprefix + str(datestring) + ".stat"
	errorfilename = runprefix + str(datestring) + ".err"
	runstatscommand = cosi_build + " -p " + paramfilename + " -n " + str(nRep) + " --drop-singletons " + str(paramDict['singrate']) + "  --custom-stats --genmapRandomRegions > " + statfilename	
	print runstatscommand
	subprocess.call(runstatscommand, shell=True)
	print "calculating error with " + statfilename
	error = calcError(statfilename, stats, pops)
	print "error: " + str(error)
	errorfile = open(errorfilename, 'w')
	errorfile.write(str(error))
	errorfile.close()
	if plot == True:
		#print "plotting..."
		plot_comparison(statfilename, "default_" + str(datestring), nRep, stats)
		print "Drawing tree.."
		draw_tree(paramDict, plotPrefix + str(datestring) + ".png")
		draw_tree(paramDict, plotPrefix + str(datestring) + "_recent.png", recent=True)
		draw_tree(paramDict, plotPrefix + str(datestring) + "_veryrecent.png", veryrecent=True)
		print "Drawing migration matrix..."
		draw_mig_matrix(paramDict, plotPrefix + str(datestring) + "_mig.png")
	return

##########
## MAIN ##
##########

def main():
	default(nReps, datestring, stats, pops, plot=plotNow)

main()
