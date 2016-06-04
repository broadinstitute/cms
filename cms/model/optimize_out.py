## {{POP GEN DATA --> DEM MODEL}}
## this script NEEDS A LOT OF FINE-TOOTHED CLEANING
## currently configured to fit 1KG P3 African populations: {YRI, ESN, MSL, GWD, LWK}
## last updated: 06.04.16 	vitti@broadinstitute.org

##############################################
## TOGGLE THIS RUN OF PIPELINE: DEFINE GRID ##
##############################################
if True:
	print "I should just put this in one spot and import it"
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

plotPipe, plotComparison, runAndPlotFinal = False, False, False
#plotPrefix = '/home/unix/vitti/public_html/dem_model/afr/trees/default_'

from draw_tree import *
from params import *
from error import *
from parse import *
from plot import *
import subprocess
import sys

gradientname = sys.argv[1]

######################
## DEFINE FUNCTIONS ##
######################

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
	draw_mig_matrix(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/'  + scenario +  "_mig.png")

	return error
def getGlobalErr(scenario):
	#in case we already ran the point 
	fileprefix = "points/" + scenario 
	errorfilename = fileprefix + ".err"
	openfile = open(errorfilename, 'r')
	error = float(openfile.readline())
	openfile.close()
	return error

##########
## MAIN ##
##########

def main():

	writefilename = gradientname + ".log"
	writefile = open(writefilename, 'w')
	localErrs, globalErrs = {}, {} #store replicate number as entry in dict
	if broad == True:
		runDir = '/idi/sabeti-scratch/jvitti/dem_model/afr/gradients/' + gradientname + "/"
	else:
		runDir = 'gradients/' + gradientname + "/"
	for i in range(startRep, numReplicates):
		if i not in skip:
			print "**************"
			print str(i)
			writefile.write('**************\n' + str(i)+'\n')
			gradientrep = gradientname + "_" + str(i)
			errorfilename = runDir + gradientrep + ".e"
			outfilename = runDir + gradientrep + ".o"

			#check to see that this gradient replicate worked
			#print errorfilename
			if readOutfile(outfilename) != False:
		#	if not checkFailed(errorfilename) and readOutfile(outfilename) != False:
				print outfilename
				writefile.write(outfilename + "\n")
				if os.path.isfile(outfilename) and sum(1 for line in open(outfilename, 'r')) > 0:
					pointParamss, errors, numIters, errorDict, keys, indices, descentmethod, stepsize, beststat, bestparam, firststat, firstparam, modifications, pops, stats = readOutfile(outfilename)
					#print modifications
					if i == 0 or i == numReplicates:
						print str(pops)
						writefile.write(str(pops) + "\t")			
						print str(stats)
						writefile.write(str(stats) + "\n")	
					
					writefile.write("num iters: " + str(numIters) + "\n")
					writefile.write("stepsize: " + str(stepsize) + "\n")
					writefile.write("descent method: " + str(descentmethod))
					writefile.write("local error: " + str(min(errors)) + "\n")
					print "local error: " + str(min(errors)) + "\n"
					localErrs[min(errors)] = gradientrep
					
					#show us what happened in the pipeline
					if plotPipe == True:
						print "\tPIPELINE:"
						plotPipeline(errors, errorDict, pointParamss, keys, indices, gradientrep, descentmethod, stepsize)

					#before and after: pop stats
					if plotComparison == True:
						print"\tBEFORE AND AFTER:"
						plot_before_and_after_comparison(firststat, beststat, gradientrep)

					#get the suggested modifications, make them, and reassess their global effect on the model.
					if runAndPlotFinal == False:
						errorfilename = gradientname + "_" + str(i) + "_best.err"
						if os.path.isfile(errorfilename):
							error = getGlobalErr(errorfilename)
							globalErrs[error] = gradientrep
							writefile.write("global error: " + str(error) + "\n")
							print "global error: " + str(error) + "\n" 
						else:
							print "must run final point"
					else:
						print"\tFINAL, nreps: " + str(bigNumReps)
						values, keys, indices = [], [], []
						for item in modifications:
							entries = item.split()
							if 'presentSizes' in entries[0] or 'admix_time' in entries[0] or 'admix_rate' in entries[0]:
								key = entries[0]
							else:
								key = eval(entries[0])
							index, value = eval(entries[1]), eval(entries[2])
							keys.append(key)
							indices.append(index)
							values.append(value)
						error = samplePointAndPlot(bigNumReps, gradientname + "_" + str(i) + "_best", values, keys, indices)
						print "global error: " + str(error) +'\n'
						writefile.write("global error: " + str(error) + "\n")
						globalErrs[error] = gradientrep

					for item in modifications:
						writefile.write("\t" + item + "\n")

			else:
				print "gradient failed; see: " + errorfilename

	#compare replicates, print best
	globalKeys = globalErrs.keys()
	localKeys = localErrs.keys()
	minGlobal = min(globalKeys)
	minLocal = min(localKeys)
	globalRep = globalErrs[minGlobal]
	localRep = localErrs[minLocal]
	writestring = "lowest local error " + str(minLocal) + " achieved in gradient # " + str(localRep) + "\nlowest global error " + str(minGlobal) + " achieved in gradient # " + str(globalRep)
	print writestring
	writefile.write(writestring)
	writefile.close()

main()