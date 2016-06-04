## {{POP GEN DATA --> DEM MODEL}}
## this script uses the specified binary of cosi 2.0 (Shlyakhter et al., 2014: https://github.com/broadinstitute/cosi2) 
## and the parameter adjustments specified herewith (relative to default: params.py)
## last updated: 06.04.16 	vitti@broadinstitute.org

cosi_build = "/idi/sabeti-scratch/ilya/cosi/builds/bnew/coalescent_joe"
scenario,nRep = "perturb_120715_430pm", 400
reusePre, pre_stat = True,  "points/default_120715_350pm.stat"
pops = [1, 2, 3, 4]
stats = ['pi', 'sfs', 'anc', 'fst'] #'r2', 'dprime', 'fst',]

#######################
## ADJUST PARAMETERS ##
#######################

keys = [('split', 2)]
indices = [0]
values = [250]

############
## IMPORT ## 	#JV: gratuitous?
############

from params import *
from plot import *
import subprocess
from error import *
from draw_tree import *

#####################
## DEFINE FUNCTION ##
#####################

def samplePoint_perturb(nRep, scenario, values, keys, indices, stats, pops, reuse=False, reusestatfile = ""):
	"""runs a modified scenario and compares to default."""
	rundir = "/idi/sabeti-scratch/jvitti/dem_model/afr/points/"
	paramfilenamepre = rundir + scenario + "_pre.par"
	statfilenamepre = rundir + scenario + "_pre.stat"
	errorfilenamepre = rundir + scenario + "_pre.err"
	paramfilenamepost = rundir + scenario + "_post.par"
	statfilenamepost = rundir + scenario + "_post.stat"
	errorfilenamepre = rundir + scenario + "_post.err"
	paramDict = generateParams()
	if reuse == True:
		print "calculating error with " + reusestatfile
		error = calcError(reusestatfile, stats, pops)
	else:
		writeParamFile(paramfilenamepre, paramDict)
		unstatscommand = cosi_build + " -p " + paramfilenamepre + " -n " + str(nRep) + " --drop-singletons " + str(paramDict['singrate']) + "  --custom-stats --genmapRandomRegions > " + statfilenamepre
		print runstatscommand
		subprocess.call(runstatscommand, shell=True)
		error=calcError(statfilenamepre, stats, pops)
	print "PRE error: " + str(error)
	
	paramDict_updated = updateParams(paramDict, keys, indices, values)
	writeParamFile(paramfilenamepost, paramDict_updated)
	runstatscommand = cosi_build + " -p " + paramfilenamepost + " -n " + str(nRep) + " --drop-singletons " + str(paramDict_updated['singrate']) + "  --custom-stats --genmapRandomRegions > " + statfilenamepost	
	print runstatscommand
	subprocess.call(runstatscommand, shell=True)
	error=calcError(statfilenamepost, stats, pops)
	print "POST error: " + str(error)
	if reuse == True:
		plot_before_and_after_comparison(reusestatfile, statfilenamepost, scenario)
	else:
		plot_before_and_after_comparison(statfilenamepre, statfilenamepost, scenario)
	draw_tree(paramDict_updated,'/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + ".png", recent=False)
	draw_tree(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + "_recent.png", recent=True)
	draw_tree(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + "_veryrecent.png", veryrecent=True)
	draw_mig_matrix(paramDict_updated, '/home/unix/vitti/public_html/dem_model/afr/trees/' + scenario + "_mig.png")
	return
	
##########
## MAIN ##
##########

def main():
	samplePoint_perturb(nRep, scenario, values, keys, indices, stats, pops, reusePre, pre_stat)

main()