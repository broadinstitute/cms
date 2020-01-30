#fills missing replicate for neut
# last updated 11.16.16

import subprocess
import sys
import os
model, pop, altpop, irep  = sys.argv[1:]

base_sim_dir = "/idi/sabeti-scratch/jvitti/sims/"
base_score_dir = "/idi/sabeti-scratch/jvitti/scores/"
model_score_dir, model_sim_dir = base_score_dir + model + "/", base_sim_dir + model + "/"
neut_score_dir, neut_sim_dir = model_score_dir + "neut/", model_sim_dir + "neut/"
tped = neut_sim_dir + "rep" + str(irep) + "_" + str(pop) + ".tped"
assert(os.path.isfile(tped))

simRecomFile = "/idi/sabeti-scratch/jvitti/test_recom.recom"
def execute(commandstring):
	#print(commandstring)
	subprocess.check_output(commandstring.split())
	return


#ihs
ihs_outfilename = model_score_dir + "neut/ihs/rep" + str(irep) + "_" + str(pop) 
ihs_commandstring = "python /idi/sabeti-scratch/jvitti/cms/cms/scans.py selscan_ihs"
argstring = tped + " " + ihs_outfilename + " --threads 7 --truncOk"
fullcmd = ihs_commandstring + " " + argstring
if not os.path.isfile(ihs_outfilename + ".ihs.out"):
	execute(fullcmd)

#delihh
delihh_commandstring = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py scores_from_sims --delIhh"
ihs_infilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/ihs/rep" + str(irep) +"_" + str(pop) + ".ihs.out"
#print(ihs_infilename)
assert os.path.isfile(ihs_infilename)
delihhwritefile = model_score_dir + "neut/delihh/rep" + str(irep) + "_" + str(pop) + ".txt"
fullcmd = delihh_commandstring + " " + ihs_infilename + " "+ delihhwritefile
if not os.path.isfile(delihhwritefile):
	execute(fullcmd)

#ihs norm
neutNormParamFilename = neut_score_dir  + "ihs/n500_" +  str(pop) + "_concat.ihs.bins"
assert os.path.isfile(neutNormParamFilename)
normedfilename = "/idi/sabeti-scratch/jvitti/scores/" + model + "/neut/ihs/rep" + str(irep) +"_" + str(pop) +"_normed.txt"
argstring = "--normalizeIhs " + ihs_infilename + " " + normedfilename + " --neutIhsNormParams " + neutNormParamFilename
commandstring = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py scores_from_sims"
fullcmd = commandstring + " "+ argstring
if not os.path.isfile(normedfilename):
	execute(fullcmd)
	renamecmd = "mv " + ihs_infilename + ".norm " + normedfilename
	execute(renamecmd)

#delihh norm
neutNormParamFilename = neut_score_dir + "delihh/delihh_n500_" +  str(pop) + ".bins"
assert os.path.isfile(neutNormParamFilename)
assert os.path.isfile(delihhwritefile)
normedfile = delihhwritefile +".norm"
argstring = "--normalizeIhs " + delihhwritefile + " " + normedfilename + " --neutIhsNormParams " + neutNormParamFilename
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py scores_from_sims"
fullcmd = cmd + " " + argstring
if not os.path.isfile(normedfile):
	execute(fullcmd)

#xpehh
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/scans.py selscan_xpehh --threads 7"
tped2 = neut_sim_dir + "rep" + str(irep) + "_" + str(altpop) + ".tped"
outfilename = neut_score_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
argumentstring = tped + " " + outfilename + " " + tped2
fullcmd = cmd + " " + argumentstring
if not os.path.isfile(outfilename + ".xpehh.out"):
	execute(fullcmd)

#norm xp
neutNormParamFilename = neut_score_dir + "xpehh/n500_" +  str(pop) + "_" + str(altpop) + "_concat.xpehh.bins"
assert os.path.isfile(neutNormParamFilename)
toNormFile = neut_score_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) +".xpehh.out"
normedFilename = neut_score_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + "_normed.txt"
argstring = "--normalizeXpehh --neutXpehhNormParams " + neutNormParamFilename + " " + toNormFile + " " + normedFilename
commandstring = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py scores_from_sims"
fullcmd = commandstring + " " + argstring
#print(toNormFile)
#print(normedFilename)
if not os.path.isfile(normedFilename):
	execute(fullcmd)
	renamecmd = "mv " + toNormFile + ".norm " + normedFilename
	execute(renamecmd)

if os.path.isfile(toNormFile + ".norm"):
	renamecmd = "mv " + toNormFile + ".norm " + normedFilename
	execute(renamecmd)

#fstdeldaf
tped2 = neut_sim_dir + "rep" + str(irep) + "_" + str(altpop) + ".tped"
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/likes_from_model.py scores_from_sims"
outfilename = neut_score_dir + "fst_deldaf/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop)
argumentstring = tped + " --fst_deldaf " + tped2 + " " + outfilename + " --recomfile " + simRecomFile
fullcmd = cmd + " " + argumentstring
if not os.path.isfile(outfilename):
	execute(fullcmd)

