#fills missing replicate for neut
# last updated 11.16.16

import subprocess
import sys
import os
#model, pop, altpop, irep  = sys.argv[1:]

simRecomFile = "/idi/sabeti-scratch/jvitti/test_recom.recom"
def execute(commandstring):
	#print(commandstring)
	subprocess.check_output(commandstring.split())
	return

models = ['nulldefault', 'default_112115_825am','default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault_constantsize']#['nulldefault', 'default_112115_825am', 'default_default_101715_12pm','nulldefault_constantsize', 'gradient_101915_treebase_6_best',]
pops = [1, 2, 3, 4]


for model in models:
	for pop in pops:
		altpops = pops[:]
		altpops.remove(pop)
		for altpop in altpops:

			for irep in range(1, 1001):

				base_sim_dir = "/idi/sabeti-scratch/jvitti/sims/"
				base_score_dir = "/idi/sabeti-scratch/jvitti/scores/"
				model_score_dir, model_sim_dir = base_score_dir + model + "/", base_sim_dir + model + "/"
				neut_score_dir, neut_sim_dir = model_score_dir + "neut/", model_sim_dir + "neut/"
				tped = neut_sim_dir + "rep" + str(irep) + "_" + str(pop) + ".tped"
				assert(os.path.isfile(tped))

				toNormFile = neut_score_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) +".xpehh.out"
				normedFilename = neut_score_dir + "xpehh/rep" + str(irep) + "_" + str(pop) + "_" + str(altpop) + "_normed.txt"

				if os.path.isfile(toNormFile + ".norm") and not os.path.isfile(normedFilename):
					renamecmd = "mv " + toNormFile + ".norm " + normedFilename
					#print(renamecmd)
					execute(renamecmd)

