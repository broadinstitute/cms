## 11.16.16 quick script to run component scores -> poppair files for extra neutral reps (#500-1000) for simulate data.
## (previously used check_reps_neut.py -> run_simstats_cms2_neut.py, but I dont need to check concordance any more)

writedir = "/idi/sabeti-scratch/jvitti/scores_composite/" #just for poppairs
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/composite.py poppair"
pwd = "/idi/sabeti-scratch/jvitti/"
models = ['nulldefault', 'default_112115_825am','default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault_constantsize']#['nulldefault', 'default_112115_825am', 'default_default_101715_12pm','nulldefault_constantsize', 'gradient_101915_treebase_6_best',]
pops = [1, 2, 3, 4]

from power_func import *


def get_score_files(model, repNum, selPop, altPop, scenario = "neut", filebase = "/idi/sabeti-scratch/jvitti/scores/"):
	if scenario == "neut":
		in_ihs_file = filebase + model + "/neut/ihs/rep" + str(repNum) + "_" + str(selPop) + "_normed.txt"
		in_delihh_file = filebase + model + "/neut/delihh/rep" + str(repNum) +  "_" + str(selPop) + ".txt.norm"
		in_xp_file = filebase + model + "/neut/xpehh/rep" + str(repNum) +  "_" + str(selPop) + "_" + str(altPop) + "_normed.txt"
		in_fst_deldaf_file = filebase + model + "/neut/fst_deldaf/rep" + str(repNum) + "_" + str(selPop) + "_" + str(altPop)

	else:
		in_ihs_file = filebase + model + "/sel" + str(selPop) + "/ihs/sel_" + str(scenario) + "/rep" + str(repNum) + "_" + str(selPop) + "_truncOk.ihs.out.norm"
		in_delihh_file = filebase + model + "/sel" + str(selPop) + '/delihh/sel_' + str(scenario) + '/rep' + str(repNum) + "_" + str(selPop) + ".txt.norm"
		in_xp_file = filebase + model + "/sel" + str(selPop) + '/xpehh/sel_' + str(scenario) +'/rep' + str(repNum) + "_" + str(selPop) + "_" + str(altPop) + ".xpehh.out.norm"
		in_fst_deldaf_file = filebase + model + "/sel" + str(selPop) + '/fst_deldaf/sel_' + str(scenario) +'/rep' + str(repNum) + "_" + str(selPop) + "_" + str(altPop)

	for returnfile in [in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file]:
		if not os.path.isfile(returnfile):
			print("Missing: " + returnfile)

	return in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file 

def main():
	
	for model in models:
		modeldir = writedir + model + "/"
		check_create_dir(modeldir)
		neutdir = modeldir + "neut/"
		check_create_dir(neutdir)
		arguments = []
		for selpop in pops:
			altpops = pops[:]
			altpops.remove(selpop)
			for altpop in altpops:
				for irep in range(1, 1001):
					in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_score_files(model, irep, selpop, altpop, "neut")
					outfile = neutdir + "rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
					if os.path.isfile(outfile) and os.path.getsize(outfile) > 0:
						pass
					else:
						argstring = in_ihs_file + " " + in_delihh_file + " " + in_xp_file + " " + in_fst_deldaf_file + " " + outfile
						arguments.append(argstring)
		writeJobsToTaskArray_fromArglist(cmd, arguments, "runadditional" + "_" + model + "_neut", discardOut = True)


main()



