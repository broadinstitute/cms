## reads in reps from a pre-matched logfile and checks if outfile already exists
## /idi/sabeti-scratch/jvitti last used/updated 10.22.16

logPrefix = 'check_neut_102716am'#"recheck_neut_102516" 
runPrefix = "check_neut_poppairs_102916"

from run_likes_func import *

writedir = "/idi/sabeti-scratch/jvitti/scores_composite/"
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/composite.py poppair"
pwd = "/idi/sabeti-scratch/jvitti/"
models = ['nulldefault', 'default_112115_825am','default_default_101715_12pm', 'gradient_101915_treebase_6_best', 'nulldefault_constantsize']#['nulldefault', 'default_112115_825am', 'default_default_101715_12pm','nulldefault_constantsize', 'gradient_101915_treebase_6_best',]


def main():
	
	for model in models:
		modeldir = writedir + model + "/"
		check_create_dir(modeldir)
		neutdir = modeldir + "neut/"
		check_create_dir(neutdir)
		arguments = []
		approved_file = "/idi/sabeti-scratch/jvitti/logs/" + logPrefix + "_" + model + "_neut_match.log"
		openfile = open(approved_file, 'r')
		for line in openfile:
			entries = line.split()
			model, selpop, altpop, irep = entries
			in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_score_files(model, irep, selpop, altpop, "neut")
			outfile = neutdir + "rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
			if os.path.isfile(outfile) and os.path.getsize(outfile) > 0:
				pass
			else:
				argstring = in_ihs_file + " " + in_delihh_file + " " + in_xp_file + " " + in_fst_deldaf_file + " " + outfile
				arguments.append(argstring)
		openfile.close()
		writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix + "_" + model + "_neut", discardOut = True)


main()



