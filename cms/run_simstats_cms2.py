## reads in reps from a pre-matched logfile and checks if outfile already exists
## /idi/sabeti-scratch/jvitti last used/updated 10.22.16

logPrefix = "check_102216"
runPrefix = "dispatch_102216" 
#"cms/cms/dispatch_102216"

from run_likes_func import *
import sys


pops = [1,2,3,4]

writedir = "/idi/sabeti-scratch/jvitti/scores_composite/"
cmd = "python /idi/sabeti-scratch/jvitti/cms/cms/composite.py poppair"
pwd = "/idi/sabeti-scratch/jvitti/"

def main():
	arguments = []

	model, pop = sys.argv[1], sys.argv[2]
	approved_file = "/idi/sabeti-scratch/jvitti/" + logPrefix + "_"+ model + "_sel" + str(pop) + "_match.log"
	openfile = open(approved_file, 'r')
	for line in openfile:
		entries = line.split()
		model, selpop, altpop, sel_freq_bin, irep = entries
		#if selpop != pop:
		#	print(selpop)
		#	print(pop)
		assert(selpop == pop)
		in_ihs_file, in_delihh_file, in_xp_file, in_fst_deldaf_file = get_score_files(model, irep, selpop, altpop, sel_freq_bin)
		outfile = writedir + model + "/sel" + str(selpop) + "/sel_" + str(sel_freq_bin) + "/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop) + ".pair"
		if os.path.isfile(outfile):
			pass
		else:
			argstring = in_ihs_file + " " + in_delihh_file + " " + in_xp_file + " " + in_fst_deldaf_file + " " + outfile
			arguments.append(argstring)
	openfile.close()
	writeJobsToTaskArray_fromArglist(cmd, arguments, runPrefix + "_" + model + "_" + str(pop), discardOut = True)

main()

