## last updated 11.15.16

basedir = "/idi/sabeti-scratch/jvitti/likes_111516_b/"

import subprocess
import os

def write_master_likesfile(writefilename, model, selpop, freq,basedir,  miss = "neut",):
	'''adapted from run_likes_func.py'''
	writefile = open(writefilename, 'w')
	for score in ['ihs', 'delihh']: 
		hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_causal.txt"#_smoothed.txt"
		misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_" + str(freq) + "_" + miss + ".txt"#"_smoothed.txt"
		#assert(os.path.isfile(hitlikesfilename) and os.path.isfile(misslikesfilename))
		writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	for score in ['xpehh', 'fst', 'deldaf']:
		hitlikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_causal.txt"#_smoothed.txt"
		misslikesfilename = basedir + model + "/" + score + "/likes_sel" + str(selpop) + "_choose_" + str(freq) + "_" + miss + ".txt"#"_smoothed.txt"
		#assert(os.path.isfile(hitlikesfilename) and os.path.isfile(misslikesfilename))
		writefile.write(hitlikesfilename + "\n" + misslikesfilename + "\n")
	writefile.close()
	print("wrote to: " + writefilename)
	return

models = ['nulldefault_constantsize', 'default_112115_825am', 'gradient_101915_treebase_6_best', 'nulldefault', 'default_default_101715_12pm']
selpops = [1, 2, 3, 4]
freqs = ['allfreq']#'hi', 'low', 'mid']
misses = ["neut"]#, "linked"]

def main():
	for model in models:
		for selpop in selpops:
			for freq in freqs:
				for miss in misses:
					
					writefiledir = basedir + model + "/master/"
					if not os.path.isdir(writefiledir):
						command = 'mkdir ' + writefiledir
						subprocess.check_output(command.split())
					writefilename = writefiledir + "likes_" + str(selpop) + "_" + str(freq) + "_vs_" + str(miss) + ".txt"
					write_master_likesfile(writefilename, model, selpop, freq, basedir, miss)

main()
