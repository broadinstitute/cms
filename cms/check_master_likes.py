## last updated 11.4.16

basedir = "/idi/sabeti-scratch/jvitti/likes_110416/"

import subprocess
import os

models = ['nulldefault_constantsize', 'default_112115_825am', 'gradient_101915_treebase_6_best', 'nulldefault', 'default_default_101715_12pm']
selpops = [1, 2, 3, 4]
freqs = ['hi', 'low', 'mid']
misses = ["neut"]#, "linked"]

def main():
	for model in models:
		for selpop in selpops:
			for freq in freqs:
				for miss in misses:
					
					masterfiledir = basedir + model + "/master/"
					masterfilename = masterfiledir + "likes_" + str(selpop) + "_" + str(freq) + "_vs_" + str(miss) + ".txt"
					#print(masterfilename)
					openfile = open(masterfilename, 'r')

					for line in openfile:
						likesfile = line.strip('\n')
						if not os.path.isfile(likesfile) or os.path.getsize(likesfile) == 0:
							print("Missing likes file: " + likesfile)
						#else:
						#	print(likesfile)
					openfile.close()

main()
