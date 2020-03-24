##	concat/normalizing neut-->bins for normalize sel--> etc.
##	01.13.2019	can naively cat files for ihs/nsl/delihh bc no header

sourcedir = "/idi/sabeti-scratch/jvitti/remodel/run/"
targetdir = "/idi/sabeti-scratch/jvitti/remodel/run/normalize_params/"
models, regimes, pops, reps = ['defdef15', 'gradient15', 'from_RC/nd', 'from_RC/ndcs'], ['soft', 'hard'], [1,2,3,4], range(1,5001)

import os
import subprocess

def main():
	all_sources = []
	for model in models:

		for regime in regimes:		
			for pop in pops:
				basedir = sourcedir + model + "_" + regime + "_sel" + str(pop) + "/xpehh/"
				all_sources.extend([basedir + item for item in os.listdir(basedir) if ".out" in item and "OUTGROUP" in item])
		print(model, str(len(all_sources)))
		
		for pop in pops:
			altpops = pops[:]
			altpops.remove(pop)
			for altpop in altpops:
				writefilename = targetdir + model + "_xpehh_OUTGROUP" + str(pop) + "_vs" + str(altpop) + ".txt"
				thisstring = str(pop) + "_vs" + str(altpop) 
				sourcefiles = [replicate for replicate in all_sources if thisstring in replicate]
				writefile = open(writefilename, 'w')
				for sourcefile in sourcefiles:
					readfile = open(sourcefile, 'r')
					readfile.readline() #strip header
					for line in readfile:
						writefile.write(line)
					readfile.close()
				writefile.close()
				print('wrote to ' + writefilename)

main()


#FREQS

#h12 for nd, ndc (AND DATA FIRST)