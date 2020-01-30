# ad hoc script to write fst mean and var from concat file for norm
# 11.23.16

rewriteBinfile, normFromBinfile = False, True

models = ['nulldefault', 'gradient_101915_treebase_6_best', 'default_default_101715_12pm', 'default_112115_825am','nulldefault_constantsize',]
pops = [1, 2, 3, 4]

from parse_func import get_concat_files
import numpy as np
import os

def rewrite_fst_binfile_from_concat(concatfilename, binfilename):
	vals = []
	concatfile = open(concatfilename, 'r')
	header = concatfile.readline()
	for line in concatfile:
		entries = line.split()
		fst = float(entries[1])
		vals.append(fst)
	concatfile.close()
	numVals = len(vals) 
	print('loaded ' + str(numVals) + " fst values...")
	mean = np.mean(vals)
	var = np.var(vals)
	binfile = open(binfilename, 'w')
	binfile.write('num\tmean\tvariance\n')
	binfile.write(str(numVals)+ "\t" + str(mean) + "\t" + str(var) + '\n')
	binfile.close()
	print("wrote to : "+ binfilename)
	return
def load_normparams_from_binfile(binfilename):
	openfile = open(binfilename, 'r')
	openfile.readline()
	dataline = openfile.readline()
	openfile.close()
	entries = dataline.split()
	numVals, mean, var = int(entries[0]), float(entries[1]), float(entries[2])
	return numVals, mean, var
def norm_fst_file(repfilename, mean, var):
	sd = np.sqrt(var)
	assert os.path.isfile(repfilename)
	writefilename = repfilename + ".norm"
	writefile = open(writefilename, 'w')
	readfile = open(repfilename, 'r')
	headerline = readfile.readline()
	newheaderline = headerline.strip('\n') + "\tfst_norm\n"
	writefile.write(newheaderline)
	for line in readfile:
		entries = line.split()
		this_unnormed_fst = float(entries[1])
		normed_fst = (this_unnormed_fst - mean) / var
		writeline = line.strip('\n') + '\t' + str(normed_fst) + '\n'
		writefile.write(writeline)
	writefile.close()
	readfile.close()
	return


def main():
	for model in models:
		for selpop in pops:
			altpops = pops[:]
			altpops.remove(selpop)
			for altpop in altpops:
				concatfilename, binfilename = get_concat_files(model, selpop, 'fst', altpop)
				if rewriteBinfile:
					rewrite_fst_binfile_from_concat(concatfilename, binfilename)
				if normFromBinfile:
					numVals, mean, var = load_normparams_from_binfile(binfilename)
					for irep in range(1, 3):#1001):
						infilename = "/idi/sabeti-scratch/jvitti/clean/scores/" + model + "/neut/fst_deldaf/rep" + str(irep) + "_" + str(selpop) + "_" + str(altpop)
						norm_fst_file(infilename, mean, var)

main()
