## ancillary functions for composite.py  	
## last updated: 	09.14.2017	vitti@broadinstitute.org

import sys
import os
import math
import numpy as np
import gzip

####################
## PREPARE INPUT ###
####################
def write_perpop_ihh_from_xp(infilename, outfilename, popNum = 1):

	outfile = open(outfilename, 'w')

	#### First iteration - ensure proper handling of file.
	if ".gz" in infilename:
		infile = gzip.open(infilename, 'rt')
	else:
		infile = open(infilename, 'r')
	infile.readline() #header
	for line in infile:
		entries = line.split()
		xpehh = float(entries[7])
		entry1, entry2, entry3, entry4 = [float(item) for item in entries[3:7]] #there are inconsistencies
		### need a flexible solution that can figure out which column is which
		if entry3 == 0 or entry4 == 0:
			pass
		else:
			testXp = np.log(float(entry2)/float(entry4))
			testXp2 = np.log(float(entry1)/float(entry3))
			if math.fabs(testXp - float(xpehh)) < .001:
				p1_ind, ihh1_ind, p2_ind, ihh2_ind = 3, 4, 5, 6
			elif math.fabs(testXp2 - float(xpehh)) < .001:
				p1_ind, ihh1_ind, p2_ind, ihh2_ind = 4, 3, 6, 5
			else:
				print(line)
				print((str(testXp)))
				print((str(testXp2)))
				print(('check indices ' + infilename))
			break
	infile.close()
	if ".gz" in infilename:
		infile = gzip.open(infilename, 'rt')
	else:
		infile = open(infilename, 'r')
	infile.readline() #header
	for line in infile:
		entries = line.split()
		xpehh = float(entries[7])
		locus = entries[0]
		pos = entries[1]
		gpos = entries[2]
		p1 = entries[p1_ind]
		p2 = entries[p2_ind]
		ihh1 = entries[ihh1_ind]
		ihh2 = entries[ihh2_ind]

		if popNum ==1:
			writeline = pos + "\t" + gpos + "\t" + p1 + "\t" + ihh1 + "\n"
		elif popNum == 2:
			writeline = pos + "\t" + gpos + "\t" + p2 + "\t" + ihh2 + "\n"
		outfile.write(writeline)
	outfile.close()
	infile.close()
	print(('wrote to: ' + outfilename))
	return outfilename
def write_pair_sourcefile(writefilename, ihsfilename, delihhfilename, nslfilename, xpehhfilename, freqsfilename):
	#if not os.path.isfile(writefilename):
	if True:
		openfile = open(writefilename, 'w')
		openfile.write(ihsfilename+ "\n")
		openfile.write(delihhfilename+ "\n")
		openfile.write(nslfilename+ "\n")
		openfile.write(xpehhfilename+ "\n")
		openfile.write(freqsfilename+ "\n")
		openfile.close()
	return writefilename
def write_run_paramfile(writefilename, ihs_master_likesfile, nsl_master_likesfile, delihh_master_likesfile, xpehh_master_likesfile,
	fst_master_likesfile, deldaf_master_likesfile, cutoffline, includeline):
	#if not os.path.isfile(writefilename):
	if True:
		openfile = open(writefilename, 'w')
		openfile.write(ihs_master_likesfile + "\n")
		openfile.write(nsl_master_likesfile + "\n") #CHANGE ORDER
		openfile.write(delihh_master_likesfile + "\n")
		openfile.write(xpehh_master_likesfile + "\n")
		openfile.write(fst_master_likesfile + "\n")
		openfile.write(deldaf_master_likesfile + "\n")	
		openfile.write(cutoffline + "\n")
		openfile.write(includeline + "\n")		
		openfile.close()
	return writefilename
def normalize(rawscore, mean, sd):
	rawscore, mean, sd = float(rawscore), float(mean), float(sd)
	normalizedvalue = (rawscore - mean) / sd
	return normalizedvalue
