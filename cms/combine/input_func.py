## ancillary functions for composite.py  	
## last updated: 	06.19.2017	vitti@broadinstitute.org

import sys
import os

####################
## PREPARE INPUT ###
####################
def write_perpop_ihh_from_xp(infilename, outfilename, popNum = 1):
	outfile = open(outfilename, 'w')
	infile = open(infilename, 'r')
	infile.readline() #header
	for line in infile:
		entries = line.split()
		locus, pos, gpos, p1, ihh1, p2, ihh2, xpehh = entries

		if popNum ==1:
			writeline = pos + "\t" + gpos + "\t" + p1 + "\t" + ihh1 + "\n"
		elif popNum == 2:
			writeline = pos + "\t" + gpos + "\t" + p2 + "\t" + ihh2 + "\n"
		outfile.write(writeline)
	outfile.close()
	infile.close()
	print('wrote to: ' + outfilename)
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
	normalizedvalue = (rawscore - mean) #/ sd
	return normalizedvalue
