## functions for composite  	
## last updated: 	04.10.2017	vitti@broadinstitute.org

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

################
## LOAD INPUT ##
################
def get_likesfiles_frommaster(masterfile, modelselpop, vs="neut"):
	"""assumes only one demographic model in master likes file, but there may be multiple pops or neut/linked comps"""
	likesfiles = []
	openfile = open(masterfile, 'r')
	for line in openfile:
		filename = line.strip('\n')
		likesfiles.append(filename)
	openfile.close()

	print('master likesfile has a total of ' + str(len(likesfiles)) + " files..")
	#print(likesfiles)
	selpop_likesfiles = [item for item in likesfiles if "sel" + str(modelselpop) in item]
	#print(selpop_likesfiles)
	ihsfiles = [item for item in selpop_likesfiles if "/ihs/" in item]
	nslfiles = [item for item in selpop_likesfiles if "/nsl/" in item]
	delihhfiles = [item for item in selpop_likesfiles if "/delihh/" in item]
	xpehhfiles = [item for item in selpop_likesfiles if "/xpehh/" in item]
	fstfiles = [item for item in selpop_likesfiles if "/fst/" in item]
	deldaffiles = [item for item in selpop_likesfiles if "/deldaf/" in item]

	ihs_selfiles = [item for item in ihsfiles if "causal" in item]
	nsl_selfiles = [item for item in nslfiles if "causal" in item]
	delihh_selfiles = [item for item in delihhfiles if "causal" in item]
	xpehh_selfiles = [item for item in xpehhfiles if "causal" in item]
	fst_selfiles = [item for item in fstfiles if "causal" in item]
	deldaf_selfiles = [item for item in deldaffiles if "causal" in item]

	ihs_neutfiles = [item for item in ihsfiles if vs in item]
	nsl_neutfiles = [item for item in nslfiles if vs in item]
	delihh_neutfiles = [item for item in delihhfiles if vs in item]
	xpehh_neutfiles = [item for item in xpehhfiles if vs in item]
	fst_neutfiles = [item for item in fstfiles if vs in item]
	deldaf_neutfiles = [item for item in deldaffiles if vs in item]

	for filelist in [ihs_selfiles, ihs_neutfiles, nsl_selfiles, nsl_neutfiles, delihh_selfiles, delihh_neutfiles, xpehh_selfiles, xpehh_neutfiles, fst_selfiles, fst_neutfiles, deldaf_selfiles, deldaf_neutfiles]:
		#print(filelist)
		if len(filelist) >1:
			print('conflicting distributions found...')
			sys.exit()
		elif len(filelist) == 0:
			print('missing distributions for population ' + str(modelselpop))
			sys.exit()

	else:
		ihs_hit_filename, ihs_miss_filename = ihs_selfiles[0], ihs_neutfiles[0]
		nsl_hit_filename, nsl_miss_filename = nsl_selfiles[0], nsl_neutfiles[0]
		delihh_hit_filename, delihh_miss_filename = delihh_selfiles[0], delihh_neutfiles[0]
		xpehh_hit_filename, xpehh_miss_filename = xpehh_selfiles[0], xpehh_neutfiles[0]	
		fst_hit_filename, fst_miss_filename = fst_selfiles[0], fst_neutfiles[0]	
		deldaf_hit_filename, deldaf_miss_filename = deldaf_selfiles[0], deldaf_neutfiles[0]
	return ihs_hit_filename, ihs_miss_filename, nsl_hit_filename, nsl_miss_filename, delihh_hit_filename, delihh_miss_filename,  xpehh_hit_filename, xpehh_miss_filename, fst_hit_filename, fst_miss_filename, deldaf_hit_filename, deldaf_miss_filename #HMMM

#REDUNDANT WITH PARSE_FUNC; fix.
def get_emp_cms_file(selpop, chrom, normed = False, basedir = "/n/regal/sabeti_lab/jvitti/clear-synth/1kg_scores/", suffix = ".model_a"): #CONNECT "MODEL" TO "RUNSUFFIX"
	""" locates CMS files for empirical data """
	#filename = basedir + "chr" + str(chrom) + "_" + str(selpop) + "_strictMask_" + model + ".cms" + suffix
	filename = basedir + "composite/chr" + str(chrom) + "_" + str(selpop) + ".cms.out" + suffix
	if normed:
		filename += ".norm"
	if not os.path.isfile(filename):
		print("MISSING empirical file : " + filename)
	return filename
def load_empscores(model, selpop, normed = False,takeIndex = -1, basedir = "/n/regal/sabeti_lab/jvitti/clear-synth/1kg_scores/", suffix = ".model_a"):
	""" for genome-wide empirical CMS data, loads a value according to takeIndex """	
	chroms = range(1,23)
	scores = []
	for chrom in chroms:
		scorefile = get_emp_cms_file(selpop, chrom, normed=normed, basedir=basedir, suffix=suffix,) #model
		print('loading from ' + scorefile)
		assert os.path.isfile(scorefile)
		openfile = open(scorefile, 'r')
		for line in openfile:
			entries=line.split()
			scores.append(float(entries[takeIndex]))
		openfile.close()
	return scores


#THIS PERHAPS MOVES? collapse others.
def normalize(rawscore, mean, sd):
	rawscore, mean, sd = float(rawscore), float(mean), float(sd)
	normalizedvalue = (rawscore - mean) / sd
	return normalizedvalue
