## functions for composite
## last updated: 	02.10.2017	vitti@broadinstitute.org

import sys

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
	return ihs_hit_filename, ihs_miss_filename, nsl_hit_filename, nsl_miss_filename, delihh_hit_filename, delihh_miss_filename,  xpehh_hit_filename, xpehh_miss_filename, fst_hit_filename, fst_miss_filename, deldaf_hit_filename, deldaf_miss_filename
