##	run new, heterogenous demography model set; generate component scores
##	02.21.2019

cmsdir = "/idi/sabeti-scratch/jvitti/cms/cms/"
writedir = "/idi/sabeti-scratch/jvitti/remodel/run3/"
simRecomFile =  "/idi/sabeti-scratch/jvitti/params/test_recom.recom"
pops = [1, 2, 3, 4]

import subprocess
import sys
import os

def execute(action):
	subprocess.check_output(action.split())
	return

def main():
	if len(sys.argv) != 6:				#'selpop' here is actually neutral outgroup in 'usepop' sel scenario
		print('call program as: python remodel_components.py {model} {regime} {selpop} {irep} {usepop}')
		sys.exit(0)
	model = sys.argv[1]
	regime = sys.argv[2]
	selpop = int(sys.argv[3])
	irep = int(sys.argv[4])
	usepop = int(sys.argv[5])

	basedir = writedir + model + "_" + regime + "/" #+ "_sel" + str(usepop) + "/"
	tped_dir = basedir + "tpeds/"
	thispop = selpop

	replicate_idstring = "rep" + str(irep) + "_pop" + str(selpop)#+ "_OUTGROUP" + str(selpop)
	
	tped_filename = tped_dir + "rep" + str(irep) + "_0_" + str(selpop) + ".tped"
	if not os.path.isfile(tped_filename):
		tped_filename += ".gz"					
	if not os.path.isfile(tped_filename):
		print(('missing: ', tped_filename))	
		sys.exit(0)				

	"""
	ihh12_commandstring = "python " + cmsdir + "scans.py selscan_ihh12" 
	ihh12_unnormedfileprefix = basedir + "ihh12/" + replicate_idstring
	ihh12_argstring = tped_filename + " " + ihh12_unnormedfileprefix + " --threads 7 "
	ihh12_fullcmd = ihh12_commandstring + " " + ihh12_argstring
	ihh12_unnormedfilename = ihh12_unnormedfileprefix + ".ihh12.out"
	print(ihh12_fullcmd)
	if not os.path.isfile(ihh12_unnormedfilename)  and not os.path.isfile(ihh12_unnormedfilename + ".gz"):
		execute(ihh12_fullcmd)


	"""
	ihs_commandstring = "python " + cmsdir + "scans.py selscan_ihs"
	ihs_outfileprefix = basedir + "ihs/" + replicate_idstring
	ihs_unnormedfile = ihs_outfileprefix + ".ihs.out"
	ihs_argstring = tped_filename + " " + ihs_outfileprefix + " --threads 7 "
	ihs_fullcmd = ihs_commandstring + " " + ihs_argstring
	#ihs_normedfile = ihs_unnormedfile + ".norm"
	#print(ihs_fullcmd)
	#if not os.path.isfile(ihs_unnormedfile) and not os.path.isfile(ihs_unnormedfile + ".gz"):
	#	execute(ihs_fullcmd)
	
	delihh_commandstring = "python " + cmsdir + "composite.py delihh_from_ihs"
	delihh_unnormedfile =  basedir + "delihh/" + replicate_idstring
	delihh_argstring = ihs_unnormedfile + " "+ delihh_unnormedfile
	delihh_fullcmd = delihh_commandstring + " " + delihh_argstring 
	delihh_normedfile = delihh_unnormedfile + ".norm"
	print(delihh_fullcmd)
	if not os.path.isfile(delihh_unnormedfile) and not os.path.isfile(delihh_unnormedfile + ".gz"):
		execute(delihh_fullcmd)		
	"""
	nsl_commandstring = "python " + cmsdir + "scans.py selscan_nsl" 
	nsl_unnormedfileprefix = basedir + "nsl/" + replicate_idstring
	nsl_argstring = tped_filename + " " + nsl_unnormedfileprefix + " --threads 7 "
	nsl_fullcmd = nsl_commandstring + " " + nsl_argstring
	nsl_unnormedfilename = nsl_unnormedfileprefix + ".nsl.out"
	print(nsl_fullcmd)
	if not os.path.isfile(nsl_unnormedfilename)  and not os.path.isfile(nsl_unnormedfilename + ".gz"):
		execute(nsl_fullcmd)


	tpeddir = tped_dir
	
	altpops = pops[:]
	altpops.remove(int(thispop))
	altpops.remove(int(usepop))
	for altpop in altpops:
		xpehh_commandstring = "python " + cmsdir + "scans.py selscan_xpehh --threads 7"
		#tped2 = tpeddir + "rep" + str(irep) + "_" + str(altpop) + ".tped"
		#tped_filename2 = get_tped_filename(selpop, irep, ancestralpop, altpop, model, tped_dir)
		tped_filename2 = tped_dir + "rep" + str(irep) + "_0_" + str(altpop) + ".tped"
		if not os.path.isfile(tped_filename2):
			tped_filename2 += ".gz"					
	
		xpehh_outfileprefix = basedir + "xpehh/" + replicate_idstring + "_vs" + str(altpop)
		xpehh_unnormedfile = basedir + "xpehh/" + replicate_idstring + "_vs" + str(altpop) + ".xpehh.out"
		xpehh_argumentstring = tped_filename + " " + xpehh_outfileprefix + " " + tped_filename2 + " --threads 7 "
		xpehh_fullcmd = xpehh_commandstring + " " + xpehh_argumentstring
		print(xpehh_fullcmd)
		if not os.path.isfile(xpehh_unnormedfile) and not os.path.isfile(xpehh_unnormedfile + ".gz"):
			execute(xpehh_fullcmd)

	
		#fstdeldaf_commandstring = "python " + cmsdir + "composite.py freqscores"
		#fstdeldaf_outfilename = basedir + "freqs/"  + replicate_idstring + "_vs" + str(altpop)
		#fstdeldaf_argumentstring = tped_filename + " " + tped_filename2 + " " + simRecomFile + " " + fstdeldaf_outfilename 
		#fstdeldaf_fullcmd = fstdeldaf_commandstring + " " + fstdeldaf_argumentstring 
		#print(fstdeldaf_fullcmd)
		#if not os.path.isfile(fstdeldaf_outfilename):
		#	execute(fstdeldaf_fullcmd)
	"""
main()
