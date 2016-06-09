## top-level script for plotting haplotypes. Adapted from Shervin's hapviz; original files in CMS_1KG_P3/viz/ 
## /idi/sabeti-scratch/jvitti/synth 01.13.16 updated 05.27.16 to pull from YRI_IPSC-filtered VCFs.

writeNow, plotNow = True, False #toggle options for this run of script
useTaskArrays = True
dataSources = ['yri_ipsc']#["filteredtped", "maskedvcf", "unmaskedvcf"]  

towriteprefix = "FADS_locus_052716_b"
chrom, startPos, endPos = 11, 61324000, 61354000 
corePos, snp_pos = False, 31095171
#corePos, snp_pos = True, 16141151
pops = ['YRI']#, 'CEU']#'CHB', 'CEU','GWD', 'LWK', 'MSL', 'ESN']#, 'ACB', 'ASW']
#pops = ['GWD', 'LWK', 'MSL',  'ACB', 'ASW', 'PUR', 'CLM', 'MXL', 'PEL','CDX',  'CHS',  'JPT', 'KHV','FIN', 'GBR',  'IBS', 'TSI', 'GIH', 'ITU',  'PJL', 'STU', 'YRI', 'BEB', 'CEU', 'CHB', 'ESN']

vers, mask, runtime, mergesuffix = 'v5a', 'strictMask', '062215_2pm_nomultialleles', "merged_121915"
scriptDir = "/idi/sabeti-scratch/jvitti/synth/scripts/viz/"
writeDir = "/idi/sabeti-scratch/jvitti/synth/hapviz/"
plotDir, pubPlotDir = "/home/unix/vitti/public_html/synth/hapviz/", "http://www.broadinstitute.org/~vitti/synth/hapviz/"

import sys
from func_vers import *

############
### MAIN ###
############

if writeNow == True:
	################################
	## PULL REGION; WRITE TO FILE ##
	################################
	argstrings = []
	for dataSource in dataSources:
		for pop in pops:
			if dataSource == 'unmaskedvcf':
				sourcefilename = "/idi/sabeti-scratch/jvitti/cms_venv/vcfs/ALL.chr" + str(chrom) + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
			elif dataSource == 'maskedvcf':
				sourcefilename = "/idi/sabeti-scratch/jvitti/cms_venv/vcfs/" + vers + "_chr" + str(chrom) + "_" + mask + "_global.vcf.gz"
			elif dataSource == "filteredtped":
				sourcefilename = "/idi/sabeti-scratch/jvitti/cms_venv/tpeds/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + ".tped.gz"
			elif dataSource == "yri_ipsc":			
				if pop == 'YRI':
					sourcefilename = "/idi/sabeti-scratch/jvitti/yri_ipsc/YRI_ipsc_indivs_chr" + str(chrom) + ".recode.vcf"
				else:
					pass
			argstring = dataSource + " " + sourcefilename + " " + str(chrom) + " " + str(startPos) + " " + str(endPos) + " " + pop + " " + towriteprefix
			argstrings.append(argstring)
	commandstring = "python pullRegion.py"
	if useTaskArrays == True:
		writeJobsToTaskArray_fromArglist(commandstring, argstrings, towriteprefix + "_pullregions", scriptDir,  UGERoptions = ['-cwd', '-q long', '-P sabeti_lab', '-M vitti@broadinstitute.org -m eas'])
	else:
		for item in argstrings:
			fullcommand = commandstring + " " + item
			print fullcommand
			#givecommand(fullcommand)

if plotNow == True:
	#############################
	## SORT HAPLOTYPES AND PLOT##
	#############################
	argstrings = []
	for dataSource in dataSources:
		for pop in pops:
			hapfilename = writeDir + towriteprefix + "_" + pop + "_" + dataSource + ".hapz"
			imagesavename = towriteprefix + "_chr" + str(chrom) + "_" + str(startPos) + "_" + str(endPos) + "_" + pop + "_" + dataSource + ".png"
			fullimageloc = plotDir + imagesavename
			argstring = "-haps " + hapfilename + " -out " + fullimageloc 
			if corePos == True:
				argstring += " -corepos chr" + str(chrom) + "." + str(snp_pos)
			argstrings.append(argstring)
	commandstring = "python shervin_hapviz.py"
	if useTaskArrays == True:
		writeJobsToTaskArray_fromArglist(commandstring, argstrings, towriteprefix + "_plothaps", scriptDir)#, UGERoptions = ['-cwd', '-q long', '-P sabeti_lab', '-M vitti@broadinstitute.org -m eas'])
	else:
		for item in argstrings:
			fullcommand = commandstring + " " + item
			print fullcommand	
