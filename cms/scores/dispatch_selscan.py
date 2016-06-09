## {{POP GEN DATA --> COMPONENT SCORES}}
## original: selscan/genloc_fixed/dispatch_selscan_cluster.py 12.30.15
# JV: should call taskarg function properly and automate poppairs
# particular to 1KG. this is a pretty basic and useless script, should probably just replace
## last updated 06.05.16 	vitti@broadinstitute.org

#############################################
## TOGGLE OPTIONS FOR THIS RUN OF PIPELINE ##
#############################################

prefix = "122915"
chroms = range(1,22)
pops = ['ESN', 'GWD', 'LWK', 'MSL',  'YRI', 'ACB', 'ASW', 'PUR', 'CLM', 'MXL', 'PEL','CDX', 'CHB', 'CHS',  'JPT', 'KHV','CEU', 'FIN', 'GBR',  'IBS', 'TSI','BEB', 'GIH', 'ITU',  'PJL', 'STU']
popPairs = [('CHS', 'CHB'), ('JPT', 'ESN'), ('LWK', 'FIN'), ('MSL', 'GBR'), ('MXL', 'PJL'), ('PEL', 'STU'), ('PUR', 'GIH'), ('CEU', 'BEB'), ('ACB', 'IBS'), ('ASW', 'TSI'), ('CDX', 'YRI'), ('CHS', 'CHB'), ('CLM', 'ITU'), ('GWD', 'KHV')]
extraArgs = ['--threads 35'] #e.g. --gapscale, --maf
mask, vers, runtime, suffix = 'strictMask', 'v5a', '062215_2pm_nomultialleles', 'full_xchr_122915'

############################################
## TOGGLE OPTIONS FOR THIS RUN OF SCRIPT  ##
##############(i.e., within run of pipeline)

convertTped = False
runIhs, runNsl, runXpehh = True, False, False
normIhs, normNsl, normXpehh = False, False, False
vcfloc, tpedloc, ihsloc, nslloc, xpehhloc = '/DATA/vcfs/', '/DATA/tpeds/', '/DATA/ihs/', '/DATA/nsl/', '/DATA/xpehh/'

import os
from func_scores import *

############
### MAIN ###
############

def main():
	####################
	### WRITE SCRIPT ###
	####################
	for chrom in chroms:
		if chrom == "X":
			vers = 'v1a'

		##################
		### PER - POP  ###
		##################
		if runIhs == True or runNsl == True or normIhs == True or normNsl == True or convertTped == True:
			for pop in pops:

				if convertTped == True:
					inputVCF = vcfloc + vers + "_chr" + str(chrom) + "_" + mask + "_global.vcf.gz" 
					genMap = '/DATA/genmaps/genetic_map_GRCh37_chr' + str(chrom) + ".txt"
					outlocation = tpedloc
					outprefix = vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + suffix
					scriptname = "convert_" + str(chrom) + "_"  + pop + "_" + prefix + ".sh"
					scriptfile = open(scriptname, 'w')
					argstring = "selscan_file_conversion"
					for arg in extraArgs:
						argstring += " " + arg
					argstring += " " + inputVCF + " " + genMap + " " + outprefix + " " + outlocation + " " + str(chrom)
					argstring += " --sampleMembershipFile /DATA/integrated_call_samples_v3.20130502.ALL.panel --filterPops " + pop + "\n"

					scriptfile.write("source /DATA/mount_venv/bin/activate\n")
					commandstring = "python2.7 /DATA/mount_venv/cms/cms/scans.py " + argstring
					scriptfile.write(commandstring)
					scriptfile.close()
					subcommand = "qsub " + scriptname
					print subcommand

				inputTped = tpedloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + suffix + ".tped.gz"
				if runIhs == True:
					scriptname = "ihs_" + str(chrom) + "_"  + pop + "_" + prefix + ".sh"
					scriptfile = open(scriptname, 'w')
					outFile = ihsloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + suffix + "_ihs"
					taskstring = "selscan_ihs"
					for arg in extraArgs:
						taskstring += " " + arg
					taskstring += " " + inputTped
					taskstring += " " + outFile
					scriptfile.write("source /DATA/mount_venv/bin/activate\n")
					commandstring = "python2.7 /DATA/mount_venv/cms/cms/scans.py " + taskstring
					scriptfile.write(commandstring)
					scriptfile.close()
					subcommand = "qsub " + scriptname
					print subcommand

				if runNsl == True:
					scriptname = "nsl_" + str(chrom) + "_"  + pop + "_" + prefix + ".sh"
					scriptfile = open(scriptname, 'w')
					outFile = nslloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + suffix + "_nsl"
					taskstring = "selscan_nsl"
					for arg in extraArgs:
						taskstring += " " + arg
					taskstring += " " + inputTped
					taskstring += " " + outFile
					scriptfile.write("source /DATA/mount_venv/bin/activate\n")
					commandstring = "python2.7 /DATA/mount_venv/cms/cms/scans.py " + taskstring
					scriptfile.write(commandstring)
					scriptfile.close()
					subcommand = "qsub " + scriptname
					print subcommand

		#######################
		### PER - POP PAIR  ###
		#######################
		elif runXpehh == True or normXpehh == True:
			for popPair in popPairs:
				pop1, pop2 = popPair[0], popPair[1]
				inputTped1 = tpedloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop1 + "_" + runtime + "_" + suffix + ".tped.gz"
				inputTped2 = tpedloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop2 + "_" + runtime + "_" + suffix + ".tped.gz"
				if runXpehh == True:
					scriptname = "xpehh_" + str(chrom) + "_"  + pop1 + "_" + pop2 + "_" + prefix + ".sh"
					scriptfile = open(scriptname, 'w')
					outFile = xpehhloc + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop1 + "_" + pop2 + "_" + runtime + "_" + suffix + "_xpehh"
					taskstring = "selscan_xpehh"
					for arg in extraArgs:
						taskstring += " " + arg
					taskstring += " " + inputTped1
					taskstring += " " + outFile
					taskstring += " " + inputTped2
					scriptfile.write("source /DATA/mount_venv/bin/activate\n")
					commandstring = "python2.7 /DATA/mount_venv/cms/cms/scans.py " + taskstring
					scriptfile.write(commandstring)
					scriptfile.close()
					subcommand = "qsub " + scriptname
					print subcommand

main()