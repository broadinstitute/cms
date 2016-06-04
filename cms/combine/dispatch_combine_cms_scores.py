## {{COMPONENT SCORES + SCORE DISTRIBUTIONS FROM SIMULATED DATA --> COMPOSITE SCORES}}
## last updated: 06.04.15   vitti@broadinstitute.org

prefix = "combine_all_pops_011416_720pm"
pops = ['GWD', 'LWK', 'MSL',  'ACB', 'ASW', 'PUR', 'CLM', 'MXL', 'PEL','CDX',  'CHS',  'JPT', 'KHV','FIN', 'GBR',  'IBS', 'TSI', 'GIH', 'ITU',  'PJL', 'STU', 'YRI', 'BEB', 'CEU', 'CHB', 'ESN']
selPops = ['YRI']
chroms = range(1,23)
UGERoptions = ['-cwd', '-q long', '-P sabeti_lab']
#from func_files import *
#from func_vers import *
import os

comparePairs = False #necessary intermediate step
compareMultiplePops = True #more powerful and flexible

############
### MAIN ###
############

def main():
	cDict = {True:'0', False:'1'}

	###############################
	#### 2-WAY POP COMPARISONS ####
	###############################

	if comparePairs == True:
		popPairs = []
		for ipop in range(len(pops)):
			for jpop in range(ipop, len(pops)):
				if ipop != jpop:
					popPairs.append((pops[ipop], pops[jpop]))

		args = []
		for chrom in chroms:
			for popPair in popPairs:		
				argstring = ""
				pop1, pop2 = popPair[0], popPair[1]
				ihs1filename, ihs2filename, delihh1filename, delihh2filename, xpfilename, fst_deldaffilename, rev_xp, rev_fstdeldaf = getFileNames(pop1, pop2, chrom)

				argstring1 = str(chrom) + "\t" + pop1 + "\t" + pop2 + "\t" + ihs1filename + "\t" + delihh1filename + "\t" + xpfilename + "\t" + str(cDict[rev_xp]) + "\t" + fst_deldaffilename + "\t" + str(cDict[rev_fstdeldaf])
				argstring2 = str(chrom) + "\t" + pop2 + "\t" + pop1 + "\t" + ihs2filename + "\t" + delihh2filename + "\t" + xpfilename + "\t" + str(cDict[not rev_xp]) + "\t" + fst_deldaffilename + "\t" + str(cDict[not rev_fstdeldaf])
				args.append(argstring1) 
				args.append(argstring2) 

		commandstring = "./combine_cms_scores_poppairs"
		writeJobsToTaskArray_fromArglist(commandstring, args, prefix + "_pairs", UGERoptions=UGERoptions)

	###########################
	#### GROUP COMPARISONS ####
	###########################

	if compareMultiplePops == True:
		args = []
		for selPop in selPops:
			otherPops = pops[:]
			otherPops.remove(selPop)
			argstring = str(chrom) + "\t" + selPop + "\t"
			for otherPop in otherPops:
				argstring += otherPop + "\t"
			args.append(argstring.strip('\t'))

		commandstring = "./combine_cms_scores_multiplepops" 
		writeJobsToTaskArray_fromArglist(commandstring, args, prefix + "_multiple", UGERoptions=UGERoptions)

main()