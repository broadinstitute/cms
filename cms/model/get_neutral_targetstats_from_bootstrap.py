## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
## call program as: python get_neutral_targetstats_frombootstrap.py {nbootstrapreps}
## last updated: 06.16.16 	vitti@broadinstitute.org

def retrievePerSnpBootstrapFilename(chrom, statType = "freqstats", datestring = "061616", altPop = ""):
	"""user can specify the output filename and location for bootstrap programs in C from the command line.
	must designate those same filenames in this function."""
	fileloc = "/idi/sabeti-scratch/jvitti/cms_venv/dry_run/popstats/"
	if statType != "fst":
		filename = fileloc + "IRN_" + statType + "_chr" + str(chrom) + "_" + datestring
	else:
		filename = fileloc + "IRN_" + altPop  + "_fst_chr" + str(chrom) + "_" + datestring
	if not os.path.isfile(filename):
		print "error: " + filename
	return filename

###############################
## TOGGLE THIS RUN OF SCRIPT ##
###############################
writefileprefix, datestring = "bootstrap_stats_gw", "061616"
getFreqs, getLD, getFst = True, True, True
fstall, altpops = False, ['BEB', 'CHB', 'YRI', 'CEU'] #if fstall == True, takes all pairs from pops
pops = ['IRN']
chroms = range(1,23)

###########################
## DEFINE BINS FOR DISTS ##
###########################
mafcutoffdprime = .2
nhist = 6
nphysdisthist, ngendisthist = 14, 17

from bootstrap_func import *
import numpy as np
import random
import sys

############
### MAIN ###
############

def main():

	if len(sys.argv) != 2:
		print "call program as: python get_neutral_targetstats_frombootstrap.py {nbootstrapreps}"
		sys.exit(1)

	nbootstraprep = int(sys.argv[1])
	writefilename = writefileprefix + "_" +  str(nbootstraprep) + "reps_" + datestring  + ".txt"
	print "writing to " + writefilename
	writefile = open(writefilename, 'w')

	################ 
	## FREQ STATS ##
	################
	if getFreqs == True:
		for pop in pops:
			print "FREQS, " + pop  + ":"
			writefile.write(pop + '\n')
			allRegionDER, allRegionANC, allRegionPI, allseqlens = [], [], [], []
			nsnps, totalregions, totallen = 0, 0, 0
			for chrom in chroms:
				freqsfilename = retrievePerSnpBootstrapFilename(chrom, "freqstats", datestring)
				if checkFileExists(freqsfilename):
					allpi, allnderiv, allnanc, nregions, seqlens = readFreqsFile(freqsfilename)
					allRegionPI.extend(allpi)
					allRegionDER.extend(allnderiv)
					allRegionANC.extend(allnanc)
					allseqlens.extend(seqlens)
					totalregions += nregions
					totallen += sum(seqlens)
					for i in range(len(allpi)):
						nsnps += len(allpi[i])
			print "TOTAL: logged frequency values for " + str(nsnps) + " SNPS across " + str(totalregions) + " regions covering " + str(totallen) + " bp.\n"
			
			###########
			#### PI ###
			###########
			pi_mean = estimatePi(allRegionPI, allseqlens)
			writefile.write(str(pi_mean)+'\t')
			estimates = []
			for j in range(nbootstraprep):
				rep_pis, rep_seqlens = [], []
				for k in range(totalregions):
					index = random.randint(0, totalregions-1)
					rep_pis.append(allRegionPI[index])
					rep_seqlens.append(allseqlens[index])
				rep_pi_mean = estimatePi(rep_pis, rep_seqlens)
				estimates.append(rep_pi_mean)
			pi_se = np.std(estimates)
			writefile.write(str(pi_se) + '\n')

			#########################################
			### SFS, ANC: MEAN ACROSS ALL REGIONS ###
			#########################################
			mafhist, anchist = estimateFreqSpectrum(allRegionDER, allRegionANC, nhist)
			npoly = sum(mafhist)
			sfs_mean = [float(x)/npoly for x in mafhist]
			anc_mean = [anchist[i]/float(mafhist[i]) for i in range(len(mafhist))]

			#########################################
			### STDERR ACROSS BOOTSTRAP ESTIMATES ###
			#########################################		
			estimates_sfs, estimates_anc = [[] for i in range(nhist)], [[] for i in range(nhist)]
			for j in range(nbootstraprep):
				rep_all_nderiv, rep_all_nanc = [], []
				flatanc = flattenList(allRegionANC)
				flatder = flattenList(allRegionDER)
				for w in range(nsnps):
					index = random.randint(0, nsnps-1)
					rep_all_nderiv.append(flatder[index])
					rep_all_nanc.append(flatanc[index])
				repmafhist, repanchist = estimateFreqSpectrum(rep_all_nderiv, rep_all_nanc, nhist)
				npoly = sum(repmafhist)
				repsfs = [float(x)/npoly for x in repmafhist]
				for ibin in range(nhist):
					estimates_sfs[ibin].append(repsfs[ibin])
				repanc = [repanchist[i]/float(repmafhist[i]) for i in range(nhist)]
				for ibin in range(nhist):
					estimates_anc[ibin].append(repanc[ibin])

			sfs_se = [np.std(x) for x in estimates_sfs]
			anc_se = [np.std(x) for x in estimates_anc]
			writefile.write(str(sfs_mean) + '\n')
			writefile.write(str(sfs_se) + '\n')
			writefile.write(str(anc_mean) + '\n')
			writefile.write(str(anc_se) + '\n')

	#########
	## LD ##
	#########
	if getLD == True:
		for pop in pops:
			print "LD, " + pop + ":"
			allRegionDists, allRegionGendists, allRegionr2, allRegionDprime = [], [], [], []
			N_r2regs, N_dprimeregs = 0, 0
			N_r2snps, N_dprimesnps = 0, 0
			for chrom in chroms:
				ldfilename = retrievePerSnpBootstrapFilename(chrom, statType = "ldstats", datestring=datestring)
				if checkFileExists(ldfilename):

					alldists, allr2, allgendists, alldprime, nr2regions, ndprimeregions = readLDFile(ldfilename, dprimecutoff = mafcutoffdprime)

					allRegionDists.extend(alldists)
					allRegionr2.extend(allr2)
					N_r2regs += nr2regions
					N_r2snps += sum([len(x) for x in allr2])

					allRegionGendists.extend(allgendists)
					allRegionDprime.extend(alldprime)
					N_dprimeregs += ndprimeregions
					N_dprimesnps += sum([len(x) for x in alldprime])

			print "TOTAL: logged r2 values for " + str(N_r2snps) + " SNP pairs.\n\tlogged D' values for " + str(N_dprimesnps) + " SNP pairs.\n"

			###################################
			### r2: MEAN ACROSS ALL REGIONS ###
			###################################
			r2sums, physDistHist = estimater2decay(allRegionr2, allRegionDists)
			r2dist = [r2sums[u]/physDistHist[u] for u in range(len(r2sums))]
			writefile.write(str(r2dist) + "\n")

			############################################
			### r2: STDERR ACROSS BOOTSTRAP ESTIMATES ##
			############################################
			estimates_r2 = [[] for i in range(nphysdisthist)]
			while len(estimates_r2[0]) < nbootstraprep:
				rep_all_r2, rep_all_physdist = [], []
				flatr2 = flattenList(allRegionr2)
				flatregions = flattenList(allRegionDists)
				nsnppairs = len(flatr2)
				for w in range(nsnppairs):
					index_r2 = random.randint(0, nsnppairs-1)
					rep_all_r2.append(flatr2[index_r2])
					rep_all_physdist.append(flatregions[index_r2])

				#add pseudocount for empty bins
				repr2sum, repphysdisthist = estimater2decay(rep_all_r2, rep_all_physdist)
				for ibin in range(len(repphysdisthist)):
					if repphysdisthist[ibin] == 0:
						repphysdisthist[ibin] = 1
				r2estimate =[repr2sum[u]/repphysdisthist[u] for u in range(len(repr2sum))]
				for ibin in range(nphysdisthist):
					estimates_r2[ibin].append(r2estimate[ibin])

			r2_se = [np.std(x) for x in estimates_r2]
			writefile.write(str(r2_se) + "\n")

			####################################
			### D': MEAN ACROSS ALL REGIONS ###
			####################################
			compLDhist, genDistHist = estimatedprimedecay(allRegionDprime, allRegionGendists)
			#add pseudocounts
			for ibin in range(len(genDistHist)):
				if genDistHist[ibin] == 0:
					genDistHist[ibin]+=1

			dprimedist = [float(compLDhist[x])/float(genDistHist[x]) for x in range(len(compLDhist))]
			writefile.write(str(dprimedist) + "\n")

			############################################
			### D': STDERR ACROSS BOOTSTRAP ESTIMATES ##
			############################################
			estimates_dprime = [[] for i in range(ngendisthist)]
			while len(estimates_dprime[0]) < nbootstraprep:
				rep_all_dprime, rep_all_gendist = [], []
				flatdprime = flattenList(allRegionDprime)
				flatgendist = flattenList(allRegionGendists)
				nsnppairs = len(flatdprime)

				for w in range(nsnppairs):
					index_dprime = random.randint(0, nsnppairs-1)
					rep_all_dprime.append(flatdprime[index_dprime])
					rep_all_gendist.append(flatgendist[index_dprime])

				repcompLDhist, repgenDistHist = estimatedprimedecay(rep_all_dprime, rep_all_gendist)
				for ibin in range(len(repgenDistHist)):
					if repgenDistHist[ibin] == 0:
						repgenDistHist[ibin] = 1
				dprimeestimate = [float(repcompLDhist[x])/float(repgenDistHist[x]) for x in range(ngendisthist)]
				for ibin in range(ngendisthist):
					estimates_dprime[ibin].append(dprimeestimate[ibin])
			dprime_se = [np.std(x) for x in estimates_dprime]
			writefile.write(str(dprime_se) + "\n")


	##########
	## Fst ##
	##########
	if getFst == True:
		print "FST:"
		for i in range(len(pops)):
			comp_pops = []
			if fstall == True:
				for j in range(i, len(pops)):
					if i != j:
						comp_pops.append(pops[j])
			else:
				comp_pops = altpops[:]

			for comppop in comp_pops:
				allRegionValues = []
				for chrom in chroms:
				
					fstfilename = retrievePerSnpBootstrapFilename(chrom, "fst", datestring, altPop = comppop)
					#fstfilename = fileloc + pop0 + "_" + pop1 + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_fst_by_snp.txt"
					#if not checkFileExists(fstfilename):
					#	fstfilename = fileloc + pop1 + "_" + pop0 + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_fst_by_snp.txt"
					if checkFileExists(fstfilename):
						allfst, nregions = readFstFile(fstfilename)
						allRegionValues.extend(allfst)
				target_mean, target_se = estimateFstByBootstrap_bysnp(allRegionValues, nrep = nbootstraprep)
				pop0, pop1 = pops[i], comppop
				writeline = pop0 + "\t" + pop1 + "\t" + str(target_mean) + "\t" + str(target_se) + '\n'
				writefile.write(writeline)
				print "TOTAL: logged Fst values for " + str(len(allRegionValues)) + " SNPs for " + pop0 + " and " + pop1 + ".\n"
				
	writefile.close()
	print "wrote to: " + writefilename

main()