## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
## call program as: python get_neutral_targetstats_frombootstrap_thinned_bysnp.py {nbootstrapreps}
## last updated: 06.04.16 	vitti@broadinstitute.org
# JV: currently configured for 1KG P3 Afr pops. (vers 12.08.15: _gapped (nb!))

dateprefix, fileprefix = "120915_3pm_chr21-22",'results_inclusive'
fileloc = "/idi/sabeti-scratch/jvitti/dem_model/target_vals/neutral/bootstrap/thinned/"

###############################
## TOGGLE THIS RUN OF SCRIPT ##
###############################

getFreqs, getLD, getFst = True, True, True
pops = ['LWK', 'GWD', 'MSL', 'YRI']
chroms = [21, 22]#range(1,23)
mafcutoffdprime = .2

if True:
	print "should store this in one spot for whole pipeline and import it here"
	num_indivs = {'GWD': 113, 'ACB': 96, 'BEB': 86, 'PEL': 85, 'LWK': 99, 'MSL': 85, 'GBR': 91, 'IBS': 107, 'ASW': 61, 'TSI': 107, 'KHV': 99, 'CEU': 99, 'YRI': 108, 'CHB': 103, 'STU': 102, 'CHS': 105, 'ESN': 99, 'FIN': 99, 'GIH': 103, 'PJL': 96, 'MXL': 64, 'PUR': 104, 'ITU': 102, 'CDX': 93, 'JPT': 104, 'CLM': 94}
	nlimiting = 2*min([num_indivs[x] for x in pops])

from bootstrap_func import *
import numpy as np
import random
import sys

############
### MAIN ###
############

def main():

	if len(sys.argv) != 2:
		print "call program as: python get_neutral_targetstats_frombootstrap_thinned_bysnp.py {nbootstrapreps}"

	nbootstraprep = int(sys.argv[1])

	writefilename = fileloc + fileprefix + "_targetstats_thinned_bootstrap_bysnp_" +  str(nbootstraprep) + "reps_" + dateprefix  + ".txt"
	print "writing to " + writefilename
	writefile = open(writefilename, 'w')

	################ 
	## FREQ STATS ##
	################
	nhist = 6
	if getFreqs == True:
		for pop in pops:
			writefile.write(pop + '\n')
			allRegionDER, allRegionANC, allRegionPI, allseqlens = [], [], [], []
			nsnps, totalregions, totallen = 0, 0, 0
			for chrom in chroms:
				freqsfilename = fileloc + pop + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_freqs_by_snp.txt"
				if chrom == 1:
					print "reading from eg: " + freqsfilename
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

			print "\tlogged frequency values for " + str(nsnps) + " SNPS across " + str(totalregions) + " regions covering " + str(totallen) + " bp."
			
			###########
			#### PI ###
			###########
			#print allRegionPI
			#print allseqlens
			pi_mean = estimatePi(allRegionPI, allseqlens)
			print str(pi_mean)
			writefile.write(str(pi_mean)+'\t')

			estimates = []
			#print str(allRegionPI)
			for j in range(nbootstraprep):
				rep_pis, rep_seqlens = [], []
				for k in range(totalregions):
					index = random.randint(0, totalregions-1)
					rep_pis.append(allRegionPI[index])
					rep_seqlens.append(allseqlens[index])
				rep_pi_mean = estimatePi(rep_pis, rep_seqlens)
				estimates.append(rep_pi_mean)
			#print str(estimates)
			pi_se = np.std(estimates)
			print str(pi_se)
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
				#flatten anc and der
				flatanc = flattenList(allRegionANC)
				flatder = flattenList(allRegionDER)
				for w in range(nsnps):
					index = random.randint(0, nsnps-1)
					rep_all_nderiv.append(flatder[index])
					rep_all_nanc.append(flatanc[index])
				#print str(rep_all_nderiv)
				#print str(rep_all_nanc)
				## CALCULATE SPECTRA FOR REGIONS	
				repmafhist, repanchist = estimateFreqSpectrum(rep_all_nderiv, rep_all_nanc, nhist)
				npoly = sum(repmafhist)
				repsfs = [float(x)/npoly for x in repmafhist]
				for ibin in range(nhist):
					estimates_sfs[ibin].append(repsfs[ibin])
				#I DONT THINK I HAVE ENOUGH REGIONS/VARIATIONS FOR THIS TO BOOTSTRAP PROPERLY. Break up 100kb regions into smaller onces and sample a greater number of them? 
				#getting division by zero error.	
				#print repmafhist
				#repderhist = [mafhist[i] - anchist[i] for i in range(nhist)]
				#repder = [repderhist[i]/float(repmafhist[i]) for i in range(nhist)]
				#repanc = [1-x for x in repder]
				repanc = [repanchist[i]/float(repmafhist[i]) for i in range(nhist)]
				for ibin in range(nhist):
					estimates_anc[ibin].append(repanc[ibin])



			sfs_se = [np.std(x) for x in estimates_sfs]
			anc_se = [np.std(x) for x in estimates_anc]
			
			print sfs_mean
			print sfs_se
			print anc_mean
			print anc_se
			
			writefile.write(str(sfs_mean) + '\n')
			writefile.write(str(sfs_se) + '\n')
			writefile.write(str(anc_mean) + '\n')
			writefile.write(str(anc_se) + '\n')

	#########
	## LD ##
	#########
	nphysdisthist, ngendisthist = 14, 17
	if getLD == True:
		for pop in pops:
			allRegionDists, allRegionGendists, allRegionr2, allRegionDprime = [], [], [], []
			N_r2regs, N_dprimeregs = 0, 0
			N_r2snps, N_dprimesnps = 0, 0
			for chrom in chroms:
				ldfilename = fileloc + pop + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_ld_by_snp_pair_gapped.txt"#_vsgw.txt"
				print "reading from: " + ldfilename
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

			print "TOTAL: "
			print "\tlogged r2 values for " + str(N_r2snps) + " SNP pairs across " + str(N_r2regs) + " regions."
			print "\tlogged D' values for " + str(N_dprimesnps) + " SNP pairs across " + str(N_dprimeregs) + " regions."

			###################################
			### r2: MEAN ACROSS ALL REGIONS ###
			###################################
			r2sums, physDistHist = estimater2decay(allRegionr2, allRegionDists)
			r2dist = [r2sums[u]/physDistHist[u] for u in range(len(r2sums))]
			#print str(r2sums)
			#print str(physDistHist)
			print str(r2dist)
			writefile.write(str(r2dist) + "\n")

			############################################
			### r2: STDERR ACROSS BOOTSTRAP ESTIMATES ##
			############################################
			estimates_r2 = [[] for i in range(nphysdisthist)]
			while len(estimates_r2[0]) < nbootstraprep:
				rep_all_r2, rep_all_physdist = [], []
				#flatten
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
			#print str(r2estimate)
			#print str(estimates_r2)

			r2_se = [np.std(x) for x in estimates_r2]
			print str(r2_se)
			writefile.write(str(r2_se) + "\n")

			####################################
			### D': MEAN ACROSS ALL REGIONS ###
			####################################
			compLDhist, genDistHist = estimatedprimedecay(allRegionDprime, allRegionGendists)
			dprimedist = [float(compLDhist[x])/float(genDistHist[x]) for x in range(len(compLDhist))]
			#print str(compLDhist)
			#print str(genDistHist)
			print str(dprimedist)
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
				#if 0 not in repgenDistHist:
				for ibin in range(len(repgenDistHist)):
					if repgenDistHist[ibin] == 0:
						repgenDistHist[ibin] = 1
				dprimeestimate = [float(repcompLDhist[x])/float(repgenDistHist[x]) for x in range(ngendisthist)]
				for ibin in range(ngendisthist):
					estimates_dprime[ibin].append(dprimeestimate[ibin])
				#print str(dprimeestimate)
			dprime_se = [np.std(x) for x in estimates_dprime]
			print str(dprime_se)
			writefile.write(str(dprime_se) + "\n")


	##########
	## Fst ##
	##########
	if getFst == True:
		for i in range(len(pops)):
			for j in range(i, len(pops)):
				if i != j:
					allRegionValues = []
					for chrom in chroms:
						pop0, pop1 = pops[i], pops[j]
						fstfilename = fileloc + pop0 + "_" + pop1 + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_fst_by_snp.txt"
						if not checkFileExists(fstfilename):
							fstfilename = fileloc + pop1 + "_" + pop0 + "_" +str(nlimiting) + "_"+ fileprefix + "_chr" + str(chrom) + "_fst_by_snp.txt"
						if chrom == 1:
							print "reading from e.g. " + fstfilename
						if checkFileExists(fstfilename):
							allfst, nregions = readFstFile(fstfilename)
							allRegionValues.extend(allfst)
					target_mean, target_se = estimateFstByBootstrap_bysnp(allRegionValues, nrep = nbootstraprep)
					writeline = pop0 + "\t" + pop1 + "\t" + str(target_mean) + "\t" + str(target_se) + '\n'
					writefile.write(writeline)
					print writeline


	writefile.close()
	print "wrote to: " + writefilename

main()