## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
## quick script; rewrites coordinates so they are easy to parse in basic C 
## last updated: 11.17.15: add in basic filters## last updated: 06.04.16 	vitti@broadinstitute.org

filterstring = "test_filter"
regionfolder = '/idi/sabeti-scratch/jvitti/dem_model/target_vals/neutral/regions/'
chroms = range(1,23)
prefixes = ['results_inclusive']
savedir = '/home/unix/vitti/public_html/dem_model/target_vals/'
publicdir = "http://broadinstitute.org/~vitti/dem_model/target_vals/"

import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as plt
import pylab
import numpy as np

#############
### MAIN ####
#############

for prefix in prefixes:
	regionlens = []
	filteredregionlens = []
	for chrom in chroms:
		sourcefilename = regionfolder + prefix + "_chr" + str(chrom) + ".tsv"
		targetfilename = regionfolder + prefix + "_" + filterstring + "_chr" + str(chrom) + ".txt"
		print "reading from: " + sourcefilename
		readfile = open(sourcefilename, 'r')
		print "writing to: " + targetfilename
		writefile = open(targetfilename, 'w')

		readfile.readline() # stripheader
		

		for line in readfile:
			entries = line.split()
			chromstring, start, end, regionid, st_d2g_cM, st_d2g_bp, st_recomb, st_diversity, st_bkgd, bp_simpleRepeat, bp_rmskRM327, bp_phastConsElements46wayPlacMammal, bp_allSoftRegions, seGenMap, perc_simpleRepeat, perc_rmskRM327, perc_phastConsElements46wayPlacMammal, perc_allSoftRegions = entries
			printchromosome = chromstring[3:]
			start = int(start)
			end = int(end)
			regionlen = (end - start)

			#####################
			### DEFINE FILTER ##
			##vvvvvvvvvvvvvvvv##
			passFilter = True

			if regionlen < 10:
				#print "discarding region with length: " + str(regionlen)
				passFilter = False

			#if float(st_d2g_cM) < .1 or int(st_d2g_bp) < 100000:
			#	passFilter = False

			#if float(perc_simpleRepeat) >.99 or float(perc_rmskRM327) > .99 or float(perc_phastConsElements46wayPlacMammal) > .99:
			#	passFilter = False

			##^^^^^^^^^^^^^^^^^##
			### DEFINE FILTER ##
			#####################

			if passFilter == True:
				regionlens.append(regionlen)
				writestring = str(printchromosome) + "\t" + str(start) + "\t" + str(end) + '\n'
				writefile.write(writestring)
			else:
				filteredregionlens.append(regionlen)
		readfile.close()	
		writefile.close()
	print "num regions: " + str(len(regionlens))
	print "total length: " + str(sum(regionlens))
	print "filtered out " + str(len(filteredregionlens)) + " regions covering " + str(sum(filteredregionlens)) + " bp."

	########################
	### PLOT DIST OF LENS ##
	########################
	weights_regionlens = np.ones_like(regionlens)/len(regionlens)

	fig = plt.figure()
	fig.add_subplot(111)
	nbins = 100
	n, bins, patches = pylab.hist(regionlens,  bins=nbins, histtype = 'step') #weights=weights_regionlens,
	plt.xlabel("region len (bp)")
	plt.ylabel("nregions")
	plt.title("distribution of region lengths in dataset: " )
	savetitle = savedir + prefix + "_regionlens.png"
	plt.savefig(savetitle)   #CHANGEABLE
	plt.close()
	print "plotted to: " + publicdir + prefix + "_regionlens.png"