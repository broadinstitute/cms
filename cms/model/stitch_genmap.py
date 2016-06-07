## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS: FILTER PUTATIVE NEUTRAL REGIONS}
## script to piece together genetic map from putative neutral regions to create a file that cosi can use. ## per IS suggestion, cuts regions before any peripheral hotspots so that rate is +- continuous
## JV 1KG P3 analysis builds on putative neutrality filters from Neutral Regions Explorer (Arbiza et al., 2012: http://nre.cb.bscb.cornell.edu/nre/)
## JV: still a lot needs be done to make this general. original source /idi/sabeti-scratch/jvitti/dem_model/target_vals/neutral/ 11.11.15
## last updated: 06.04.16 	vitti@broadinstitute.org

writefilename = "NRE_111115_10kb_plus_filtered_recom_regions_concat.txt"
genmapfile = "NRE_111115_10kb_plus_filtered_recom_regions.txt"
bufferBP = 1#500 #for stitching together 
hotspotcutoff = 100 #cM/Mb

import random

###################
## HELPER METHOD ##
###################

def scaleAndCut(positions, rates, hotspotcutoff = 5, scaleRate = 100000000):
	"""shifts region coordinates to zero, truncates if rate at edge > hotspotcutoff.
	scalerate converts HapMap cM/Mb rates to cosi rates."""
	assert len(positions) == len(rates)
	startindex, endindex = 0, len(positions)-1
	#check for hotspots
	for i in range(len(positions)):
		if rates[i] < hotspotcutoff:
			startindex = i
			break
	for j in range(len(positions)-1, -1, -1):
		if rates[j] < hotspotcutoff:
			endindex = j
			break
	toTakePos = positions[startindex:endindex]
	toTakeRates = rates[startindex:endindex]
	print str(toTakePos)
	startPos = toTakePos[0]
	toTakePos = [x - startPos for x in toTakePos]
	toTakeRates = [x / scaleRate for x in toTakeRates]
	return toTakePos, toTakeRates

##########
## MAIN ##
##########

def main():
	print "loading from... " + genmapfile
	openfile = open(genmapfile, 'r')
	positions_by_region, rates_by_region = [], []
	positions, rates = [], []
	debug = False
	for line in openfile:
		#finished a region
		if line == "\n" and len(rates) != 0:
			if len(rates) > 1:
				modifiedPositions, modifiedRates = scaleAndCut(positions, rates, hotspotcutoff = hotspotcutoff)
				if debug == True:
					print str(positions[-1])
					print str(rates[-1])
					print str(modifiedPositions[-1])
					print str(modifiedRates[-1])
					print str(len(positions)) + " vs " + str(len(modifiedPositions))
					print "*********************"
				positions_by_region.append(modifiedPositions)
				rates_by_region.append(modifiedRates)
				positions, rates = [], []
		#same region
		elif line != "\n":
			entries = line.split('\t')
			if len(entries) ==  3:
				pos, rate = int(entries[0]), float(entries[1])
				positions.append(pos)
				rates.append(rate)

	openfile.close()

	print "logged " + str(len(positions_by_region)) + " regions..."
	#print "randomizing order.."
	#choose random order to merge regions
	regions = range(len(positions_by_region))
	#print str(regions)
	random.shuffle(regions)
	#print str(regions)
	print "writing to file: " + writefilename
	writefile = open(writefilename, 'w')
	lastpos = 0
	for iregion in regions:
		thesePos, theseRates = positions_by_region[iregion], rates_by_region[iregion]
		for i in range(len(thesePos)):
			pos = thesePos[i]
			if pos == 0:
				pos += bufferBP
			pos += lastpos
			writeline = str(pos) + "\t" + str(theseRates[i]) + "\n"
			writefile.write(writeline)
			lastpos = pos
	writefile.close()

main()