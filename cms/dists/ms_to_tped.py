## {{DEM MODEL --> COMPONENT SCORE PROBABILITY DISTS}} -- script ib.
## last updated 05.31.16; update to fixStrictlyIncreasing() needs to be validated, should fix simulations that get 'hung up' (presumably due to scenario where three adjacent SNPs exist?)

import sys

def rewriteGenotypes(genotypes, pop_index):
	#in the tped, 1 should be anc and 0 should be derived
	#in ms output, ancestral is 0 and the derived is 1
	#in cosi full output, ancestral is 2 and derived is 1
	rewrite = []
	for item in genotypes:
		if item == '0':
			rewrite.append('1')
		if item == '1':
			rewrite.append('0')
	assert len(genotypes) == len(rewrite)
	return rewrite
def checkStrictlyIncreasing(numberlist):
	"""ensures no dups and monotonic increase, for writing locs to tped from ms"""
	return all(x < y for x, y in zip(numberlist, numberlist[1:]))
def fixStrictlyIncreasing(numberlist, causal_index, check = True):
	"""cheats out any duplicate positions that arise due to rounding"""
	if numberlist[causal_index] != 500000 and check == True:
		print "ERROR! got the wrong causal index"
	while not checkStrictlyIncreasing(numberlist):
		#spread out in both directions from causal variant
		for ivar in range(causal_index, len(numberlist)-1):
			interval = 1
			while numberlist[ivar] == numberlist[ivar+1]: #former IF
				if (numberlist[ivar+1] + 1) not in numberlist:
					numberlist[ivar+1] +=interval
				elif (numberlist[ivar] - 1) not in numberlist:
					numberlist[ivar] -=interval
				interval +=1

		for ivar in range(causal_index, 0, -1):
			interval = 1
			while numberlist[ivar] == numberlist[ivar-1]: #former IF
				if (numberlist[ivar] + 1) not in numberlist:
					numberlist[ivar] += interval
				elif (numberlist[ivar-1] - 1) not in numberlist:
					numberlist[ivar-1] -= interval
				interval +=1
	return numberlist
def MS_TO_TPED(sourcefilename, writefileprefix, regionlen, checkCausalVar = True, numPops = 4, numIndivsPerSamp=172, writeJson = True, causalPos = 0.5):
	openfile = open(sourcefilename, 'r')
	print "\tReading cosi output from "+ sourcefilename
	for i in range(3): #strip cosi header
		openfile.readline()
	seedline = openfile.readline()
	assert "//" in seedline
	genlenline = openfile.readline()
	segsitesline = openfile.readline()
	segsites = int(segsitesline.split()[-1])
	posline = openfile.readline()
	genposline = openfile.readline()

	genlen = float(genlenline.split()[-1])
	genpos_frac = [float(x) for x in genposline.split()[1:]]
	physpos_frac = [float(x) for x in posline.split()[1:]]
	genpos = [x * genlen for x in genpos_frac]
	physpos = [(x * float(regionlen)) for x in physpos_frac]
	physpos = [int(y) for y in physpos]
	assert segsites == len(physpos)

	#print str(bool(checkCausalVar))
	if "T" in str(checkCausalVar) or 't' in str(checkCausalVar): #is this the best way to handle this?
	#if checkCausalVar == "True":
	#if bool(checkCausalVar) == True:
		assert causalPos in physpos_frac
		causal_index = physpos_frac.index(causalPos)
		#print "found causal index: " + str(causal_index)
		#print str(physpos[causal_index])
		#print str(genpos[causal_index])

		#CHECK NO DUPLICATES
		#check to make sure phys positions are strictly increasing
	else:
		causal_index = 0
	if checkStrictlyIncreasing(physpos) == False:
		#print "error with phys positions"
		physpos = fixStrictlyIncreasing(physpos,causal_index, check = checkCausalVar)

	#median = len(physpos)/2
	#physpos = fixStrictlyIncreasing(physpos,median)
	assert (checkStrictlyIncreasing(physpos))
	assert (type(x) == int for x in physpos)


	haplotypes = []
	for j in range(numPops):
		for k in range(numIndivsPerSamp):
			haplotypeline = openfile.readline()
			haplotypes.append(haplotypeline.strip('\n'))
	openfile.readline() #trailer

	#transpose haplotypes #(should use izip?)
	genotypes = zip(*haplotypes)
	if len(genotypes) != segsites:
		print str(segsites)
		print str(len(genotypes))
		print str(genotypes)
		print str(i)
	assert len(genotypes) == segsites
	assert len(genotypes[0]) == (numPops * numIndivsPerSamp)

	#write a tped for each pop in each rep!! 
	for j in range(numPops):

		fileprefix = writefileprefix + "_pop" + str(j+1) 
		tpedfilename = fileprefix + ".tped"
		jsonfilename = writefileprefix + ".metadata.json"
		tpedfile = open(tpedfilename, 'w') 
		for s in range(segsites):	
			rewrittenGenotypes = rewriteGenotypes(genotypes[s], j)

			towritestring = "1 " + str(s) + " " + str(genpos[s]) + " "+ str(physpos[s]) + " "
			#str(i) + " " + str(s) + " " + str(genpos[s]) + " "+ str(physpos[s]) + " "
			for w in range(len(rewrittenGenotypes)):
				#
				start = j * numIndivsPerSamp
				end = start + numIndivsPerSamp
				#print str(start) + " - " + str(end)
				if w >= start and w < end:
					towritestring += str(rewrittenGenotypes[w]) + " "
				#towritestring.strip(' ')
			tpedfile.write(towritestring[:-1] + "\n")
			#tpedfile.write('\n')
		tpedfile.close()
		if writeJson == True:
			if j == 1:
				jsonfile = open(jsonfilename, 'w')
				jsonfile.write("{\n\t\"cosi_ms_file\": \"" + sourcefilename + "\",\n}\n")
				#jsonfile.write("{\n\t\"causal_variant\": \"" + str(physpos[causal_index]) + "\",\n}\n")
				
				jsonfile.close()
		#tpedfile.close()
	openfile.close()
	print "\twrote output to " + writefileprefix + "_{rep}_{pop}.tped"

###########
## MAIN ###
###########

def main():
	sourcefilename, writefileprefix, regionlen, checkCausal = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	MS_TO_TPED(sourcefilename, writefileprefix, regionlen, checkCausal)

main()