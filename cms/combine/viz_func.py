## functions for viewing haplotypes original Shervin Tabrizi update Joe Vitti
## last updated 07.27.17 	vitti@broadinstitute.org

import matplotlib
import numpy as np
import sys
import os

def load_from_hap(inputfilename, filter_maf=None, corePos = None):
	#openfile = open(inputfilename)
	#header = openfile.readline()
	haplotypes = []
	physpositions = []
	h = open(inputfilename, 'r')
	index = 0
	coreindex = -1
	if h.readline().strip() == "##format=hapmap2transposed":
		transpose = True
	else:
		transpose = False
		h.seek(0)

	if filter_maf == None:
		filter_maf = -1.
	else:
		filter_maf = float(filter_maf)

	for line in h: 
		genotypes = line.strip().split()[1]
		int_genotypes = [int(x) for x in genotypes]
		n_a0 = float(int_genotypes.count(0))
		n_a1 = float(int_genotypes.count(1))
		n = float(len(int_genotypes))
		f_a0 = n_a0 / n
		f_a1 = n_a1 / n
		if f_a0 < filter_maf or f_a1 < filter_maf:
			pass
		else:
			haplotypes.append(''.join(genotypes))
			entries = line.split()
			label = entries[0]
			physpos = int(label.split('.')[1])
			physpositions.append(physpos)
			if physpos == corePos:
				coreindex = index
			index +=1
	h.close()

	if transpose:
		haplotypes2 = ['']*len(haplotypes[0])
		for i in haplotypes:
			for j in range(len(i)):
				haplotypes2[j] += i[j]
	haplotypes = haplotypes2
	#print(str(len(haplotypes)))
	return haplotypes, coreindex, physpositions
def pullRegion(inputfilename, startpos, endpos, maf=None, corePos = None, transpose = True, saveHapFile = None):
	coreindex = -1
	varids, genpositions, physpositions, all_genotypes = [], [], [], []
	if ".gz" in inputfilename:
		print("unzip!") #should implement gzip
		sys.exit(0)
		#infile = gzip.open(inputfilename, 'rb')
	infile = open(inputfilename, 'r')
	i = 0
	for line in infile:
		entries = line.split()
		chrom, varid, genpos, physpos, genotypes = entries[0], entries[1], entries[2], entries[3], entries[4:]
		physpos = int(physpos)
		if corePos is not None:
			if physpos == int(corePos):
				coreindex = i
		if int(physpos) >= startpos and physpos <= endpos:
			varids.append(varid)
			genpositions.append(float(genpos))
			physpositions.append(physpos)
			all_genotypes.append(genotypes)
			i+=1
		elif physpos > endpos:
			break
		
	infile.close()

	if maf is not None:
		filtered_haps = []
		filtered_pos = []
		filter_maf = float(maf)
		for ilist in range(len(all_genotypes)):
			genotypelist = all_genotypes[ilist]
			thisPhys = physpositions[ilist]
		#for genotypelist in all_genotypes:
			genotypes = [int(x) for x in genotypelist]
			n_a0 = float(genotypes.count(0))
			n_a1 = float(genotypes.count(1))
			n = float(len(genotypes))

			f_a0 = n_a0 / n
			f_a1 = n_a1 / n
			if f_a0 < filter_maf or f_a1 < filter_maf:
				pass
			else:
				filtered_haps.append(genotypelist)
				filtered_pos.append(thisPhys)
		if corePos is not None:
			coreindex = filtered_pos.index(corePos)
		all_genotypes = filtered_haps


	haplotypes = []
	for genotypelist in all_genotypes: 
		haplotype = ''.join(genotypelist)
		haplotypes.append(haplotype)

	if transpose:
		haplotypes2 = zip(*haplotypes)
		haplotypes = list(haplotypes2)

	if saveHapFile is not None:
		print('come back to this')
		pass

	return haplotypes, coreindex, physpositions
def difference(hap1,hap2):
	dif = 0
	for i in range(len(hap1)):
		if hap1[i] != hap2[i]: dif += 1
	return dif
def hapSort(haplotypes):
	"""
	sorts haplotypes to better visualize haplotype blocks
	"""

	haplotypes_sorted = []	
	
	print("Now getting differences")
	diffArray = np.zeros( (len(haplotypes), len(haplotypes)) )
	
	
	for i in range(len(haplotypes)):
		if i % 10 == 0: print("now haplotype " + str(i) + " out of " + str(len(haplotypes)))
		hapsToCompare = range(i,len(haplotypes))
		a = [difference(haplotypes[i],haplotypes[x]) for x in hapsToCompare]
		for j in range(len(hapsToCompare)):
			diffArray[i][hapsToCompare[j]] = a[j]
			diffArray[hapsToCompare[j]][i] = a[j]
	
	sums = []
	for row in diffArray:
		sums.append(sum(row))
	
	indices = []

	hapIndex = sums.index(min(sums))
	indices.append(hapIndex)
		
	for i in range(len( diffArray )): diffArray[i][hapIndex] = 1e6

	while len(indices) < len(haplotypes):
		diffList = list(diffArray[hapIndex])
		hapIndex = diffList.index(min(diffList))
		indices.append(hapIndex)
		for i in range(len( diffArray )): diffArray[i][hapIndex] = 1e6
		
	for index in indices:
		haplotypes_sorted.append( haplotypes[index] )
	
	assert len(haplotypes) == len(haplotypes_sorted)
	#print("length of orig haps: " +str(len(haplotypes)))
	#print("length of sorted haps: " + str(len(haplotypes_sorted)))

	return [haplotypes_sorted,indices]
def hapSort_coreallele(haplotypes, coreindex):
	"""sorts haplotypes to better visualize haplotype blocks
	THIS VERS: separates based on core allele specified as string"""

	if coreindex == -1:
		#print("did not find core snp in " + infilename)
		return hapSortSub(haplotypes)

	else:

		zerohaps, onehaps, twohaps = [], [], [] #partition based on allele at core site

		for i in range(len(haplotypes)):
			if haplotypes[i][coreindex] == '0':
				zerohaps.append(haplotypes[i])
				#zerocount +=1
			elif haplotypes[i][coreindex] == '1':
				onehaps.append(haplotypes[i])
				#onecount +=1
			elif haplotypes[i][coreindex] == '2':
				twohaps.append(haplotypes[i])
				#twocount +=1   

		zerosort = hapSortSub(zerohaps)
		onesort = hapSortSub(onehaps)
		#twosort = hapSortSub(twohaps)

		allhaps, allindices = [], []
		allhaps.extend(zerosort[0])
		allhaps.extend(onesort[0])
		#allhaps.extend(twosort[0])
		allindices.extend(zerosort[1])
		allindices.extend(onesort[1])
		#allindices.extend(twosort[1])
		return [allhaps, allindices, coreindex]
def hapSortSub(haplotypes):
	"""helper to hapSort_coreallele"""
	haplotypes_sorted = []	
	
	print("Now getting differences")
	diffArray = np.zeros( (len(haplotypes), len(haplotypes)) )
	
	
	for i in range(len(haplotypes)):
		if i % 10 == 0: print("now haplotype " + str(i) + " out of " + str(len(haplotypes)))
		hapsToCompare = range(i,len(haplotypes))
		a = [difference(haplotypes[i],haplotypes[x]) for x in hapsToCompare]
		for j in range(len(hapsToCompare)):
			diffArray[i][hapsToCompare[j]] = a[j]
			diffArray[hapsToCompare[j]][i] = a[j]
	
	sums = []
	for row in diffArray:
		sums.append(sum(row))
	
	indices = []

	hapIndex = sums.index(min(sums))
	indices.append(hapIndex)
		
	for i in range(len( diffArray )): diffArray[i][hapIndex] = 1e6

	while len(indices) < len(haplotypes):
		diffList = list(diffArray[hapIndex])
		hapIndex = diffList.index(min(diffList))
		indices.append(hapIndex)
		for i in range(len( diffArray )): diffArray[i][hapIndex] = 1e6
		
	for index in indices:
		haplotypes_sorted.append( haplotypes[index] )
	
	assert len(haplotypes) == len(haplotypes_sorted)
	#print("length of orig haps: " +str(len(haplotypes)))
	#print("length of sorted haps: " + str(len(haplotypes_sorted)))

	return [haplotypes_sorted,indices]
def hapViz(ax, haplotypes,out,coreindex = None,markers=None):
	"""
	main function for viewing haplotypes
	takes list with one haplotype per entry (as string of 1's and 0's)
	makes haplotype figure as output
	"""	

	# colors for haplotypes
	col0 = '#ABDDA4' #green
	col1 = '#2B83BA' #blue
	col2 = '#d7191c' #red
	
	my_cmap = matplotlib.colors.ListedColormap([col0, col1])#[col0, col1, col2])
	my_cmap.set_bad(color='w', alpha=0)

	num_haps = len(haplotypes)
	num_snps = len(haplotypes[0]) 
	
	hap_spacer = 1
	if num_snps > 500:
		hap_height = 10 
	else:
		hap_height = 1
	snp_width = 1
	xCoord = [i*snp_width for i in range(num_snps)]
	yCoord = [j*hap_height for j in range(num_haps)] 

	ax.axis([0, num_snps*snp_width, 0, num_haps * hap_height + (num_haps - 1 )* hap_spacer])

	spacer = [float('nan') for x in range(num_snps)]
	spacer_array = np.array(spacer)

	data1 = []#np.array([])
	for haplotype in haplotypes:
		genotype_chars = list(haplotype)
		genotype_ints = [int(x) for x in genotype_chars]
		genotype_array = np.array(genotype_ints)
		for iPixel in range(hap_height):
			data1.append(genotype_array)
		for jPixel in range(hap_spacer):
			data1.append(spacer_array)


	ax.imshow(data1, interpolation='none', cmap=my_cmap)#,  zorder=0) #extent=[0, num_snps, 0, num_haps],

	for tick in ax.xaxis.get_major_ticks():
		tick.tick2On = False
		tick.label1On = False
		tick.tick1On = False
	for tick in ax.yaxis.get_major_ticks():
		tick.tick2On = False
		tick.tick1On = False
		tick.label1On = False
	
	for border in ['right', 'left', 'top', 'bottom']:
		ax.spines[border].set_color('none')

	return ax
def readAnnotations(annotationfilename):
	if not os.path.isfile(annotationfilename):
		print("missing file: " + annotationfilename)
		return [], []
	else:
		positions, annotations = [], []
		openfile = open(annotationfilename, 'r')
		for line in openfile:
			entries = line.split()
			positions.append(entries[0])
			annotations.append(' '.join(entries[1:]))
		openfile.close()
	return positions, annotations
def find_snp_index(infilename, snppos, startpos):
	h = open(infilename, 'r')##
	coreindex = -1
	indexcounter = 0
	for line in h: 
		pos = line.strip().split()[3]
		if int(pos) == int(snppos):
			#print("found the core snp")
			coreindex = indexcounter
			break
		if int(pos) < startpos:
			pass
		else:
			indexcounter +=1
	#print(str(coreindex))
	return indexcounter
