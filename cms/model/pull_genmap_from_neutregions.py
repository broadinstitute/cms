## {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS: FILTER PUTATIVE NEUTRAL REGIONS}
## JV 1KG P3 analysis builds on putative neutrality filters from Neutral Regions Explorer (Arbiza et al., 2012: http://nre.cb.bscb.cornell.edu/nre/)
## JV: still a lot needs be done to make this general.
## last updated: 06.04.16 	vitti@broadinstitute.org

genmapfileprefix, genmapfilesuffix = "/idi/sabeti-scratch/ilya/gsvn/Data/Ilya_Data/genmaps/hm2/genetic_map_chr", "_b37.txt"
writefileprefix = "results_inclusive_test_filter"
chroms = range(1,23)
	
def main():
	writefilename = writefileprefix + "_recom_regions.txt"
	print "writing to " + writefilename
	writefile = open(writefilename, "w")
	for chrom in chroms:
		regionsfilename = "regions/" + writefileprefix + "_chr" + str(chrom) + ".txt"
		regionsfile = open(regionsfilename, 'r')
		print "reading from: " + regionsfilename
		firstline = regionsfile.readline()
		entries = firstline.split()
		chrom, start, end = int(entries[0]), int(entries[1]), int(entries[2])
		writefile.write("\n\n" + str(chrom) + "\n")
		genmapfilename = genmapfileprefix + str(chrom) + genmapfilesuffix
		genmapfile = open(genmapfilename, 'r')
		genmapfile.readline() #strip header

		for line in genmapfile:
			entries = line.split()
			pos, rate, mappos = int(entries[0]), float(entries[1]), float(entries[2])
			#found a region
			if pos >= start and pos < end:
				writefile.write(line)
				entries = line.split()
				thispos, thisrate, thismappos = int(entries[0]), float(entries[1]), float(entries[2])
				thiscosirate = thisrate / 10**8 # /10^8 to get cosi rate (cM/Mb -> M/b)

			if pos >= end: #next region
				writefile.write('\n')
				line2 = regionsfile.readline()
				if not line2:
					break
				entries = line2.split()
				chrom, start, end = int(entries[0]), int(entries[1]), int(entries[2])
		genmapfile.close()
		regionsfile.close()
	writefile.close()

main()