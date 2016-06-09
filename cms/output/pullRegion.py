# /idi/sabeti-scratch/jvitti/synth
# called in parallel by visualize_haps.py
# 01.13.16 updated 05.27.16 to pull from YRI_IPSC-filtered VCFs.

from func_viz import writeHapFile
from func_parse import *
import sys

if len(sys.argv) != 8:
	print "call program as: python pullRegion.py {dataSource} {infilename} {chrom} {startpos} {endpos} {pop} {prefix}"
	sys.exit()

dataSource = sys.argv[1]
infilename = sys.argv[2]
chrom = sys.argv[3]
startPos = int(sys.argv[4])
endPos = int(sys.argv[5])
pop = sys.argv[6]
prefix = sys.argv[7]

if dataSource == "unmaskedvcf" or dataSource == "maskedvcf" or datasource == "yri_ipsc":
	varids, genpositions, physpositions, all_genotypes = readVcfRegion(infilename, [pop], startPos, endPos)
elif dataSource == "filteredtped":
	varids, genpositions, physpositions, all_genotypes = readTpedRegion(infilename, startPos, endPos)
else:
	print "dataSource must be set to one of: {unmaskedvcf}, {maskedvcf}, {filteredtped}, {yri_ipsc}"
	sys.exit()

print "loaded genotypes for " + str(len(all_genotypes)) + " sites... "	
writeHapFile(all_genotypes, varids, chrom, prefix + "_" + str(pop) + "_" + dataSource)