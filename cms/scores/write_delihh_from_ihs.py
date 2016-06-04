## {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
## original /idi/sabeti-scratch/jvitti/cms_venv/write_delIHH_file.py 12.26.15
## call program as: python write_delihh_from_ihs.py {POP}
## last updated 06.04.16 	vitti@broadinstitute.org

suffix = "merged_121915"
mask, runtime, vers = "strictMask", "062215_2pm_nomultialleles", "v5a" 
chroms = range(1,23)

from math import fabs
import sys

def write_delIHH_file(readfilename, writefilename):
	"""given a selscan iHS file, parses it and writes an analogous file with delIHH information"""
	readfile = open(readfilename, 'r')
	writefile = open(writefilename, 'w')
	for line in readfile:
		entries = line.split()
		locus, phys, freq_1, ihh_1, ihh_0, ihs, ihh_0_l, ihh_0_r, ihh_1_l, ihh_1_r = entries #[0:10]
					#ancestral - derived
		unstand_delIHH = fabs(float(ihh_1) - float(ihh_0))
		writeline = locus + "\t" + phys + "\t" + freq_1 + "\t" + str(ihs) + "\t" + str(unstand_delIHH) +  "\n"
		writefile.write(writeline)
	writefile.close()
	readfile.close()

def main():
	pop = sys.argv[1]
	for chrom in chroms:
		readfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime +"_" + suffix  + ".ihs.out"
		writefilename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + suffix +  ".delihh.out"
		write_delIHH_file(readfilename, writefilename)
		print "wrote file to: " +  writefilename
main()