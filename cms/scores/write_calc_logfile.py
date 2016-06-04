## {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
## script to generate an excel spreadsheet with information on calculations performed
## lifted from /idi/sabeti-scratch/jvitti/synth/ 12.29.15; last updated 01.14.16
#-> JV: should use func_scores.py function to get filename
## 06.04.16 	vitti@broadinstitute.org

runLabel = "122915_fixed_genloc"
savename = "calc_log_" + runLabel + ".xls"
chroms = range(1,23)
chroms = [str(z) for z in chroms]
print "MUST GET NEW HAP VALS FOR X CHROMS; EXCLUDING FROM LOGFILE."
#chroms.append('X') #i.e. haploid X, excluding PARs
pops = [ 'YRI', 'ESN', 'GWD', 'LWK', 'MSL', 'ACB', 'ASW', 'PUR', 'CLM', 'MXL', 'PEL', 'CHB', 'CDX', 'CHS', 'JPT', 'KHV', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'BEB', 'GIH', 'ITU', 'PJL', 'STU']
masks = ['strictMask'] #'pilot_style_mask'
runtimes = ['062215_2pm_nomultialleles'] #'062115_3pm_multirow' 
mergesuffix = "merged_121915"
superpops = ['AFR', 'EAS', 'SAS', 'EUR', 'AMR', 'AFRAM']
superpopDict = {'BEB':'SAS', 'GIH':"SAS", 'ITU':'SAS', 'PJL':'SAS', 'STU':'SAS', 'CDX':'EAS', 'CHB':'EAS', 'CHS':'EAS', 'JPT':'EAS', 'KHV':'EAS', 'CEU':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR', 'TSI':'EUR', 'ESN':'AFR', 'GWD':'AFR', 'LWK':'AFR', 'MSL':'AFR', 'YRI':'AFR', 'ACB':'AFRAM', 'ASW':'AFRAM', 'PUR':'AMR', 'CLM':'AMR', 'MXL':'AMR', 'PEL':'AMR'}
superpopColorDict = {'AMR':'red', 'AFR':'gold', 'AFRAM':'light_orange', 'SAS':'plum', 'EAS':'lime', 'EUR':'turquoise'}
num_X = {'ACB':145, 'ASW':96, 'BEB':130, 'CDX':142, 'CEU':149, 'CHB':160, 'CHS':158, 'CLM':145, 'ESN':145, 'FIN':160, 'GBR':136, 'GIH':150, 'GWD':171, 'IBS':160, 'ITU':145, 'JPT':152, 'KHV':152, 'LWK':154, 'MSL':128, 'MXL':96, 'PEL':129, 'PJL':144, 'PUR':154, 'STU':149, 'TSI':161, 'YRI':164}
# taken from: /idi/sabeti-scratch/jvitti/CMS_1KG_P3/vcfs/integrated_call_samples_v3.20130502.ALL.panel
num_indivs = {'GWD': 113, 'ACB': 96, 'BEB': 86, 'PEL': 85, 'LWK': 99, 'MSL': 85, 'GBR': 91, 'IBS': 107, 'ASW': 61, 'TSI': 107, 'KHV': 99, 'CEU': 99, 'YRI': 108, 'CHB': 103, 'STU': 102, 'CHS': 105, 'ESN': 99, 'FIN': 99, 'GIH': 103, 'PJL': 96, 'MXL': 64, 'PUR': 104, 'ITU': 102, 'CDX': 93, 'JPT': 104, 'CLM': 94}

from func_files import *
from tempfile import TemporaryFile
from xlwt import Workbook, easyxf
import subprocess
import os.path
import random
import gzip

###################
## DEFINE COLORS ##
###################

boldstyle = easyxf('font: bold 1;')
wrapstyle = easyxf('alignment: wrap on, vert center, horiz center')
AFR_style = easyxf("pattern: pattern solid, fore_colour red; alignment: wrap on, vert center, horiz center")
EAS_style = easyxf("pattern: pattern solid, fore_colour lime; alignment: wrap on, vert center, horiz center")
EUR_style = easyxf("pattern: pattern solid, fore_colour turquoise; alignment: wrap on, vert center, horiz center")
SAS_style = easyxf("pattern: pattern solid, fore_colour plum; alignment: wrap on, vert center, horiz center")
AFRAM_style = easyxf("pattern: pattern solid, fore_colour light_orange; alignment: wrap on, vert center, horiz center")
AMR_style = easyxf("pattern: pattern solid, fore_colour gold; alignment: wrap on, vert center, horiz center")

##########
## MAIN ##
##########

def main():
	book = Workbook()
	sheet1 = book.add_sheet('datafiles') #vcf, tped
	sheet2 = book.add_sheet('single-pop calcs') #iHS, delDAF, delIHH
	sheet3 = book.add_sheet('xpop calcs') #XP-EHH, Fst

	#########################
	## sheet 1 - datafiles ##
	#########################

	headers = ['superpop', 'pop', 'chr', 'nchroms', 'nsnps', 'mask', 'runtime', 'vcf', 'tped']
	for i in range(len(headers)):
		sheet1.write(0, i, headers[i], boldstyle)

	#tpedlocs, vcflocs = {}, {}
	for mask in masks:
		for runtime in runtimes:

			if runtime == "062215_2pm_nomultialleles" and mask == "strictMask":
				#to save time
				num_snps = {'X': 1700999, '15': 1448337, '1': 3995280, '2': 4620766, '16': 1530805, '3': 3881176, '17': 1246051, '4': 3741568, '18': 1555683, '5': 3436942, '19': 671024, '6': 3272818, '20': 1170888, '7': 2823663, '21': 695112, '22': 596753, '8': 3036057, '9': 2240022, '10': 2477005, '11': 2482524, '12': 2412835, '13': 1926052, '14': 1694284} #taken from wc -l from metadata tped files
			else:
				num_snps = {}
				autosomes = range(1,23)
				for chrom in autosomes:
					sitefile = "/idi/sabeti-scratch/jvitti/cms_venv/tpeds/v5a_chr" + str(chrom) + "_" + mask + "_" + random.choice(pops) + "_" + runtime + ".tped.allele_metadata.gz"
					with gzip.open(sitefile) as f:
						for i, line in enumerate(f): #one for header cancels out one ignored at end of for-loop
							pass
					num_snps[chrom] = i

			irow = 1
			for pop in pops:
				superpop = superpopDict[pop]
				for chrom in chroms:
					if chrom == "X":
						vers = "v1a"
						nindivs = num_X[pop]
					else:
						vers = "v5a"
						chrom = int(chrom) #for excel
						nindivs = num_indivs[pop]

					nsnps = num_snps[str(chrom)]

					#OPTIONAL: can store in dict to make specific retrieval easier
					#identifier = (pop, chrom, mask, runtime)
					
					#vcf
					vcffilename = "/idi/sabeti-scratch/jvitti/cms_venv/vcfs/" + vers + "_chr" + str(chrom) + "_" + mask + "_global.vcf.gz"
					vcfloc = checkFileExists(vcffilename)
					
					#tped
					tpedfilename = "/idi/sabeti-scratch/jvitti/cms_venv/tpeds/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + ".tped.gz"
					tpedloc = checkFileExists(tpedfilename)

					#tpedlocs[identifier] = tpedloc
					sheet1.write(irow, 0, superpop, eval(superpop + "_style"))
					sheet1.write(irow, 1, pop, wrapstyle)
					sheet1.write(irow, 2, chrom, wrapstyle)
					#num indivs
					sheet1.write(irow, 3, nindivs, wrapstyle)
					sheet1.write(irow, 4, nsnps, wrapstyle)
					#num snps
					sheet1.write(irow, 5, mask, wrapstyle)
					sheet1.write(irow, 6, runtime, wrapstyle)
					sheet1.write(irow, 7, vcfloc, wrapstyle)
					sheet1.write(irow, 8, tpedloc, wrapstyle)
					irow += 1
	#format cells for locations
	sheet1.col(6).width = 256 * 24
	sheet1.col(7).width = 256 * 40 #~x char wide
	sheet1.col(8).width = 256 * 50

	################################
	## sheet 2 - single-pop calcs ##
	################################

	headers = ['pop', 'chr', 'ihs', 'delIHH', "pop_IHH", "ihs_norm", 'delIHH_norm']#, "nsl", "nsl_norm"] #delIHH_norm is from XP-EHH calcs
	for i in range(len(headers)):
		sheet2.write(0, i, headers[i], boldstyle)
	mask, runtime = "strictMask", "062215_2pm_nomultialleles" #only calculations for this dataset
	iHSlocs, delDAFlocs, delIHHlocs = {}, {}, {}
	irow = 1
	for pop in pops:
		for chrom in chroms:
			if chrom == "X":
				vers = "v1a"
				nindivs = num_X[pop]
				ihsfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + str(nindivs) + ".ihs.out"
				delihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + str(nindivs) +".delihh.out"
				ihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/xp_run/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + str(num_X[pop]) + "_" + runtime + "_xpehh.out.ihh"
				ihsnormfilename = ihsfilename + "_6col.100bins.norm"
				delihhnormfilename = '/idi/sabeti-scratch/jvitti/cms_venv/delIHH_normed/' + pop + "_" + mask + "_" + runtime + "_xchr_" + str(num_X[pop]) + ".delihh.100bins.norm"
				#nslfilename = "MISSING"
				#nslnormfilename = "MISSING"

			else:
				vers = "v5a"
				chrom = int(chrom) #for excel
				ihsfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/"+ vers + "_autosomes_" + mask + "_" + pop + "_" + runtime + "_" + mergesuffix + ".ihs.out"
				ihsnormfilename = ihsfilename + ".100bins.norm.chr" + str(chrom) + ".txt"
				delihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/merged/"+ vers + "_autosomes_" + mask + "_" + pop + "_" + runtime + "_" + mergesuffix + ".delihh.out"
				delihhnormfilename = delihhfilename + ".100bins.norm.chr" + str(chrom) + ".txt"
				if  chrom <= 20:
					ihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/merged/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_" + runtime + "_" + mergesuffix + ".xpehh.out.ihh"
				elif chrom > 20: #did not have to remake these
					ihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/xp_run/" + vers + "_chr" + str(chrom) + "_" + mask + "_" + pop + "_062215_2pm_nomultialleles_xpehh.out.ihh"
				#nslfilename = "MISSING"
				#nslnormfilename = "MISSING"

			#identifier = (pop, chrom)
			#ihslocs[identifier] = ihsloc

			ihsloc = checkFileExists(ihsfilename)
			delihhloc = checkFileExists(delihhfilename)
			ihhloc = checkFileExists(ihhfilename)
			ihsnormloc = checkFileExists(ihsnormfilename)
			delihhnormloc = checkFileExists(delihhnormfilename)
			#nslloc = checkFileExists(nslfilename)
			#nslnormloc = checkFileExists(nslnormfilename)

			sheet2.write(irow, 0, pop, eval(superpopDict[pop] + "_style"))
			sheet2.write(irow, 1, chrom, wrapstyle)
			sheet2.write(irow, 2, ihsloc, wrapstyle)
			sheet2.write(irow, 3, delihhloc, wrapstyle)
			sheet2.write(irow, 4, ihhloc, wrapstyle)
			sheet2.write(irow, 5, ihsnormloc, wrapstyle)
			sheet2.write(irow, 6, delihhnormloc, wrapstyle)
			#sheet2.write(irow, 7, nslloc, wrapstyle)
			#sheet2.write(irow, 8, nslnormloc, wrapstyle)


			irow +=1
	#format cells for locations
	sheet2.col(2).width = 256 * 50
	sheet2.col(3).width = 256 * 50 #~x char wide
	sheet2.col(4).width = 256 * 50
	sheet2.col(5).width = 256 * 50 #~x char wide
	sheet2.col(6).width = 256 * 50
	sheet2.col(7).width = 256 * 50 #~x char wide
	sheet2.col(8).width = 256 * 50

	##########################
	## sheet 3 - xpop calcs ##
	##########################

	headers = ['pop1', 'pop2', 'chrom', 'fst, delDAF', 'xp-ehh', 'xp-ehh_norm']
	for i in range(len(headers)):
		sheet3.write(0, i, headers[i], boldstyle)
	mask, runtime = "strictMask", "062215_2pm_nomultialleles" #only calculations for this dataset
	#xpehhlocs, Fstlocs = {}, {}
	popPairs = []
	for ipop in range(len(pops)):
		for jpop in range(ipop, len(pops)):
			if ipop != jpop:
				popPairs.append((pops[ipop], pops[jpop]))

	irow = 1
	for popPair in popPairs:
		for chrom in chroms:
			if chrom == "X":
				vers = "v1a"
			else:
				vers = "v5a"
				chrom = int(chrom) #for excel formatting

			identifier = (popPair, chrom, mask, runtime)

			pop1, pop2 = popPair[0], popPair[1] #must catch possibility of misordered
			nsamp1, nsamp2 = num_X[pop1], num_X[pop2]
			#must also account for haploid X samplesizes 

			#Fst
			if chrom == "X":
				fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/xchr" + "/" + pop1 + "_" + str(nsamp1) + "_" + pop2 + "_" + str(nsamp2) + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
			else:
				fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop1 + "_" + pop2 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
			fstloc = checkFileExists(fstfilename)
			if fstloc == "MISSING":
				pop1, pop2 = popPair[1], popPair[0] #must catch possibility of misordered
				nsamp1, nsamp2 = num_X[pop1], num_X[pop2]
				if chrom == "X":
					fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/xchr" + "/" + pop1 + "_" + str(nsamp1) + "_" + pop2 + "_" + str(nsamp2) + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
				else:
					fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop1 + "_" + pop2 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
				fstloc = checkFileExists(fstfilename)


			#XP-EHH
			if chrom == 'X':
				xpehhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/" + vers + "_chrX_" + mask + "_" + runtime + "_" + pop1 + "_" + str(num_X[pop1]) + "_" + pop2 + "_" + str(num_X[pop2]) + ".xpehh"
				xpehhloc = checkFileExists(xpehhfilename)
				if xpehhloc == "MISSING":
					xpehhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/" + vers + "_chrX_" + mask + "_" + runtime + "_" + pop2 + "_" + str(num_X[pop2]) + "_" + pop1 + "_" + str(num_X[pop1]) + ".xpehh"
					xpehhloc = checkFileExists(xpehhfilename)

			else:
				xpehhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/" + vers + "_autosomes_" + mask  + "_" + runtime +  "_" + pop1 + "_" + pop2 + ".xpehh"
				xpehhloc = checkFileExists(xpehhfilename)
				if xpehhloc == "MISSING":
					xpehhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/" + vers + "_autosomes_" + mask  + "_" + runtime +  "_" + pop2 + "_" + pop1 + ".xpehh"
					xpehhloc = checkFileExists(xpehhfilename)

			if chrom == "X":
				xpehh_normedfile = xpehhloc + ".norm"

			else:
				xpehh_normedfile = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask + "_" + runtime + "_" + mergesuffix + "_" + pop1 + "_" + pop2 + ".xpehh.norm.chr" + str(chrom) + ".txt"
				#xpehh_normedfileloc = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask + "_" + runtime + "_" + mergesuffix + "_" + pop1 + "_" + pop2 + ".xpehh.norm" + 
				#xpehh_normedfile = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/bychrom/" + vers + "_chr" + str(chrom) + "_" + mask +"_" + runtime +"_" + pop1 +"_" + pop2 +".xpehh.norm"
				xpehh_normedfileloc = checkFileExists(xpehh_normedfile)
				if xpehh_normedfileloc == "MISSING":
					xpehh_normedfile = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask + "_" + runtime + "_" + mergesuffix + "_" + pop2 + "_" + pop1 + ".xpehh.norm.chr" + str(chrom) + ".txt"
					xpehh_normedfileloc = checkFileExists(xpehh_normedfile)

			sheet3.write(irow, 0, pop1, eval(superpopDict[pop1] + "_style"))
			sheet3.write(irow, 1, pop2, eval(superpopDict[pop2] + "_style"))		
			sheet3.write(irow, 2, chrom, wrapstyle)
			sheet3.write(irow, 3, fstloc, wrapstyle)
			sheet3.write(irow, 4, xpehhloc, wrapstyle)
			sheet3.write(irow, 5, xpehh_normedfileloc, wrapstyle)
			irow += 1

	#format cells for locations
	sheet3.col(0).width = 256 * 10
	sheet3.col(1).width = 256 * 10
	sheet3.col(2).width = 256 * 10
	sheet3.col(3).width = 256 * 50
	sheet3.col(4).width = 256 * 50
	sheet3.col(5).width = 256 * 80

	book.save(savename)
	book.save(TemporaryFile())

	logfileloc = "~/public_html/synth/"
	mvcommand = "mv " + savename + " " + logfileloc
	print mvcommand
	#not sure why it's not letting me automatically save to my public html and must move manually.
	#givecommand(mvcommand)
	print "Saved file: http://www.broadinstitute.org/~vitti/synth/" + savename

main()