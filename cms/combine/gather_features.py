## Given a region, pulls values for selection summary statistics (e.g. iHS, Fst etc) to populate  files for use with CMS_ML. See e.g. AJin documentation. 
## /idi/sabeti-scratch/jvitti/synth/ original 12.29.15; last updated 1.13.16

#towriteprefix = "control_LCT_010816"
#chrom = 2
#startPos, endPos = 135650000, 136650000
towriteprefix = "test_snpdropout_ERLIN1_010816"
chrom = 10
startPos, endPos = 101614433, 102614433
pops = ['YRI', 'CEU', 'CHB'] 
writedir = "/idi/sabeti-scratch/jvitti/synth/ml_regions/"
vers, mask, runtime, mergesuffix = 'v5a', 'strictMask', '062215_2pm_nomultialleles', "merged_121915"
num_indivs = {'GWD': 113, 'ACB': 96, 'BEB': 86, 'PEL': 85, 'LWK': 99, 'MSL': 85, 'GBR': 91, 'IBS': 107, 
'ASW': 61, 'TSI': 107, 'KHV': 99, 'CEU': 99, 'YRI': 108, 'CHB': 103, 'STU': 102, 'CHS': 105, 'ESN': 99, 
'FIN': 99, 'GIH': 103, 'PJL': 96, 'MXL': 64, 'PUR': 104, 'ITU': 102, 'CDX': 93, 'JPT': 104, 'CLM': 94}

from func import *

############
### MAIN ###
############

for pop in pops:
	towritefilename = writedir + towriteprefix + "_chr" + str(chrom) + "_" + str(startPos) + "_" + str(endPos) + "_" + str(pop) + ".txt"
	print "Writing to file: " + towritefilename

	############################
	### LOAD SCORES - PER POP ##
	############################

	if chrom <= 20:
		ihsfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/"+ vers + "_chr" + str(chrom) + "_" + mask +  "_" + pop + "_" + runtime + "_" + mergesuffix + ".ihs.out"
	else:
		ihsfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/"+ vers + "_chr" + str(chrom) + "_" + mask + "_" + runtime + ".ihs.out"
	ihsnormfilename = "/idi/sabeti-scratch/jvitti/cms_venv/ihs/merged/"+ vers + "_autosomes_" + mask + "_" + pop + "_" + runtime + "_" + mergesuffix + ".ihs.out.100bins.norm.chr" + str(chrom) + ".txt"
	delihhfilename = "/idi/sabeti-scratch/jvitti/cms_venv/delIHH/merged/"+ vers + "_autosomes_" + mask + "_" + pop + "_" + runtime + "_" + mergesuffix + ".delihh.out"
	delihhnormfilename = delihhfilename + ".100bins.norm.chr" + str(chrom) + ".txt"

	for filename in [ihsfilename, ihsnormfilename, delihhfilename, delihhnormfilename]:
		if checkFileExists(filename) == "MISSING":
			print "Error: missing file " 
			print filename

	all_locus_a, all_phys_a, all_freq_1_a, all_ihh_1_a, all_ihh_0_a, all_ihs_a, all_ihh_0_l, all_ihh_0_r, all_ihh_1_l, all_ihh_1_r = readihsfile_fullinfo(ihsfilename, startPos, endPos)
	all_locus_b, all_phys_b, all_freq_1_b, all_ihh_1_b, all_ihh_0_b, all_ihs_b, all_normalized_ihs, all_finalvals_b = readihsfile_normed(ihsnormfilename, startPos, endPos)
	all_locus_c, all_phys_c, all_freq_1_c, all_ihs_c, all_unstand_delihh, all_normalized_delihh, all_finalvals_c =readdelihhfile_normed(delihhnormfilename, startPos, endPos)

	assert all_phys_a == all_phys_b and all_phys_b == all_phys_c
	assert all_ihs_b == all_ihs_a and all_ihs_b == all_ihs_c
	print "\n"

	#####################################
	### LOAD SCORES - POP COMPARISONS ###
	#####################################
	otherPops = pops[:]
	otherPops.remove(pop)
	comparisons = {}
	FstCompOrders, XpCompOrders = [], [] #keep track of whether putative sel pop is encoded as 'A' or 'B'
	for otherPop in otherPops:

		pop1, pop2 = pop, otherPop
		fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop1 + "_" + pop2 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
		xpehh_normedfile = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask + "_" + runtime + "_" + mergesuffix + "_" + pop1 + "_" + pop2 + ".xpehh.norm.chr" + str(chrom) + ".txt"
		
		if not os.path.isfile(fstfilename):
			FstCompOrders.append((pop2, pop1))
			fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop2 + "_" + pop1 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
		else:
			FstCompOrders.append((pop1, pop2))

		if not os.path.isfile(xpehh_normedfile): 
			XpCompOrders.append((pop2, pop1))
			xpehh_normedfile = "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/" + vers + "_autosomes_" + mask + "_" + runtime + "_" + mergesuffix + "_" + pop2 + "_" + pop1 + ".xpehh.norm.chr" + str(chrom) + ".txt"
		else:
			XpCompOrders.append((pop1, pop2))

		for filename in [fstfilename, xpehh_normedfile]:
			if not checkFileExists(filename):
				print "Error: missing file " 
				print filename

		all_locus_d, all_phys_d, all_genpos, all_freq_A_1, all_ihhA, all_freq_B_1, all_ihhB, allunstandXP, allstandXP, allfinalvals = readxpehhfile(xpehh_normedfile, startPos, endPos)
		all_phys_e, all_fst, all_deldaf = readfstdeldaffile(fstfilename, startPos, endPos)
		for item in ['all_locus_d', 'all_phys_d', 'all_genpos', 'all_freq_A_1', 'all_ihhA', 'all_freq_B_1', 'all_ihhB', 'allunstandXP', 'allstandXP', 'allfinalvals', 'all_phys_e', 'all_fst', 'all_deldaf']:
			thiskey = (otherPop, item)
			comparisons[thiskey] = eval(item)

		assert len(otherPops) == 2 #else, LSBL needs reconfigured
		pop1, pop2 = otherPops
		fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop1 + "_" + pop2 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
		if not os.path.isfile(fstfilename):
			fstfilename = "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr" + str(chrom) + "/" + pop2 + "_" + pop1 + "_chr" + str(chrom) + "_" + runtime + "_" + mask + "_fst_daf.txt"
		outgroups_phys, outgroups_fst, outgroups_deldaf = readfstdeldaffile(fstfilename, startPos, endPos)

	########## SHOULD I MAYBE MOVE THIS ELSEWHERE?
	#calc_TajimasD(filename, startpos, endpos, calcPop = 'A', windowSize = 1000):
	print "\n\n\t******D*******:\n"
	nIndivs = num_indivs[pop]
	return_phys, return_d = calc_TajimasD(xpehh_normedfile, startPos, endPos, nIndivs)
	##########

	######################
	### COLLATE SCORES ###
	######################
	common_phys_pos = set(all_phys_a) & set(all_phys_b) & set(all_phys_c) & set(all_phys_d) & set(all_phys_e) & set(return_phys)
	common_phys_pos = list(set(common_phys_pos))
	common_phys_pos.sort()

	######################
	### WRITE TO FILE ###
	######################

	print "Writing to: " + towritefilename
	towritefile = open(towritefilename, 'w')
	header = "rs\tPos\tgdPos\tiHS\tdeltaIHH\tmeanFst\tXPEHH\tdeltaDAF\tderFreq\tmeanAnc\tLeft iHH_A\tLeft iHH_D\tRight iHH_A\tRight iHH_D\tBoth iHH_A\tBoth iHH_D\tLSBL\tTajimasD\n"
	towritefile.write(header)
	missing_fst = 0 # ValueError; index_fst = outgroups_phys.index(common_phys_pos[i])
	for i in range(len(common_phys_pos)):

		index_a = all_phys_a.index(common_phys_pos[i])
		index_b = all_phys_b.index(common_phys_pos[i])
		index_c = all_phys_c.index(common_phys_pos[i])		
		index_d = all_phys_d.index(common_phys_pos[i])
		index_e = all_phys_e.index(common_phys_pos[i])
		index_tajima = return_phys.index(common_phys_pos[i])
		locus_id = str(chrom) + "_" + str(common_phys_pos[i])

		fst_comps, xpehh_comps = [], []
		der_freqs, anc_freqs = {}, {}
		for j in range(len(otherPops)):
			otherPop = otherPops[j]
			FstCompOrder = FstCompOrders[j]
			XpCompOrder = XpCompOrders[j]

			index_otherpop = comparisons[(otherPop, 'all_phys_e')].index(common_phys_pos[i])
			this_fst = comparisons[(otherPop, 'all_fst')][index_otherpop]
			fst_comps.append(this_fst)

			xp_index_otherpop = comparisons[(otherPop, 'all_phys_d')].index(common_phys_pos[i])
			this_xp = comparisons[(otherPop, 'allstandXP')][xp_index_otherpop]
			if otherPop == XpCompOrder[0]: #selpop is here coded as outgroup; must reverse.
				print "double-check xpehh reclaculation."
				temp = math.exp(float(this_xp))
				temp2 = temp**-1
				this_xp = math.log(temp2)
			xpehh_comps.append(this_xp)

			index_A_1 = comparisons[(otherPop, 'all_phys_d')].index(common_phys_pos[i])
			this_freq_A_1 = float(comparisons[(otherPop, 'all_freq_A_1')][index_A_1])
			this_freq_B_1 = float(comparisons[(otherPop, 'all_freq_B_1')][index_A_1])

			if otherPop == FstCompOrder[0]: #pop A is the outgroup
				der_freqs[otherPop] = (1 - this_freq_A_1)
				anc_freqs[otherPop] = this_freq_A_1
				der_freqs[pop] = (1 - this_freq_B_1)
				anc_freqs[pop] = this_freq_B_1

			else: #pop B is the outgroup
				der_freqs[otherPop] = (1 - this_freq_B_1)
				anc_freqs[otherPop] = this_freq_B_1
				der_freqs[pop] = (1 - this_freq_A_1)
				anc_freqs[pop] = this_freq_A_1

		meanAnc = np.average([anc_freqs[population] for population in otherPops])
		meanDer = np.average([der_freqs[population] for population in otherPops])
		DerFreq = 1 - all_freq_1_a[index_a]
		delDaf = DerFreq - meanDer
		maxxp = max(xpehh_comps)
		meanFst = np.average(fst_comps)

		if common_phys_pos[i] in outgroups_phys:
			index_fst = outgroups_phys.index(common_phys_pos[i])
			outgroup_fst_comp = float(outgroups_fst[index_fst])
			lsbl = (sum(fst_comps) - outgroup_fst_comp)/2

			towriteline = locus_id + "\t" +  str(common_phys_pos[i]) + "\t" + str(all_genpos[index_d]) + "\t" + str(all_normalized_ihs[index_b]) + "\t" + str(all_normalized_delihh[index_c]) + "\t" + str(meanFst) + "\t" + str(maxxp) + "\t"
			towriteline += str(delDaf) + "\t" + str(DerFreq) + "\t" + str(meanAnc) + "\t" + str(all_ihh_1_l[index_a]) + "\t" + str(all_ihh_0_l[index_a]) + "\t"
			towriteline += str(all_ihh_1_r[index_a]) + "\t" + str(all_ihh_0_r[index_a]) + "\t" + str(all_ihh_1_a[index_a]) + "\t" + str(all_ihh_0_a[index_a]) + "\t"
			towriteline += str(lsbl) + "\t" + str(return_d[index_tajima]) + "\n"
			#print towriteline
			towritefile.write(towriteline)
		else:
			# ask andrew what he did here?for now, skip.		
			missing_fst +=1

	print "passed on lsbl calcs for " + str(missing_fst)
	print "this is worth investigating some more."
	towritefile.close()
