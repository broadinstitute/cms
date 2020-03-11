##	/idi/sabeti-scratch/jvitti/remodel/run/ ## last used 11.17.2019
##	02.20.2019	## CONFIGURED FOR THREE POPS NB, n=200 #3.3 ->> 4pops

stemdir = "/idi/sabeti-scratch/jvitti/remodel/run4/"

#DO GRAVEL SEPARATELY

import os
import sys
import gzip
import subprocess

def get_headerline(nPerPop, tped_base):
	headerline = "##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
	if "gravel" in tped_base:
		pops = [1,2,3]
	else:
		pops = [1,2,3,4]
	for pop in pops:
		for i in range(int(nPerPop/2)):
			this_indiv ="pop" + str(pop) + "_" + str(i)
			headerline += this_indiv + "\t"
	headerline = headerline[:-1] 
	return headerline
def rewrite_from_tped_add_og(tped_base, vcf_filename):
	pops = [1, 2, 3, 4]
	if "gravel" in tped_base:
		pops.remove(4)
	#print(pops)
	for pop in pops:
		sample_file = tped_base + str(pop) + ".tped"
		bgzip_cmd = "bgzip " + sample_file
		if not os.path.isfile(sample_file + ".gz"):
			subprocess.check_output(bgzip_cmd.split())

	if os.path.isfile(vcf_filename + ".gz"):
		rm_cmd = "rm " + vcf_filename + ".gz"
		subprocess.check_output( rm_cmd.split() )

	if True: #rewrite add pop 4
	#if not os.path.isfile(vcf_filename + ".gz"):
		vcf_file = open(vcf_filename, 'w')
		sample_file = tped_base +"1.tped.gz"
		tped = gzip.open(sample_file, 'rt')
		for line in tped:
			entries = line.split()
			physpos = int(entries[3])
			genotypes = entries[4:]
			these_reformatted_genotypes, nIndivs = reformat_genotypes_as_diploid_string(genotypes)
			break
		tped.close()
		nPerPop = nIndivs
		
		headerline = get_headerline(nPerPop, tped_base)	
		vcf_file.write(headerline + "\n")
		writelines = {}
		for pop in pops:
			tped_filename = tped_base + str(pop) + ".tped.gz"
			#tped = open(tped_filename, 'r')
			tped = gzip.open(tped_filename, 'rt')
			for line in tped:
				entries = line.split()
				physpos = int(entries[3])
				genotypes = entries[4:]
				these_reformatted_genotypes, nIndivs = reformat_genotypes_as_diploid_string(genotypes)
				if physpos in writelines.keys():
					writelines[physpos] +=  "\t" + these_reformatted_genotypes
				else:
					writelines[physpos] = these_reformatted_genotypes
			tped.close()
		all_lines = list(writelines.keys())
		all_lines.sort()
		for item in all_lines:
			full_writeline = "1\t" + str(item) + "\t" + str(item) + "\t0\t1\t100\tPASS\t.\tGT\t" + str(writelines[item]) + "\n"
			vcf_file.write(full_writeline)
		vcf_file.close()
		print('wrote to ' + vcf_filename)
		bgzip_cmd = "bgzip " + vcf_filename
		print(bgzip_cmd)
		subprocess.check_output(bgzip_cmd.split())
		tabix_cmd = "tabix " + vcf_filename + ".gz"
		print(tabix_cmd)
		subprocess.check_output(tabix_cmd.split())
	return
def reformat_genotypes_as_diploid_string(genotypes):
	these_reformatted_genotypes = ""
	assert(len(genotypes) == 200 or len(genotypes) == 172)
	for i in range(0,len(genotypes),2):
		writegeno = ''
		this_A, this_B = genotypes[i], genotypes[i+1]
		if this_A == "1":
			writegeno +="0"
		elif this_A == "0":
			writegeno += "1"
		writegeno += "|"
		if this_B == "1":
			writegeno +="0"
		elif this_B == "0":
			writegeno += "1"
		these_reformatted_genotypes += writegeno + "\t"
	these_reformatted_genotypes = these_reformatted_genotypes[:-1]
	return these_reformatted_genotypes, len(genotypes)
def get_rank(causalPos, positions, values):
	if int(causalPos) not in positions:
		causal_rank = 99999999
	else:
		causalIndex = positions.index(int(causalPos))
		causalValue = values[causalIndex]
		values.sort()
		values.reverse()
		causal_rank = values.index(causalValue)
	return causal_rank
def load_isafe_local(filename):
	positions, dafs, scores = [], [], []
	openfile = open(filename, 'r')
	openfile.readline()
	for line in openfile:
		entries = line.split('\t')
		position = int(entries[0])
		daf = float(entries[2])
		score = float(entries[1])
		positions.append(position)
		dafs.append(daf)
		scores.append(score)
	openfile.close()
	return positions, dafs, scores
def get_analytics(this_filename, causalPos = 1500000):
	positions, dafs, scored = load_isafe_local(this_filename)
	n_scored = len(positions)
	this_rank = get_rank(causalPos, positions, scored)
	return this_rank, n_scored


def main():
	model, regime, pop, irep = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
	indivsdir = "/idi/sabeti-scratch/jvitti/remodel/run/indivfiles/"
	if model == "gravel":
		indivsdir += "threepop_model/"
	basedir =stemdir + model + "_" + regime +"/"#+ "_sel" + pop + "/"
	this_tpeddir = basedir + "tpeds/"
	this_isafedatadir = basedir + "isafe_data/"
	this_isafedir = basedir + "isafe/"
	
	##############################
	## CONVERT TO VCF AND INDEX ##
	##############################
	vcf_filename = this_isafedatadir + "rep" + str(irep) + "_0_.vcf"
	#if True: #REDO, ADD POP 4!!
	#if not os.path.isfile(vcf_filename + ".gz.tbi") or not os.path.isfile(vcf_filename + '.gz'):
	if True:
		tped_base = this_tpeddir + "rep" + str(irep) + "_0_"
		rewrite_from_tped_add_og(tped_base, vcf_filename)
	vcf_filename = vcf_filename +".gz"

	###############
	## RUN ISAFE ## 
	###############
	
	sample_case_filename, sample_cont_filename = indivsdir + "dummy_sim_indivs_pop" + str(pop) + "_n100.tsv.gz", indivsdir + "dummy_sim_indivs_pop" + str(pop) + "_out_n100.tsv.gz"
	savefilename = this_isafedir + "OG_rep" + str(irep) + "_0_pop" + str(pop) 
	#iSafe_cmd = #"python2.7 /idi/sabeti-scratch/jvitti/iSAFE/src/isafe.py --input " + vcf_filename + " --region 1:0-5000000 --output " + savefilename + " --MaxRank 301 --MaxFreq 1 --AA /idi/sabeti-scratch/jvitti/dummy_sim_ancestral.fa --vcf-cont " + vcf_filename + " --sample-cont " + sample_cont_filename + " --sample-case " + sample_case_filename
	iSafe_cmd = "python /idi/sabeti-scratch/jvitti/iSAFE/src/isafe.py --input " + vcf_filename + " --region 1:0-5000000 --output " + savefilename + " --MaxRank 301 --MaxFreq 1 --AA /idi/sabeti-scratch/jvitti/dummy_sim_ancestral.fa --vcf-cont " + vcf_filename + " --sample-cont " + sample_cont_filename + " --sample-case " + sample_case_filename
	if not os.path.isfile(savefilename + ".iSAFE.out"):
		print(iSafe_cmd)
		subprocess.check_output( iSafe_cmd.split() )
		
	##############
	## LOG RANK ##
	##############
	#rankfilename = savefilename + ".rank"
	#if not os.path.isfile(rankfilename):
	#	thisRank, nTotal = get_analytics(savefilename + ".iSAFE.out")
	#	writefile = open(rankfilename, 'w')
	#	writeline = str(thisRank) + "\t" + str(nTotal) + "\n"
	#	writefile.write(writeline)
	#	writefile.close()
	#	print('wrote to ' + rankfilename)

main()
