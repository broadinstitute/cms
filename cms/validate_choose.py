#pull values from .cms files (choose) and view score distribution **PASS .CMS FILES TO LIKES_FROM_SCORES
# 11.4.16
import os
import matplotlib as mp 
mp.use('TKAgg') 
import matplotlib.pyplot as plt
import numpy as np

examinepos = 916842

pairfiles = ['/idi/sabeti-scratch/jvitti/scores_composite/nulldefault/neut/rep1_1_2.pair', '/idi/sabeti-scratch/jvitti/scores_composite/nulldefault/neut/rep1_1_3.pair', '/idi/sabeti-scratch/jvitti/scores_composite/nulldefault/neut/rep1_1_4.pair']
#'/idi/sabeti-scratch/jvitti/scores_composite2/nulldefault/neut/rep1_1.cms.out'

def read_cms_repfile(infilename):
	physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = [], [], [], [], [], [], [], [], []
	if not os.path.isfile(infilename):
		return physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed
	openfile = open(infilename, 'r')
	for line in openfile:
		entries = line.split()
		physpos.append(int(entries[0]))
		genpos.append(float(entries[1]))
		ihs_normed.append(float(entries[2]))
		delihh_normed.append(float(entries[3]))
		xpehh_normed.append(float(entries[4]))
		fst.append(float(entries[5]))
		deldaf.append(float(entries[6]))
		cms_unnormed.append(float(entries[7]))
		if ".norm" in infilename:
			cms_normed.append(float(entries[8]))
	openfile.close()
	return physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed

def read_cms_pairfile(infilename):
	locus, physpos, genpos, DAF_selpop, delDAF,fst,xp_normed,ihs1_normed,delihh1_normed= [], [], [], [], [], [], [], [], []
	if not os.path.isfile(infilename):
		return locus, physpos, genpos, DAF_selpop, delDAF,fst,xp_normed,ihs1_normed,delihh1_normed
	openfile = open(infilename, 'r')
	openfile.readline()
	for line in openfile:
		entries = line.split()
		locus.append(int(entries[0]))
		physpos.append(int(entries[1]))
		genpos.append(float(entries[2]))
		DAF_selpop.append(float(entries[3]))
		delDAF.append(float(entries[4]))
		fst.append(float(entries[5]))
		xp_normed.append(float(entries[6]))
		ihs1_normed.append(float(entries[7]))
		delihh1_normed.append(float(entries[8]))
	openfile.close()
	return locus, physpos, genpos, DAF_selpop, delDAF,fst,xp_normed,ihs1_normed,delihh1_normed

def main():

	for irep in range(1,90):
		compositefile = "/idi/sabeti-scratch/jvitti/scores_composite2/default_default_101715_12pm/sel1/sel_0.10/rep" + str(irep) + "_1_vsneut.cms.out"

		physpos, genpos, ihs_normed, delihh_normed, xpehh_normed, fst, deldaf, cms_unnormed, cms_normed = read_cms_repfile(compositefile)

		#scoreindex = physpos.index(examinepos)
		#print('xp: ' + str(xpehh_normed[scoreindex]))
		#print('fst: ' + str(fst[scoreindex]))
		#print('deldaf: ' + str(deldaf[scoreindex]))
		
		#for pairfile in pairfiles:
		#	locus, pairphyspos, pairgenpos, DAF_selpop, delDAF,pairfst,xp_normed,ihs1_normed,delihh1_normed = read_cms_pairfile(pairfile)

		#	if examinepos in pairphyspos:
		#		compscoreindex = pairphyspos.index(examinepos)
				#print('xp: ' + str(xp_normed[compscoreindex]))
		#		print('fst: ' + str(pairfst[compscoreindex]))
				#print('deldaf: ' + str(delDAF[compscoreindex]))
		if len(cms_unnormed) > 0:
			cms_unnormed2 = [np.log(item) for item in cms_unnormed]
			print(min(cms_unnormed2))
		#fig, ax = plt.subplots()
		#ax.hist(cms_unnormed)
		#plt.savefig('/web/personal/vitti/test.png')

main()
