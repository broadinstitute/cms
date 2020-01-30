## check distribution of p3 composite score Bayes Factors using output from combine_scores_multiplepops_decompose
## 11.05.06

infilename = "/idi/sabeti-scratch/jvitti/scores_composite2/gradient_101915_treebase_6_best/sel2/sel_0.70/rep23_2_vsneut.cms.out.decompose"
#"/idi/sabeti-scratch/jvitti/scores_composite2/default_default_101715_12pm/neut/rep459_1_vsneut.cms.out"
#"/idi/sabeti-scratch/jvitti/scores_composite2/default_default_101715_12pm/neut/rep459_4_vsneut.cms.out"
#"/idi/sabeti-scratch/jvitti/scores_composite2/default_default_101715_12pm/sel4/sel_0.60/rep459_4_vsneut.cms.out"
#"/idi/sabeti-scratch/jvitti/scores_composite2/nulldefault/sel1/sel_0.30/rep459_1_vsneut.cms.out"
#filebase = "/idi/sabeti-scratch/jvitti/synth/cms_composite/"
nbins = 500
import numpy as np
import matplotlib as mp
mp.use('agg')
import matplotlib.pyplot as plt

def get_log10_bf(scorearray):
	#print(scorearray)
	scorearray = [float(item) for item in scorearray if not np.isnan(float(item))]
	log10_array = [np.log10(item) for item in scorearray]
	return log10_array
def readvals_bayesfactors(filename):
	physpos, delihh_bfs, ihs_bfs, fst_bfs, deldaf_bfs, xpehh_bfs = [], [], [], [], [], []
	openfile = open(filename, 'r')
	for line in openfile:
		entries = line.split()
		pos, delihh_bf, ihs_bf, fst_bf, deldaf_bf, xpehh_bf = entries
		physpos.append(pos)
		delihh_bfs.append(delihh_bf)
		ihs_bfs.append(ihs_bf)
		fst_bfs.append(fst_bf)
		deldaf_bfs.append(deldaf_bf)
		xpehh_bfs.append(xpehh_bf)
	openfile.close()
	delihh_bfs = get_log10_bf(delihh_bfs)
	ihs_bfs = get_log10_bf(ihs_bfs)
	fst_bfs = get_log10_bf(fst_bfs)
	deldaf_bfs = get_log10_bf(deldaf_bfs)
	xpehh_bfs = get_log10_bf(xpehh_bfs)

	return delihh_bfs, ihs_bfs, fst_bfs, deldaf_bfs, xpehh_bfs

def main():
	

	#print(infilename)
	delihh_bfs, ihs_bfs, fst_bfs, deldaf_bfs, xpehh_bfs = readvals_bayesfactors(infilename)
	#for scorearray in [delihh_bfs, ihs_bfs, fst_bfs, deldaf_bfs, xpehh_bfs]:
		#print(scorearray[2])

	f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1)
	ax1.hist(delihh_bfs, bins = nbins)
	#ax1.set_xlim([-15, 15])
	ax2.hist(ihs_bfs, bins = nbins)
	ax3.hist(fst_bfs, bins = nbins)
	ax4.hist(deldaf_bfs, bins = nbins)
	ax5.hist(xpehh_bfs, bins = nbins)

	savefilename = "/web/personal/vitti/" + "test_decompose.png" #examinefile.split('.')[0] + ".png"

	plt.savefig(savefilename)
	print('plotted to ' + savefilename)

main()
