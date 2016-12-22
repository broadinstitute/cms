## structures input to various functions of the script, power.py
## last updated 12.17.16

import argparse

def full_parser_power():
	parser=argparse.ArgumentParser(description="This script contains command-line utilities for calculating CMS 2.0 power from simulated data and significance for CMS scores from empirical data.")
	subparsers = parser.add_subparsers(help="sub-commands")
	
	#####################
	## SIMULATED DATA ###
	#####################
	run_neut_sims_parser = subparsers.add_parser('run_neut_sims', help='run additional neutral simulations to create test set')
	run_sel_sims_parser = subparsers.add_parser('run_sel_sims')
	for run_sims_parser in [run_neut_sims_parser, run_sel_sims_parser]:
		run_sims_parser.add_argument('--startrep', type=int, action="store", help='replicate number at which to start', default=1)
		run_sims_parser.add_argument('--paramfolder', type=str, action="store", help='location to read model parameters', default="/idi/sabeti-scratch/jvitti/params/")
	run_sel_sims_parser.add_argument('--trajdir', action="store")


	run_neut_repscores_parser = subparsers.add_parser('run_neut_repscores', help ='run per-replicate component scores')
	run_sel_repscores_parser = subparsers.add_parser('run_sel_repscores', help='run per-replicate component scores')
	for run_repscores_parser in [run_neut_repscores_parser, run_sel_repscores_parser]:
		run_repscores_parser.add_argument('--irep', type=int, action="store", help="replicate number", default=1)
		run_repscores_parser.add_argument('--tpedfolder', type=str, action="store", help="location of input tped data", default="/idi/sabeti-scratch/jvitti/clean/sims/")
		run_repscores_parser.add_argument('--simRecomFile', type=str, action="store", help="location of input recom file", default="/idi/sabeti-scratch/jvitti/params/test_recom.recom")

	run_norm_neut_repscores_parser = subparsers.add_parser('run_norm_neut_repscores')
	run_norm_neut_repscores_parser.add_argument('--edge', type=int, action="store", help="use interior of replicates; define per-end bp. (e.g. 1.5Mb -> 1Mb: 25000)")
	run_norm_neut_repscores_parser.add_argument('--chromlen', type=int, action="store", help="per bp (1.5mb = 1500000)")
	norm_from_binfile_parser = subparsers.add_parser('norm_from_binfile')
	sel_norm_from_binfile_parser = subparsers.add_parser('sel_norm_from_binfile')
	for norm_parser in [run_norm_neut_repscores_parser,norm_from_binfile_parser, sel_norm_from_binfile_parser]:
		norm_parser.add_argument('--score', action='store', default='')
		norm_parser.add_argument('--altpop', action='store', default='')

	for selbin_parser in [run_sel_sims_parser, sel_norm_from_binfile_parser]:
		selbin_parser.add_argument('--selbin', action="store")
	
	run_poppair_parser = subparsers.add_parser('run_poppair', help='collate component scores by population pair')
	run_poppair_parser.add_argument('--altpop', action='store')
	run_poppair_parser.add_argument('--scorebase', action='store', default="/n/regal/sabeti_lab/jvitti/clear/scores/")

	composite_sims_parser = subparsers.add_parser('composite_sims', help='calculate composite scores for simulations')
	normsims_parser = subparsers.add_parser('normsims', help="normalize simulated composite scores to neutral")

	#############################
	## LIKELIHOODS FROM SCORES ##
	#############################
	run_likes_parser = subparsers.add_parser('run_likes') # HPASE OUT
	run_likes_comp_parser = subparsers.add_parser('run_likes_comp')	 #PHASE OUT
	write_master_likes_parser = subparsers.add_parser('write_master_likes', help="write file for demographic scenario and population pointing to files for all score likelihoods")

	#####################
	## EMPIRICAL DATA ###
	#####################
	composite_emp_parser = subparsers.add_parser('composite_emp')	
	composite_emp_parser.add_argument('--basedir', type=str, action='store', default='/idi/sabeti-scratch/jvitti/scores_composite4_b/"')
	normemp_parser = subparsers.add_parser('normemp', help="normalize CMS scores to genome-wide") #norm emp REGIONS?

	#######################
	## VISUALIZE OUTPUT ###
	#######################
	repviz_parser = subparsers.add_parser('repviz', help="visualize component and combined scores")
	repviz_parser.add_argument('--vizRep', type=int, help="rep number", default=1, action="store")
	distviz_parser = subparsers.add_parser('distviz', help="visualize distribution of CMS scores")
	distviz_parser.add_argument('--oneFile', action="store", default=None, help="quick viz score dist from a single file")
	distviz_parser.add_argument('--simsDist', action="store_true", default=False, help="visualize from all simulated replicates")
	distviz_parser.add_argument('--empDist', action="store_true", default=False, help="visualize from empirical scores")
	distviz_parser.add_argument('--normed_cms', action="store_true", help="normalized values", default = False)



	#####################
	## QUANTIFY POWER ###
	#####################
	cdf_parser = subparsers.add_parser('cdf', help = 'plot cumulative density function of causal rank')
	fpr_parser = subparsers.add_parser('fpr', help='calculate false positive rate for CMS_gw based on neutral simulations')
	tpr_parser = subparsers.add_parser('tpr', help='calculate false positive rate for CMS_gw based on simulations with selection')
	roc_parser = subparsers.add_parser('roc', help="calculate receiving operator characteristic curve given false and true positive rates")
	roc_parser.add_argument('--plot_curve', action="store_true", default=False)
	roc_parser.add_argument('--find_opt', action="store_true", default=False)
	roc_parser.add_argument('--maxFPR', type=float, action="store", default=.001)
	#roc_parser.add_argument('--fprloc', type=str, action="store", default="/idi/sabeti-scratch/jvitti/cms2_power/fpr_4c/")
	#roc_parser.add_argument('--tprloc', type=str, action="store", default="/idi/sabeti-scratch/jvitti/cms2_power/tpr_4c/")

	#############################
	## EMPIRICAL SIGNIFICANCE ###
	#############################
	gw_regions_parser = subparsers.add_parser('gw_regions', help="pull designated significant regions from genome-wide normalized results")
	gw_regions_parser.add_argument('--geneFile')
	regionlog_parser = subparsers.add_parser('regionlog', help='write regions to excel sheet with gene overlap')
	regionlog_parser.add_argument('--gene_bedfile', help="name of file", type = str, action='store', default = "/idi/sabeti-scratch/jvitti/knownGenes_110116.txt")
	manhattan_parser = subparsers.add_parser('manhattan', help='generate manhattan plot of p-values of empirical results.')	
	manhattan_parser.add_argument('--zscores', action = 'store_true', help="plot -log10(p-values) estimated from neutral simulation")
	manhattan_parser.add_argument('--maxSkipVal', help="expedite plotting by ignoring anything obviously insignificant", default=-10e10)
	manhattan_parser.add_argument('--poolModels', help="experimental - combine neut distributions from multiple demographic models")
	extended_manhattan_parser = subparsers.add_parser('extended_manhattan', help = "generate per-chrom plots as one fig")
	extended_manhattan_parser.add_argument('--regionsfile', help="input file of regions designated as above threshhold")

	##################
	## SHARED ARGS ###
	##################

	for write_parser in [run_neut_sims_parser, run_neut_repscores_parser, run_norm_neut_repscores_parser, norm_from_binfile_parser, 
			run_poppair_parser, composite_sims_parser, run_sel_sims_parser, run_sel_repscores_parser, sel_norm_from_binfile_parser,
			normsims_parser, fpr_parser, tpr_parser, roc_parser]:
		write_parser.add_argument('--writedir', type =str, help='where to write output', default = "/idi/sabeti-scratch/jvitti/")
		write_parser.add_argument('--checkOverwrite', action="store_true", default=False)

	for makelikes_parser in [run_likes_parser, run_likes_comp_parser]:
		makelikes_parser.add_argument('numPerBin', type=int, action='store', default=1000, help="number of replicates in each bin")
		makelikes_parser.add_argument('selPos', type=int, action='store', default=500000, help="position of the causal allele in simulates")
		makelikes_parser.add_argument('numLikesBins', type=int, action='store', default=60, help='number of histogram bins')
	
	for uselikes_parser in [write_master_likes_parser, composite_sims_parser, composite_emp_parser]:
		uselikes_parser.add_argument('--likes_basedir', default="/idi/sabeti-scratch/jvitti/likes_111516_b/", help="location of likelihood tables ")

	for regions_parser in [fpr_parser, gw_regions_parser, regionlog_parser, tpr_parser]:
		regions_parser.add_argument('regionlen', type = int, action='store', help='length of region to query', default="100000")
		regions_parser.add_argument('thresshold', type = float, action='store', help='percentage of region to exceed cutoff', default="30")
		regions_parser.add_argument('cutoff', type = float, action='store', help='minimum significant value for region definition', default="3.0")
		regions_parser.add_argument('--saveLog', type =str, help="save results as text file", )


	for sim_parser in [normsims_parser, repviz_parser, distviz_parser, fpr_parser, tpr_parser, run_neut_repscores_parser,run_norm_neut_repscores_parser, norm_from_binfile_parser, composite_sims_parser, 
						run_neut_sims_parser, run_poppair_parser, run_sel_sims_parser, run_sel_repscores_parser, sel_norm_from_binfile_parser]:
		sim_parser.add_argument('--simpop', action='store', help='simulated population', default=1)

	for emp_parser in [normemp_parser, manhattan_parser, extended_manhattan_parser, gw_regions_parser, distviz_parser]:
		emp_parser.add_argument('--emppop', action='store', help='empirical population', default="YRI")

	for plot_parser in [repviz_parser, distviz_parser, manhattan_parser, extended_manhattan_parser, cdf_parser]:
		plot_parser.add_argument('--savefilename', action='store', help='path of file to save', default="/web/personal/vitti/test.png")

	for run_cms_parser in [run_neut_repscores_parser, run_norm_neut_repscores_parser, norm_from_binfile_parser, run_poppair_parser, composite_sims_parser, run_sel_repscores_parser, sel_norm_from_binfile_parser]:
		run_cms_parser.add_argument('--cmsdir', help='TEMPORARY, will become redundant with conda packaging', action = 'store', default= "/idi/sabeti-scratch/jvitti/cms/cms/")

	for sel_sim_parser in [normsims_parser]:
		sel_sim_parser.add_argument('--nrep_sel', type= int, action='store', default='500')
		sel_sim_parser.add_argument('--nrep_neut', type= int, action='store', default='1000')


	for sel_sim_parser in [normsims_parser, fpr_parser, tpr_parser]:
		sel_sim_parser.add_argument('--suffix', type= str, action='store', default='')



	for commonparser in [normsims_parser, repviz_parser, distviz_parser, fpr_parser, normemp_parser, gw_regions_parser, manhattan_parser, run_norm_neut_repscores_parser,norm_from_binfile_parser,
						regionlog_parser, cdf_parser, tpr_parser, extended_manhattan_parser, run_neut_sims_parser, run_neut_repscores_parser, composite_sims_parser, run_poppair_parser, 
						run_sel_sims_parser, run_sel_repscores_parser, sel_norm_from_binfile_parser]:
		commonparser.add_argument('--model', type=str, default="nulldefault")
		commonparser.add_argument('--nrep', type=int, default=1000)
		commonparser.add_argument('--likessuffix', action='store', help='neut or linked', default="neut")

	return parser
