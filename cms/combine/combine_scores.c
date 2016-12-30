// last updated 12.30.16: new top-level program for compositing 		vitti@broadinstitute.org
// gcc -O0 -ggdb3 -lm -Wall -o combine_scores combine_scores.c cms_data.c

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"

/**********/
/***MAIN***/
/**********/

int main(int argc, char **argv) {
	popComp_data_multiple score_data;
	likes_data_multiple ihs_likes_data, nsl_likes_data, delihh_likes_data;
	likes_data_multiple xpehh_likes_data, fst_likes_data, deldaf_likes_data;	
	int nComparisons;
	//FILE *outf=NULL;
	//char outfilename[256]; 
	int ibin;
	int isnp, iComp;
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double compLikeRatio;
	char cms_param_filename[528];
    char ihs_master_likesfilename[256], nsl_master_likesfilename[256], delihh_master_likesfilename[256];
    char xpehh_master_likesfilename[256], fst_master_likesfilename[256], deldaf_master_likesfilename[256];    
	const int line_size = 15000000; 
	FILE *inf=NULL;

	if (argc <= 2) {
		fprintf(stderr, "Usage: ./combine_scores <cms_run_paramfile> <input_pair_file1> ...\n");
		exit(0);
	}
	nComparisons = argc - 2;
	//<cms_run_paramfile> has selpop and altpops and dem model and other options
	//	* some number of (pop-pairs), with a file for each (which in turn lists score files for each poppair to popComp)
	//	* likes: options (likesFreqs?) and top-level files
	//	* tosavefile? any other options? (maf filter? decompose? bf?)
	// 	* REGION? maybe this is separate/secondary? 
			//CMS_RUN_PARAMFILE: first six lines are six master_likesfiles that each have four lines: hit_hi, hit_mid, hit_lo, miss


	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "Preparing to load component scores...\n");
	get_popComp_data_multiple(&score_data, nComparisons, argc, argv); 
	fprintf(stderr, "\tloaded data object with %d snps and %d population comparisons.\n", score_data.nsnps, score_data.ncomp);
	//for (isnp = 0; isnp < score_data.nsnps; isnp++ ){fprintf(stderr, "%f\t", score_data.ihs_normed[1][isnp]);} // DEBUG

	/////////////////////////////
	// LOAD SCORE LIKELIHOODS ///
	/////////////////////////////
	fprintf(stderr, "Preparing to load score likelihoods...\n");
	sprintf(cms_param_filename, "%s", argv[1]);
	inf = fopen(cms_param_filename, "r"); 
	fgets(ihs_master_likesfilename, line_size, inf);
	strtok(ihs_master_likesfilename, "\n");
	fgets(nsl_master_likesfilename, line_size, inf);
	strtok(nsl_master_likesfilename, "\n");
	fgets(delihh_master_likesfilename, line_size, inf);
	strtok(delihh_master_likesfilename, "\n");
	fgets(xpehh_master_likesfilename, line_size, inf);
	strtok(xpehh_master_likesfilename, "\n");
	fgets(fst_master_likesfilename, line_size, inf);
	strtok(fst_master_likesfilename, "\n");
	fgets(deldaf_master_likesfilename, line_size, inf);
	strtok(deldaf_master_likesfilename, "\n");
	fclose(inf);

	get_likes_data_multiple(&ihs_likes_data, ihs_master_likesfilename); 
	get_likes_data_multiple(&nsl_likes_data, nsl_master_likesfilename); 
	get_likes_data_multiple(&delihh_likes_data, delihh_master_likesfilename); 
	get_likes_data_multiple(&xpehh_likes_data, xpehh_master_likesfilename); 
	get_likes_data_multiple(&fst_likes_data, fst_master_likesfilename); 
	get_likes_data_multiple(&deldaf_likes_data, deldaf_master_likesfilename); 
	for (ibin = 0; ibin < ihs_likes_data.nbins; ibin++){
		fprintf(stderr, "%f\t%f\t%f\t%f\t%f\t%f\n", ihs_likes_data.start_bin[ibin], ihs_likes_data.end_bin[ibin], ihs_likes_data.miss_probs[ibin], ihs_likes_data.hit_probs_hi[ibin], ihs_likes_data.hit_probs_mid[ibin], ihs_likes_data.hit_probs_low[ibin]);

	}

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	//outf = fopen(outfilename, "w");
	//assert(outf != NULL);
	//fprintf(stderr, "Preparing to write to: %s\n", outfilename);
	for (isnp = 0; isnp < score_data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < score_data.ncomp; iComp++){
			if (score_data.physpos[iComp][isnp] != 0){break;}
		} //advance to the first comparison for which we have any data
		thisihs = score_data.ihs_normed[iComp][isnp];
		thisihh = score_data.delihh_normed[iComp][isnp];
		thisnsl = score_data.nsl_normed[iComp][isnp];
		thisxpehh = compareXp(&score_data, isnp); //determine others by comparison. XP: take maximum.
		thisfst = compareFst(&score_data, isnp);	//Do I want to rewrite this to calculate LSBL? Should be easy since I'm passing the whole data object.
		thisdelDaf = comparedelDaf(&score_data, isnp);	//Similarly, it would be easy to rewrite this function to give us (daf - AVE(outgroup daf)). 
		//Will also need to redo likelihood tables similarly. For now, leave as-is.
		compLikeRatio = 1;
		thisdaf = score_data.daf_selpop[iComp][isnp];

		//thisdaf -- > determines which index we use for likes_data_multiple

		//compLikeRatio = delihh_bf * nsl_bf  * fst_bf * deldaf_bf * xpehh_bf; //* ihs_bf
		//fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf);
		//fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf, compLikeRatio);
	} // end isnp
	//fclose(outf);
	fprintf(stderr, "well, that's a wrap.\n");
	free_popComp_data_multiple(&score_data);
	free_likes_data_multiple(&ihs_likes_data);
	free_likes_data_multiple(&nsl_likes_data);
	free_likes_data_multiple(&delihh_likes_data);
	free_likes_data_multiple(&xpehh_likes_data);
	free_likes_data_multiple(&fst_likes_data);				
	free_likes_data_multiple(&deldaf_likes_data);		
	return 0;
} // end main