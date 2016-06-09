// {{COMPONENT SCORES + SCORE DISTRIBUTIONS FROM SIMULATED DATA --> COMPOSITE SCORES}}
// for a putative selpop and an arbitrary number of comparison populations, pulls and collates all component score statistics. 
// --> COMPILE: gcc -o combine_cms_scores_multiplepops -O0 -ggdb3 -lm -Wall combine_cms_scores_multiplepops.c cms_data.c
// last updated: 06.04.15   vitti@broadinstitute.org

//need to build in options to turn the knobs!, weighting of component scores, ways of handling comparisons
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"
 
 /* I SHOULD JUST MOVE THESE TO CMS_DATA BECAUSE THEY ARE SHARED NO NEED TO REDECLARE
/* does it matter if it's region or nonregion? have i differentiated in the functions
/* it might even be easier just to keep it here, s.t. we calculate CMS and CMS_GW separately...

/*********************/
/***DEFINE FUNCTIONS**/
/*********************/

float getProb(likes_data* data, double value);
float compareXp(popComp_data_multiple* data, int isnp);
float compareFst(popComp_data_multiple* data, int isnp);
float comparedelDaf(popComp_data_multiple* data, int isnp);

float getProb(likes_data* data, double value){
  int ibin;
  for (ibin = 0; ibin < data->nbins; ibin++){
    if (value >= data->start_bin[ibin] && value <= data->end_bin[ibin]){return data->probs[ibin];}
  }
  return 0;
}
float compareXp(popComp_data_multiple* data, int isnp){//currently: takes max val
  double xp;
  int iComp;
  xp = -100;
  for (iComp = 0; iComp < data->ncomp; iComp++){
    if (data->xp_normed[iComp][isnp] > xp){xp = data->xp_normed[iComp][isnp];}
  }
  return xp;
}
float compareFst(popComp_data_multiple* data, int isnp){//currently: takes max val
  double fst;
  int iComp;
  fst =  -100;
  //fprintf(stderr, "%d\n", data->ncomp);
  for (iComp = 0; iComp < data->ncomp; iComp++){
    if (data->fst[iComp][isnp] > fst){fst = data->fst[iComp][isnp];}
    //else{fprintf(stderr, "blimey what the eff: %f\n", data->fst[iComp][isnp]);}
  }

  return fst;
}
float comparedelDaf(popComp_data_multiple* data, int isnp){//currently: takes max val
  double deldaf;
  int iComp;
  deldaf = -100;
  for (iComp = 0; iComp < data->ncomp; iComp++){
    if (data->delDAF[iComp][isnp] > deldaf){deldaf = data->delDAF[iComp][isnp];}
  }
  return deldaf;
}

/**********/
/***MAIN***/
/**********/

int main(int argc, char **argv) {
	popComp_data_multiple data;
	FILE *outf=NULL;
	int nComparisons, iComp, isnp;
	char outfilename[256]; 
	double thisihs, thisihh; // per-pop
	double thisfst, thisxpehh, thisdelDaf;
	double compLikeRatio, numerator, denominator;

	if (argc < 4) {
		fprintf(stderr, "Usage: ./combine_cms_scores_multiplepops <chrom> <selpop> <poplikes> <otherpop1> <...> \n");
		exit(0);
	}
	fprintf(stderr, "Preparing to load component scores...\n");
	get_popComp_data_multiple(&data, argc, argv); 
	fprintf(stderr, "nsnps: %d\n", data.nsnps);

	strcpy(outfilename, "/idi/sabeti-scratch/jvitti/synth/cms_combined/chr");
	strcat(outfilename, argv[1]);
	strcat(outfilename, "_sel");
	strcat(outfilename, argv[2]);
	strcat(outfilename, "_vs_");
	strcat(outfilename, "TEST_addlikes_otherpops_etc.outcmsgw");
	//strcat(outfilename, "TEST_addlikes_otherpops_etc.outcmsgw");

	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: ");
	fprintf(stderr, outfilename);
	fprintf(stderr, "\n");

	/////////////////
	// LOAD LIKES ///
	/////////////////
	char delihh_hit_filename[256], delihh_miss_filename[256];
	char ihs_hit_filename[256], ihs_miss_filename[256];
	char xpehh_hit_filename[256], xpehh_miss_filename[256];
	char fst_hit_filename[256], fst_miss_filename[256];
	char deldaf_hit_filename[256], deldaf_miss_filename[256];
	char likeloc[64] = "/idi/sabeti-scratch/jvitti/synth/likes/";
	//"/idi/sabeti-scratch/jvitti/cms_venv/likes_p1/";

	likes_data delihh_hit, delihh_miss, ihs_hit, ihs_miss, xpehh_hit, xpehh_miss, fst_hit, fst_miss, deldaf_hit, deldaf_miss;

	fprintf(stderr, "loading likes tables from /idi/sabeti-scratch/jvitti/cms_venv/likes_p1/ (!!!)\n");
	fprintf(stderr, "loading neut dists from neut sims (CMS_GW)...\n");
	sprintf(delihh_hit_filename, likeloc);
	strcat(delihh_hit_filename, "delihh_default_112115_825am_");
		//"StdDiff_hits_");
	strcat(delihh_hit_filename, argv[3]);
	strcat(delihh_hit_filename, "_likes_causal.txt");

	sprintf(delihh_miss_filename, likeloc);
	strcat(delihh_miss_filename, "delihh_default_112115_825am_");
	//strcat(delihh_miss_filename, "StdDiff_miss_");
	strcat(delihh_miss_filename, argv[3]);
	strcat(delihh_miss_filename, "_likes_neut.txt");

	sprintf(ihs_hit_filename, likeloc);
		strcat(ihs_hit_filename, "ihs_default_112115_825am_");
	//strcat(ihs_hit_filename, "iHS_hits_");
	strcat(ihs_hit_filename, argv[3]);
	strcat(ihs_hit_filename, "_likes_causal.txt");

	sprintf(ihs_miss_filename, likeloc);
			strcat(ihs_miss_filename, "ihs_default_112115_825am_");
	//strcat(ihs_miss_filename, "iHS_miss_");
	strcat(ihs_miss_filename, argv[3]);
	strcat(ihs_miss_filename, "_likes_neut.txt");

	sprintf(xpehh_hit_filename, likeloc);
		strcat(xpehh_hit_filename, "xp_default_112115_825am_");
	//strcat(xpehh_hit_filename, "max_xpop_hits_");
	strcat(xpehh_hit_filename, "1_2");//argv[3]);
	strcat(xpehh_hit_filename, "_likes_causal.txt");

	sprintf(xpehh_miss_filename, likeloc);
			strcat(xpehh_miss_filename, "xp_default_112115_825am_");
	//strcat(xpehh_miss_filename, "max_xpop_miss_");
	strcat(xpehh_miss_filename, "1_2");//argv[3]);
	strcat(xpehh_miss_filename, "_likes_neut.txt");

	sprintf(deldaf_hit_filename, likeloc);
		strcat(deldaf_hit_filename, "deldaf_default_112115_825am_");
	//strcat(deldaf_hit_filename, "freqDiff_hits_");
	strcat(deldaf_hit_filename, "1_2");//argv[3]);
	strcat(deldaf_hit_filename, "_likes_causal.txt");

	sprintf(deldaf_miss_filename, likeloc);
		strcat(deldaf_miss_filename, "deldaf_default_112115_825am_");
	//strcat(deldaf_miss_filename, "freqDiff_miss_");
	strcat(deldaf_miss_filename, "1_2");//argv[3]);
	strcat(deldaf_miss_filename, "_likes_neut.txt");

	sprintf(fst_hit_filename, likeloc);
		strcat(fst_hit_filename, "fst_default_112115_825am_");
	//strcat(fst_hit_filename, "meanFst_hits_");
	strcat(fst_hit_filename, "1_2");//argv[3]);
	strcat(fst_hit_filename, "_likes_causal.txt");

	sprintf(fst_miss_filename, likeloc);
	strcat(fst_miss_filename, "fst_default_112115_825am_");
	//strcat(fst_miss_filename, "meanFst_miss_");
	strcat(fst_miss_filename, "1_2");//argv[3]);
	strcat(fst_miss_filename, "likes_neut.txt");
	
	get_likes_data(&delihh_hit, delihh_hit_filename);
	get_likes_data(&delihh_miss, delihh_miss_filename);
	get_likes_data(&ihs_hit, ihs_hit_filename);
	get_likes_data(&ihs_miss, ihs_miss_filename);
	get_likes_data(&xpehh_hit, xpehh_hit_filename);
	get_likes_data(&xpehh_miss, xpehh_miss_filename);
	get_likes_data(&fst_hit, fst_hit_filename);
	get_likes_data(&fst_miss, fst_miss_filename);
	get_likes_data(&deldaf_hit, deldaf_hit_filename);
	get_likes_data(&deldaf_miss, deldaf_miss_filename);

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////

	nComparisons = data.ncomp;
	for (isnp = 0; isnp < data.nsnps; isnp++){

		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////

		iComp = 0; //arbitrary; should confirm that values are consistent  //NVM, some will be zero if they arent seg in other pop.
		//fprintf(stderr, "%f\t%f\n", data.ihs_normed[iComp][isnp], data.ihs_normed[iComp + 1][isnp]);
		//assert(data.ihs_normed[iComp][isnp] == data.ihs_normed[iComp + 1][isnp]);
		thisihs = data.ihs_normed[iComp][isnp];
		thisihh = data.delihh_normed[iComp][isnp];

		thisxpehh = compareXp(&data, isnp);
		thisfst = compareFst(&data, isnp);
		thisdelDaf = comparedelDaf(&data, isnp);

	compLikeRatio = 1;
	numerator = 1;
	numerator *= getProb(&delihh_hit, thisihh);
	numerator *= getProb(&ihs_hit, thisihs);
	numerator *= getProb(&fst_hit, thisfst);
	numerator *= getProb(&deldaf_hit, thisdelDaf);
	numerator *= getProb(&xpehh_hit, thisxpehh);
	denominator = 1;
	denominator *= getProb(&delihh_miss, thisihh);
	denominator *= getProb(&ihs_miss, thisihs);
	denominator *= getProb(&fst_miss, thisfst);
	denominator *= getProb(&deldaf_miss, thisdelDaf);
	denominator *= getProb(&xpehh_miss, thisxpehh);
		

	compLikeRatio = numerator / denominator;

	fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", data.physpos[iComp][isnp], data.genpos[iComp][isnp], thisihs, thisihh, thisxpehh, thisfst, thisdelDaf, compLikeRatio);



	}

	fclose(outf);
	free_likes_data(&delihh_hit);
	free_likes_data(&delihh_miss);
	free_likes_data(&ihs_hit);
	free_likes_data(&ihs_miss);
	free_likes_data(&xpehh_hit);
	free_likes_data(&xpehh_miss);
	free_likes_data(&fst_hit);
	free_likes_data(&fst_miss);
	free_likes_data(&deldaf_hit);
	free_likes_data(&deldaf_miss);
	free_popComp_data_multiple(&data);
	return 0;
} // end main