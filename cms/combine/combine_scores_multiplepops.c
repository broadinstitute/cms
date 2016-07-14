// for a set of likelihood tables, together with collated CMS comparison scores for a putative selected population vs. any number of outgroups, pulls and collates all component score statistics. 
// last updated: 07.11.16   vitti@broadinstitute.org

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"
#include "pop_comparison.h" 

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

	likes_data delihh_hit, delihh_miss, ihs_hit, ihs_miss, xpehh_hit, xpehh_miss, fst_hit, fst_miss, deldaf_hit, deldaf_miss;
	char delihh_hit_filename[256], delihh_miss_filename[256]; // LIKES
	char ihs_hit_filename[256], ihs_miss_filename[256];
	char xpehh_hit_filename[256], xpehh_miss_filename[256];
	char fst_hit_filename[256], fst_miss_filename[256];
	char deldaf_hit_filename[256], deldaf_miss_filename[256];

	if (argc < 13) {
		fprintf(stderr, "Usage: ./combine_cms_scores_multiplepops <outfilename> <delihh_hit_filename> <delihh_miss_filename> <ihs_hit_filename> <ihs_miss_filename> <xpehh_hit_filename> <xpehh_miss_filename> <fst_hit_filename> <fst_miss_filename> <deldaf_hit_filename> <deldaf_miss_filename> <popPair file 1> <popPair file 2...>\n");
		exit(0);
	}
	fprintf(stderr, "Preparing to load component scores...\n");
	get_popComp_data_multiple(&data, argc, argv[12:]); //argv[12:] or thereabouts, need to configure function. #possibly define likes in one file? 
	fprintf(stderr, "nsnps: %d\n", data.nsnps);

	strcpy(outfilename, argv[1]);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: ");
	fprintf(stderr, outfilename);
	fprintf(stderr, "\n");
	strcpy(delihh_hit_filename, argv[2]);
	strcpy(delihh_miss_filename, argv[3]);
	strcpy(ihs_hit_filename, argv[4]);
	strcpy(ihs_miss_filename, argv[5]);
	strcpy(xpehh_hit_filename, argv[6]);
	strcpy(xpehh_miss_filename, argv[7]);
	strcpy(fst_hit_filename, argv[8]);
	strcpy(fst_miss_filename, argv[9]);
	strcpy(deldaf_hit_filename, argv[10]);
	strcpy(deldaf_miss_filename, argv[11]);
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