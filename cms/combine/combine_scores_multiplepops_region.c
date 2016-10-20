// for a set of likelihood tables, together with collated CMS comparison scores for a putative selected population vs. any number of outgroups, pulls and collates all component score statistics. 
// last updated: 10.20.16   vitti@broadinstitute.org

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
	int iComp, isnp; //
	char outfilename[256]; 
	double thisihs, thisihh; // per-pop
	double thisfst, thisxpehh, thisdelDaf;
	double compLike, numerator, denominator;
	double prior;
	likes_data delihh_hit, delihh_miss, ihs_hit, ihs_miss, xpehh_hit, xpehh_miss, fst_hit, fst_miss, deldaf_hit, deldaf_miss;
	char delihh_hit_filename[256], delihh_miss_filename[256]; // LIKES
	char ihs_hit_filename[256], ihs_miss_filename[256];
	char xpehh_hit_filename[256], xpehh_miss_filename[256];
	char fst_hit_filename[256], fst_miss_filename[256];
	char deldaf_hit_filename[256], deldaf_miss_filename[256];


	float delihh_hitprob, delihh_missprob;
	float ihs_hitprob, ihs_missprob;
	float xpehh_hitprob, xpehh_missprob;
	float fst_hitprob, fst_missprob;
	float deldaf_hitprob, deldaf_missprob;

	if (argc < 15) {
		fprintf(stderr, "Usage: ./combine_cms_scores_multiplepops_region <startbp> <endbp> <outfilename> <ihs_hit_filename> <ihs_miss_filename> <delihh_hit_filename> <delihh_miss_filename> <xpehh_hit_filename> <xpehh_miss_filename> <fst_hit_filename> <fst_miss_filename> <deldaf_hit_filename> <deldaf_miss_filename> <popPair file 1> <popPair file 2...>\n");
		exit(0);
	}
	fprintf(stderr, "Preparing to load component scores...\n");

	get_popComp_data_multiple_region(&data, argc, argv); 
	fprintf(stderr, "nsnps: %d\n", data.nsnps);

	strcpy(outfilename, argv[3]);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: %s\n", outfilename);
	strcpy(ihs_hit_filename, argv[4]);
	strcpy(ihs_miss_filename, argv[5]);
	strcpy(delihh_hit_filename, argv[6]);
	strcpy(delihh_miss_filename, argv[7]);

	strcpy(xpehh_hit_filename, argv[8]);
	strcpy(xpehh_miss_filename, argv[9]);
	strcpy(fst_hit_filename, argv[10]);
	strcpy(fst_miss_filename, argv[11]);
	strcpy(deldaf_hit_filename, argv[12]);
	strcpy(deldaf_miss_filename, argv[13]);
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
	prior = 1. / data.nsnps;
	for (isnp = 0; isnp < data.nsnps; isnp++){

		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < data.ncomp; iComp++){
			if (data.physpos[iComp][isnp] != 0){break;}
		}
		thisihs = data.ihs_normed[iComp][isnp];
		thisihh = data.delihh_normed[iComp][isnp];

		thisxpehh = compareXp(&data, isnp);
		thisfst = compareFst(&data, isnp);
		thisdelDaf = comparedelDaf(&data, isnp);

	delihh_hitprob = getProb(&delihh_hit, thisihh);
	ihs_hitprob = getProb(&ihs_hit, thisihs);
	fst_hitprob = getProb(&fst_hit, thisfst);
	deldaf_hitprob = getProb(&deldaf_hit, thisdelDaf);
	xpehh_hitprob = getProb(&xpehh_hit, thisxpehh);

	delihh_missprob = getProb(&delihh_miss, thisihh); 
	ihs_missprob = getProb(&ihs_miss, thisihs);
	fst_missprob = getProb(&fst_miss, thisfst);
	deldaf_missprob = getProb(&deldaf_miss, thisdelDaf);
	xpehh_missprob = getProb(&xpehh_miss, thisxpehh);


	compLike = 1;
	numerator = 1.;
	numerator *= getProb(&delihh_hit, thisihh) * prior;
	numerator *= getProb(&ihs_hit, thisihs) * prior;
	numerator *= getProb(&fst_hit, thisfst)* prior;
	numerator *= getProb(&deldaf_hit, thisdelDaf) * prior;
	numerator *= getProb(&xpehh_hit, thisxpehh) * prior;
	denominator = 1.;
	denominator *= ((getProb(&delihh_miss, thisihh) * (1-prior)) + (getProb(&delihh_hit, thisihh) * prior));
	denominator *= ((getProb(&ihs_miss, thisihs)* (1-prior)) + (getProb(&ihs_hit, thisihh) * prior));
	denominator *= ((getProb(&fst_miss, thisfst)* (1-prior)) + (getProb(&fst_hit, thisfst) * prior));
	denominator *= ((getProb(&deldaf_miss, thisdelDaf)* (1-prior)) + (getProb(&deldaf_hit, thisdelDaf) * prior));
	denominator *= ((getProb(&xpehh_miss, thisxpehh)* (1-prior)) + (getProb(&xpehh_hit, thisxpehh) * prior));
		
	compLike = numerator / denominator;

	fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", data.physpos[iComp][isnp], data.genpos[iComp][isnp], thisihs, thisihh, thisxpehh, thisfst, thisdelDaf, compLike);
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