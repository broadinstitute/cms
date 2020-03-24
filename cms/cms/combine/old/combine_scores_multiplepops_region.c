// for a set of likelihood tables, together with collated CMS comparison scores for a putative selected population vs. any number of outgroups, pulls and collates all component score statistics. 
// last updated: 11.01.16   vitti@broadinstitute.org

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
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double compLike, numerator, denominator;
	double prior;
	likes_data delihh_hit_low, delihh_hit_mid, delihh_hit_hi, delihh_miss_low, delihh_miss_mid, delihh_miss_hi;
	likes_data ihs_hit_low, ihs_hit_mid, ihs_hit_hi, ihs_miss_low, ihs_miss_mid, ihs_miss_hi;
	likes_data xpehh_hit_low, xpehh_hit_mid, xpehh_hit_hi, xpehh_miss_low, xpehh_miss_mid, xpehh_miss_hi;
	likes_data fst_hit_low, fst_hit_mid, fst_hit_hi, fst_miss_low, fst_miss_mid, fst_miss_hi;
	likes_data deldaf_hit_low, deldaf_hit_mid, deldaf_hit_hi, deldaf_miss_low, deldaf_miss_mid, deldaf_miss_hi;
	char delihh_hit_low_filename[256], delihh_miss_low_filename[256];
	char delihh_hit_mid_filename[256], delihh_miss_mid_filename[256]; 	
	char delihh_hit_hi_filename[256], delihh_miss_hi_filename[256]; 
	char ihs_hit_low_filename[256], ihs_miss_low_filename[256];
	char ihs_hit_mid_filename[256], ihs_miss_mid_filename[256];
	char ihs_hit_hi_filename[256], ihs_miss_hi_filename[256];
	char xpehh_hit_low_filename[256], xpehh_miss_low_filename[256];
	char xpehh_hit_mid_filename[256], xpehh_miss_mid_filename[256];
	char xpehh_hit_hi_filename[256], xpehh_miss_hi_filename[256];
	char fst_hit_low_filename[256], fst_miss_low_filename[256];
	char fst_hit_mid_filename[256], fst_miss_mid_filename[256];
	char fst_hit_hi_filename[256], fst_miss_hi_filename[256];
	char deldaf_hit_low_filename[256], deldaf_miss_low_filename[256];
	char deldaf_hit_mid_filename[256], deldaf_miss_mid_filename[256];
	char deldaf_hit_hi_filename[256], deldaf_miss_hi_filename[256];

	float delihh_hitprob, delihh_missprob;
	float ihs_hitprob, ihs_missprob;
	float xpehh_hitprob, xpehh_missprob;
	float fst_hitprob, fst_missprob;
	float deldaf_hitprob, deldaf_missprob;

	if (argc < 34) {
		fprintf(stderr, "Usage: ./combine_cms_scores_multiplepops_region <startbp> <endbp> <outfilename> <ihs_hit_filenames {low, mid, hi}> <ihs_miss_filenames {low, mid, hi}> <delihh_hit_filenames {low, mid, hi}> <delihh_miss_filenames {low, mid, hi}> <xpehh_hit_filenames {low, mid, hi}> <xpehh_miss_filenames {low, mid, hi}> <fst_hit_filenames {low, mid, hi}> <fst_miss_filenames {low, mid, hi}> <deldaf_hit_filenames {low, mid, hi}> <deldaf_miss_filenames {low, mid, hi}> <popPair file 1> <popPair file 2...>\n");
		exit(0);
	}
	fprintf(stderr, "Preparing to load component scores...\n");

	get_popComp_data_multiple_region(&data, argc, argv); 
	fprintf(stderr, "nsnps: %d\n", data.nsnps);

	strcpy(outfilename, argv[3]);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: %s\n", outfilename);
	strcpy(ihs_hit_low_filename, argv[2]);
	strcpy(ihs_hit_mid_filename, argv[3]);
	strcpy(ihs_hit_hi_filename, argv[4]);
	strcpy(ihs_miss_low_filename, argv[5]);
	strcpy(ihs_miss_mid_filename, argv[6]);
	strcpy(ihs_miss_hi_filename, argv[7]);
	strcpy(delihh_hit_low_filename, argv[8]);
	strcpy(delihh_hit_mid_filename, argv[9]);
	strcpy(delihh_hit_hi_filename, argv[10]);
	strcpy(delihh_miss_low_filename, argv[11]);
	strcpy(delihh_miss_mid_filename, argv[12]);
	strcpy(delihh_miss_hi_filename, argv[13]);
	strcpy(xpehh_hit_low_filename, argv[14]);
	strcpy(xpehh_hit_mid_filename, argv[15]);
	strcpy(xpehh_hit_hi_filename, argv[16]);
	strcpy(xpehh_miss_low_filename, argv[17]);
	strcpy(xpehh_miss_mid_filename, argv[18]);
	strcpy(xpehh_miss_hi_filename, argv[19]);
	strcpy(fst_hit_low_filename, argv[20]);
	strcpy(fst_hit_mid_filename, argv[21]);
	strcpy(fst_hit_hi_filename, argv[22]);
	strcpy(fst_miss_low_filename, argv[23]);
	strcpy(fst_miss_mid_filename, argv[24]);
	strcpy(fst_miss_hi_filename, argv[25]);
	strcpy(deldaf_hit_low_filename, argv[26]);
	strcpy(deldaf_hit_mid_filename, argv[27]);
	strcpy(deldaf_hit_hi_filename, argv[28]);
	strcpy(deldaf_miss_low_filename, argv[29]);
	strcpy(deldaf_miss_mid_filename, argv[30]);
	strcpy(deldaf_miss_hi_filename, argv[31]);
	//separate probability distributions based on observed derived allele frequency
	get_likes_data(&delihh_hit_low, delihh_hit_low_filename);
	get_likes_data(&delihh_miss_low, delihh_miss_low_filename);
	get_likes_data(&ihs_hit_low, ihs_hit_low_filename);
	get_likes_data(&ihs_miss_low, ihs_miss_low_filename);
	get_likes_data(&xpehh_hit_low, xpehh_hit_low_filename);
	get_likes_data(&xpehh_miss_low, xpehh_miss_low_filename);
	get_likes_data(&fst_hit_low, fst_hit_low_filename);
	get_likes_data(&fst_miss_low, fst_miss_low_filename);
	get_likes_data(&deldaf_hit_low, deldaf_hit_low_filename);
	get_likes_data(&deldaf_miss_low, deldaf_miss_low_filename);

	get_likes_data(&delihh_hit_mid, delihh_hit_mid_filename);
	get_likes_data(&delihh_miss_mid, delihh_miss_mid_filename);
	get_likes_data(&ihs_hit_mid, ihs_hit_mid_filename);
	get_likes_data(&ihs_miss_mid, ihs_miss_mid_filename);
	get_likes_data(&xpehh_hit_mid, xpehh_hit_mid_filename);
	get_likes_data(&xpehh_miss_mid, xpehh_miss_mid_filename);
	get_likes_data(&fst_hit_mid, fst_hit_mid_filename);
	get_likes_data(&fst_miss_mid, fst_miss_mid_filename);
	get_likes_data(&deldaf_hit_mid, deldaf_hit_mid_filename);
	get_likes_data(&deldaf_miss_mid, deldaf_miss_mid_filename);

	get_likes_data(&delihh_hit_hi, delihh_hit_hi_filename);
	get_likes_data(&delihh_miss_hi, delihh_miss_hi_filename);
	get_likes_data(&ihs_hit_hi, ihs_hit_hi_filename);
	get_likes_data(&ihs_miss_hi, ihs_miss_hi_filename);
	get_likes_data(&xpehh_hit_hi, xpehh_hit_hi_filename);
	get_likes_data(&xpehh_miss_hi, xpehh_miss_hi_filename);
	get_likes_data(&fst_hit_hi, fst_hit_hi_filename);
	get_likes_data(&fst_miss_hi, fst_miss_hi_filename);
	get_likes_data(&deldaf_hit_hi, deldaf_hit_hi_filename);
	get_likes_data(&deldaf_miss_hi, deldaf_miss_hi_filename);
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

		thisdaf = data.daf_selpop[iComp][isnp];

		if (thisdaf <= .35){
			delihh_hitprob = getProb(&delihh_hit_low, thisihh);
			ihs_hitprob = getProb(&ihs_hit_low, thisihs);
			fst_hitprob = getProb(&fst_hit_low, thisfst);
			deldaf_hitprob = getProb(&deldaf_hit_low, thisdelDaf);
			xpehh_hitprob = getProb(&xpehh_hit_low, thisxpehh);

			delihh_missprob = getProb(&delihh_miss_low, thisihh); 
			ihs_missprob = getProb(&ihs_miss_low, thisihs);
			fst_missprob = getProb(&fst_miss_low, thisfst);
			deldaf_missprob = getProb(&deldaf_miss_low, thisdelDaf);
			xpehh_missprob = getProb(&xpehh_miss_low, thisxpehh);
		}

		else if(thisdaf > .35 && thisdaf <= .65){
			delihh_hitprob = getProb(&delihh_hit_mid, thisihh);
			ihs_hitprob = getProb(&ihs_hit_mid, thisihs);
			fst_hitprob = getProb(&fst_hit_mid, thisfst);
			deldaf_hitprob = getProb(&deldaf_hit_mid, thisdelDaf);
			xpehh_hitprob = getProb(&xpehh_hit_mid, thisxpehh);

			delihh_missprob = getProb(&delihh_miss_mid, thisihh); 
			ihs_missprob = getProb(&ihs_miss_mid, thisihs);
			fst_missprob = getProb(&fst_miss_mid, thisfst);
			deldaf_missprob = getProb(&deldaf_miss_mid, thisdelDaf);
			xpehh_missprob = getProb(&xpehh_miss_mid, thisxpehh);
		}

		else{
			delihh_hitprob = getProb(&delihh_hit_hi, thisihh);
			ihs_hitprob = getProb(&ihs_hit_hi, thisihs);
			fst_hitprob = getProb(&fst_hit_hi, thisfst);
			deldaf_hitprob = getProb(&deldaf_hit_hi, thisdelDaf);
			xpehh_hitprob = getProb(&xpehh_hit_hi, thisxpehh);

			delihh_missprob = getProb(&delihh_miss_hi, thisihh); 
			ihs_missprob = getProb(&ihs_miss_hi, thisihs);
			fst_missprob = getProb(&fst_miss_hi, thisfst);
			deldaf_missprob = getProb(&deldaf_miss_hi, thisdelDaf);
			xpehh_missprob = getProb(&xpehh_miss_hi, thisxpehh);
		}

	compLike = 1;
	numerator = 1.;
	numerator *= delihh_hitprob * prior;
	numerator *= ihs_hitprob * prior;
	numerator *= fst_hitprob * prior;
	numerator *= deldaf_hitprob * prior;
	numerator *= xpehh_hitprob * prior;
	denominator = 1.;
	denominator *= ((delihh_missprob * (1-prior)) + (delihh_hitprob * prior));
	denominator *= ((ihs_missprob* (1-prior)) + (ihs_hitprob * prior));
	denominator *= ((fst_missprob* (1-prior)) + (fst_hitprob * prior));
	denominator *= ((deldaf_missprob* (1-prior)) + (deldaf_hitprob * prior));
	denominator *= ((xpehh_missprob* (1-prior)) + (xpehh_hitprob * prior));
		
	compLike = numerator / denominator;

	fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", data.physpos[iComp][isnp], data.genpos[iComp][isnp], thisihs, thisihh, thisxpehh, thisfst, thisdelDaf, compLike);
	}

	fclose(outf);
	free_likes_data(&delihh_hit_low);
	free_likes_data(&delihh_hit_mid);
	free_likes_data(&delihh_hit_hi);		
	free_likes_data(&delihh_miss_low);
	free_likes_data(&delihh_miss_mid);
	free_likes_data(&delihh_miss_hi);
	free_likes_data(&ihs_hit_low);
	free_likes_data(&ihs_hit_mid);
	free_likes_data(&ihs_hit_hi);
	free_likes_data(&ihs_miss_low);
	free_likes_data(&ihs_miss_mid);
	free_likes_data(&ihs_miss_hi);
	free_likes_data(&xpehh_hit_low);
	free_likes_data(&xpehh_hit_mid);
	free_likes_data(&xpehh_hit_hi);
	free_likes_data(&xpehh_miss_low);
	free_likes_data(&xpehh_miss_mid);
	free_likes_data(&xpehh_miss_hi);
	free_likes_data(&fst_hit_low);
	free_likes_data(&fst_hit_mid);
	free_likes_data(&fst_hit_hi);		
	free_likes_data(&fst_miss_low);
	free_likes_data(&fst_miss_mid);	
	free_likes_data(&fst_miss_hi);
	free_likes_data(&deldaf_hit_low);
	free_likes_data(&deldaf_hit_mid);
	free_likes_data(&deldaf_hit_hi);		
	free_likes_data(&deldaf_miss_low);
	free_likes_data(&deldaf_miss_mid);
	free_likes_data(&deldaf_miss_hi);		

	free_popComp_data_multiple(&data);
	return 0;
} // end main