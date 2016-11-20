// for a set of likelihood tables, together with collated CMS comparison scores for a putative selected population vs. any number of outgroups, pulls and collates all component score statistics. 
// last updated: 11.15.16   vitti@broadinstitute.org //phasing out.

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
	double compLikeRatio;//, numerator, denominator;

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

	float delihh_hitprob, delihh_missprob, delihh_bf, delihh_minbf; //bayes factor
	float ihs_hitprob, ihs_missprob, ihs_bf, ihs_minbf;
	float xpehh_hitprob, xpehh_missprob, xpehh_bf, xpehh_minbf;
	float fst_hitprob, fst_missprob, fst_bf, fst_minbf;
	float deldaf_hitprob, deldaf_missprob, deldaf_bf, deldaf_minbf;

	if (argc < 32) {
		fprintf(stderr, "Usage: ./combine_cms_scores_multiplepops <outfilename> <ihs_hit_filenames {low, mid, hi}> <ihs_miss_filenames {low, mid, hi}> <delihh_hit_filenames {low, mid, hi}> <delihh_miss_filenames {low, mid, hi}> <xpehh_hit_filenames {low, mid, hi}> <xpehh_miss_filenames {low, mid, hi}> <fst_hit_filenames {low, mid, hi}> <fst_miss_filenames {low, mid, hi}> <deldaf_hit_filenames {low, mid, hi}> <deldaf_miss_filenames {low, mid, hi}> <popPair file 1> <popPair file 2...>\n");
		exit(0);
	}
	fprintf(stderr, "Preparing to load component scores...\n");

	get_popComp_data_multiple(&data, argc, argv); 
	fprintf(stderr, "nsnps: %d\n", data.nsnps);

	strcpy(outfilename, argv[1]);
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

		compLikeRatio = 1;

		thisdaf = data.daf_selpop[iComp][isnp];
		//fprintf(stderr, "%d\t%d\t%f\t", data.physpos[iComp][isnp], isnp, thisdaf); //debug

		if (thisdaf <= .35){
			//fprintf(stderr, "lowfreq\n"); //debug
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

			delihh_minbf = getMinBf(&delihh_miss_low, &delihh_hit_low);
			ihs_minbf = getMinBf(&ihs_miss_low, &ihs_hit_low);
			fst_minbf = getMinBf(&fst_miss_low, &fst_hit_low);
			deldaf_minbf = getMinBf(&deldaf_miss_low, &deldaf_hit_low);
			xpehh_minbf = getMinBf(&xpehh_miss_low, &xpehh_hit_low);

		}

		else if(thisdaf > .35 && thisdaf <= .65){
			//fprintf(stderr, "midfreq\n"); //debug
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
		
			delihh_minbf = getMinBf(&delihh_miss_mid, &delihh_hit_mid);
			ihs_minbf = getMinBf(&ihs_miss_mid, &ihs_hit_mid);
			fst_minbf = getMinBf(&fst_miss_mid, &fst_hit_mid);
			deldaf_minbf = getMinBf(&deldaf_miss_mid, &deldaf_hit_mid);
			xpehh_minbf = getMinBf(&xpehh_miss_mid, &xpehh_hit_mid);
		}

		else{
			//fprintf(stderr, "hifreq\n"); //debug
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

			delihh_minbf = getMinBf(&delihh_miss_hi, &delihh_hit_hi);
			ihs_minbf = getMinBf(&ihs_miss_hi, &ihs_hit_hi);
			fst_minbf = getMinBf(&fst_miss_hi, &fst_hit_hi);
			deldaf_minbf = getMinBf(&deldaf_miss_hi, &deldaf_hit_hi);
			xpehh_minbf = getMinBf(&xpehh_miss_hi, &xpehh_hit_hi);

		}

		//catch pseudocounts per SG/IS CMS 1.0 implementation
		if (delihh_hitprob < 2e-10){ 
			delihh_bf = delihh_minbf;
		}
		else{
			delihh_bf = delihh_hitprob / delihh_missprob;
		}
		if (ihs_hitprob < 2e-10){
			ihs_bf = ihs_minbf;
		}
		else{
			ihs_bf = ihs_hitprob / ihs_missprob;
		}
		if (fst_hitprob < 2e-10){
			fst_bf = fst_minbf;
		}
		else{
			fst_bf = fst_hitprob / fst_missprob;
		}
		if (deldaf_hitprob < 2e-10){
			deldaf_bf = deldaf_minbf;
		}
		else{
			deldaf_bf = deldaf_hitprob / deldaf_missprob;			
		}
		if (xpehh_hitprob < 2e-10){
			xpehh_bf = xpehh_minbf;
		}
		else{
			xpehh_bf = xpehh_hitprob / xpehh_missprob;
		}

		compLikeRatio = delihh_bf * ihs_bf * fst_bf * deldaf_bf * xpehh_bf;

		//fprintf(stderr, "ihs %f\t hit %e\tmiss %e\tbf %e\n", thisihs, ihs_hitprob, ihs_missprob, ihs_bf); //debug
		//fprintf(stderr, "delihh %f\t hit %e\tmiss %e\tbf %e\n", thisihh, delihh_hitprob, delihh_missprob, delihh_bf); //debug
		//fprintf(stderr, "fst %f\t hit %e\tmiss %e\tbf %e\n", thisfst, fst_hitprob, fst_missprob, fst_bf); //debug
		//fprintf(stderr, "deldaf %f\t hit %e\tmiss %e\tbf %e\n", thisdelDaf, deldaf_hitprob, deldaf_missprob, deldaf_bf); //debug
		//fprintf(stderr, "xp %f\t hit %e\tmiss %e\tbf %e\n", thisxpehh, xpehh_hitprob, xpehh_missprob, xpehh_bf); //debug
		//fprintf(stderr, "clr: %e\n", compLikeRatio);

		fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", data.physpos[iComp][isnp], data.genpos[iComp][isnp], thisihs, thisihh, thisxpehh, thisfst, thisdelDaf, compLikeRatio);
	} // end isnp

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