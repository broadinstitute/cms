// 	last updated 06.12.2017: furnishes a composite scores as a summary statistic of Bayes factors: P( score | sel v unsel)	
//	vitti@broadinstitute.org
// 		CMS_RUN_PARAMFILE: first six lines are six master_likesfiles that each have four lines: hit_hi, hit_mid, hit_lo, miss; 
// 		optional next line: (minPos, maxPos, minDaf, writeLikes); optional next line 0T 1F 6x for ihs ihh nsl fst deldaf xpehh

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
	const int line_size = 15000000; 
	popComp_data_multiple score_data;
	likes_data_multiple ihs_likes_data, nsl_likes_data, delihh_likes_data;
	likes_data_multiple xpehh_likes_data, fst_likes_data, deldaf_likes_data;	
	FILE *inf=NULL, *outf=NULL, *outf2=NULL;	
	char *token, *running;
	char cms_param_filename[528], paramline[528], outfilename[256], outfilename_likes[256];
	char ihs_master_likesfilename[256], nsl_master_likesfilename[256], delihh_master_likesfilename[256];
    char xpehh_master_likesfilename[256], fst_master_likesfilename[256], deldaf_master_likesfilename[256];
	float delihh_hitprob, delihh_missprob, delihh_bf, delihh_minbf, delihh_maxbf; 
	float nsl_hitprob, nsl_missprob, nsl_bf, nsl_minbf, nsl_maxbf; //bayes factor
	float ihs_hitprob, ihs_missprob, ihs_bf, ihs_minbf, ihs_maxbf;
	float xpehh_hitprob, xpehh_missprob, xpehh_bf, xpehh_minbf, xpehh_maxbf;
	float fst_hitprob, fst_missprob, fst_bf, fst_minbf, fst_maxbf;
	float deldaf_hitprob, deldaf_missprob, deldaf_bf, deldaf_minbf, deldaf_maxbf;
	int isnp, iComp, itoken, thisPos, likesFreqIndex, nComparisons, maxPos, minPos;
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double compLikeRatio, minDaf;
	//int ibin;  //for debug
	int proceed; //Boolean used to log whether each SNP passes filter 0T 1F
	int takeIhs, takeDelihh, takeNsl, takeXpehh, takeFst, takeDeldaf; //Bools as above
	int writeLikes;
	//char takeScoreString[6];

	if (argc <= 3) {
		fprintf(stderr, "Usage: ./combine_scores <savefilename> <cms_run_paramfile> <input_pair_file1> ...\n");
		exit(0);
	}
	nComparisons = argc - 3;
	
	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "Preparing to load component scores...\n");
	get_popComp_data_multiple(&score_data, nComparisons, argc, argv); 
	fprintf(stderr, "\tloaded data object with %d snps and %d population comparisons.\n", score_data.nsnps, score_data.ncomp);
	//for (isnp = 0; isnp < score_data.nsnps; isnp++ ){fprintf(stderr, "%f\t", score_data.ihs_normed[1][isnp]);} // DEBUG

	////////////////////////////////////////
	// LOAD SCORE LIKELIHOODS (DEM MODEL) //
	// AND OTHER RUN PARAMETERS ////////////
	////////////////////////////////////////
	fprintf(stderr, "Preparing to load score likelihoods...\n");
	sprintf(cms_param_filename, "%s", argv[2]);
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
	//set defaults
	minPos = -1;		
	maxPos = 2147483647;
	minDaf = 0;
	takeIhs = takeDelihh = takeNsl = takeXpehh = takeFst = takeDeldaf = 0; //all T by default
	//if additional lines is included, parse 
	if (fgets(paramline, line_size, inf) != NULL){
		for (running = paramline, itoken=0; (token = strsep(&running, " \t")) != NULL; itoken++){
			if (itoken == 0) {minPos = atoi(token);}
			else if (itoken == 1){maxPos = atoi(token);}
			else if (itoken == 2){minDaf = atof(token);}
			else if (itoken == 3){writeLikes = atoi(token);}			
		} // end for running
	}  //end if fgets paramline
	if (fgets(paramline, line_size, inf) != NULL){
		for (running = paramline, itoken=0; (token = strsep(&running, " \t")) != NULL; itoken++){
			if (itoken == 0) {takeIhs = atoi(token);}
			else if (itoken == 1){takeDelihh = atoi(token);}
			else if (itoken == 2){takeNsl = atoi(token);}			
			else if (itoken == 3){takeFst = atoi(token);}
			else if (itoken == 4){takeDeldaf = atoi(token);}			
			else if (itoken == 5){takeXpehh = atoi(token);}		
		} // end for running
	}  //end if fgets paramline
	fclose(inf);
	fprintf(stderr, "loaded parameters: minPos %d maxPos %d minDaf %f\n", minPos, maxPos, minDaf);		
	get_likes_data_multiple(&ihs_likes_data, ihs_master_likesfilename); 
	get_likes_data_multiple(&nsl_likes_data, nsl_master_likesfilename); 
	get_likes_data_multiple(&delihh_likes_data, delihh_master_likesfilename); 
	get_likes_data_multiple(&xpehh_likes_data, xpehh_master_likesfilename); 
	get_likes_data_multiple(&fst_likes_data, fst_master_likesfilename); 
	get_likes_data_multiple(&deldaf_likes_data, deldaf_master_likesfilename); 
	//for (ibin = 0; ibin < ihs_likes_data.nbins; ibin++){fprintf(stderr, "%f\t%f\t%f\t%f\t%f\t%f\n", ihs_likes_data.start_bin[ibin], ihs_likes_data.end_bin[ibin], ihs_likes_data.miss_probs[ibin], ihs_likes_data.hit_probs_hi[ibin], ihs_likes_data.hit_probs_mid[ibin], ihs_likes_data.hit_probs_low[ibin]);} // DEBUG

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	strcpy(outfilename, argv[1]);
	//fprintf(stderr, "Preparing to write to: %s\n", outfilename);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(outf, "physPos\tgenPos\tpopDAF\tnormed_iHS\tnormed_deliHH\tnormed_nsl\tnormed_xp-ehh\tfst\tdelDAF\tcompLikeRatio_CMS");

	if (writeLikes == 0){
		strcpy(outfilename_likes, argv[1]);
		strcat(outfilename_likes, ".likes");
		outf2 = fopen(outfilename_likes, "w");
		assert(outf2 != NULL);
		fprintf2(outf, "physPos\tgenPos\tpopDAF\tlike_iHS\tlike_deliHH\tlike_nsl\tlike_xp-ehh\tlikefst\tlikedelDAF\tcompLikeRatio_CMS");	
	} //end if write likes
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
		thisxpehh = compareXp(&score_data, isnp);
		thisfst = compareFst_PBS(&score_data, isnp);
		thisdelDaf = comparedelDaf_outgroup_ave(&score_data, isnp);
		if (thisihs < 0){thisihs*=-1.;}				//ABS VAL
		if (thisnsl < 0){thisnsl*=-1.;} 			//ABS VAL				
		if (thisdelDaf < 0){thisdelDaf*=-1.;} 		//ABS VAL -- FOLDED LIKELIHOOD DIST FOR CMS_GW
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = score_data.daf_selpop[iComp][isnp];
		if (thisdaf < minDaf){proceed=1;} 
		compLikeRatio = 0;
		//if still a go...
		if(proceed == 0){
			compLikeRatio = 1;
			/////////////////////////////////////
			//LIKESFREQS (current default)
			//thisdaf -- > determines which index we use for likes_data_multiple
			if (thisdaf <= .35){likesFreqIndex = 0;}
			else if(thisdaf > .35 && thisdaf <= .65){likesFreqIndex =1;}
			else{likesFreqIndex = 2;}

			delihh_hitprob = getHitProb(&delihh_likes_data, likesFreqIndex, thisihh);
			nsl_hitprob = getHitProb(&nsl_likes_data, likesFreqIndex, thisnsl);			
			ihs_hitprob = getHitProb(&ihs_likes_data, likesFreqIndex, thisihs);
			fst_hitprob = getHitProb(&fst_likes_data, likesFreqIndex, thisfst);
			deldaf_hitprob = getHitProb(&deldaf_likes_data, likesFreqIndex, thisdelDaf);
			xpehh_hitprob = getHitProb(&xpehh_likes_data, likesFreqIndex, thisxpehh);
		
			delihh_missprob = getMissProb(&delihh_likes_data, thisihh);
			nsl_missprob = getMissProb(&nsl_likes_data, thisnsl);			
			ihs_missprob = getMissProb(&ihs_likes_data, thisihs);
			fst_missprob = getMissProb(&fst_likes_data, thisfst);
			deldaf_missprob = getMissProb(&deldaf_likes_data, thisdelDaf);
			xpehh_missprob = getMissProb(&xpehh_likes_data, thisxpehh);

			delihh_minbf = getMinBf(&delihh_likes_data, likesFreqIndex);
			nsl_minbf = getMinBf(&nsl_likes_data, likesFreqIndex);
			ihs_minbf = getMinBf(&ihs_likes_data, likesFreqIndex);	
			fst_minbf = getMinBf(&fst_likes_data, likesFreqIndex);
			deldaf_minbf = getMinBf(&deldaf_likes_data, likesFreqIndex);
			xpehh_minbf = getMinBf(&xpehh_likes_data, likesFreqIndex);			
			
			delihh_maxbf = getMaxBf(&delihh_likes_data, likesFreqIndex);
			nsl_maxbf = getMaxBf(&nsl_likes_data, likesFreqIndex);
			ihs_maxbf = getMaxBf(&ihs_likes_data, likesFreqIndex);	
			fst_maxbf = getMaxBf(&fst_likes_data, likesFreqIndex);
			deldaf_maxbf = getMaxBf(&deldaf_likes_data, likesFreqIndex);
			xpehh_maxbf = getMaxBf(&xpehh_likes_data, likesFreqIndex);			
			
			///////////////////////////////////////////////////////
			//catch pseudocounts per SG/IS CMS 1.0 implementation// make this toggleable as well?
			///////////////////////////////////////////////////////
			delihh_bf = 0;
			if (delihh_missprob > 2e-10 && delihh_hitprob > 2e-10){delihh_bf = delihh_hitprob / delihh_missprob;}
			if (delihh_missprob < 2e-10 && delihh_hitprob > 2e-10){delihh_bf = delihh_maxbf;}
			if (delihh_hitprob < 2e-10 && delihh_missprob > 2e-10){delihh_bf = delihh_minbf;}
			if (delihh_hitprob < 2e-10 && delihh_missprob < 2e-10){delihh_bf = 1;} // no data
			//HMM, SHOULD THIS EVEN HAPPEN THOUGH?? IT SHOULD GET BINNED...

			nsl_bf = 0;
			if (nsl_missprob > 2e-10 && nsl_hitprob > 2e-10){nsl_bf = nsl_hitprob / nsl_missprob;}
			if (nsl_missprob < 2e-10 && nsl_hitprob > 2e-10){nsl_bf = nsl_maxbf;}
			if (nsl_hitprob < 2e-10 && nsl_missprob > 2e-10){nsl_bf = nsl_minbf;}

			ihs_bf = 0;
			if (ihs_missprob > 2e-10 && ihs_hitprob > 2e-10){ihs_bf = ihs_hitprob / ihs_missprob;}
			if (ihs_missprob < 2e-10 && ihs_hitprob > 2e-10){ihs_bf = ihs_maxbf;}
			if (ihs_hitprob < 2e-10 && ihs_missprob > 2e-10){ihs_bf = ihs_minbf;}

			fst_bf = 0;
			if (fst_missprob > 2e-10 && fst_hitprob > 2e-10){fst_bf = fst_hitprob / fst_missprob;}
			if (fst_missprob < 2e-10 && fst_hitprob > 2e-10){fst_bf = fst_maxbf;}
			if (fst_hitprob < 2e-10 && fst_missprob > 2e-10){fst_bf = fst_minbf;}

			deldaf_bf = 0;
			if (deldaf_missprob > 2e-10 && deldaf_hitprob > 2e-10){deldaf_bf = deldaf_hitprob / deldaf_missprob;}
			if (deldaf_missprob < 2e-10 && deldaf_hitprob > 2e-10){deldaf_bf = deldaf_maxbf;}
			if (deldaf_hitprob < 2e-10 && deldaf_missprob > 2e-10){deldaf_bf = deldaf_minbf;}
			
			xpehh_bf = 0;
			if (xpehh_missprob > 2e-10 && xpehh_hitprob > 2e-10){xpehh_bf = xpehh_hitprob / xpehh_missprob;}			
			if (xpehh_missprob < 2e-10 && xpehh_hitprob > 2e-10){xpehh_bf = xpehh_maxbf;}
			if (xpehh_hitprob < 2e-10 && xpehh_missprob > 2e-10){xpehh_bf = xpehh_minbf;}

			
			/////////////////////
			/// GET CMS SCORE ///
			/////////////////////		
			if(takeIhs == 0){compLikeRatio *= ihs_bf;}			//;fprintf(stderr, "ihs\t");}
			if(takeDelihh == 0){compLikeRatio *= delihh_bf;}	//;fprintf(stderr, "delihh\t");}
			if(takeNsl == 0){compLikeRatio *= nsl_bf;}			//;fprintf(stderr, "nsl\t");}
			if(takeFst == 0){compLikeRatio *= fst_bf;}			//;fprintf(stderr, "fst\t");}
			if(takeDeldaf == 0){compLikeRatio *= deldaf_bf;}	//;fprintf(stderr, "deldaf\t");}
			if(takeXpehh == 0){compLikeRatio *= xpehh_bf;}		//;fprintf(stderr, "xp\n");}

			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf, compLikeRatio);
			if (writeLikes == 0){fprintf(outf2, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, ihs_bf, delihh_bf, nsl_bf, xpehh_bf, fst_bf, deldaf_bf, compLikeRatio);} //end if write likes
		
			/*
			//DEBUG 
			fprintf(stderr, "ihs %f\t hit %e\tmiss %e\tbf %e\n", thisihs, ihs_hitprob, ihs_missprob, ihs_bf); //debug
			fprintf(stderr, "delihh %f\t hit %e\tmiss %e\tbf %e\n", thisihh, delihh_hitprob, delihh_missprob, delihh_bf); //debug
			fprintf(stderr, "fst %f\t hit %e\tmiss %e\tbf %e\n", thisfst, fst_hitprob, fst_missprob, fst_bf); //debug
			fprintf(stderr, "deldaf %f\t hit %e\tmiss %e\tbf %e\n", thisdelDaf, deldaf_hitprob, deldaf_missprob, deldaf_bf); //debug
			fprintf(stderr, "xp %f\t hit %e\tmiss %e\tbf %e\n", thisxpehh, xpehh_hitprob, xpehh_missprob, xpehh_bf); //debug
			fprintf(stderr, "clr: %e\n", compLikeRatio);
			fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf);
			*/
		}//end if-a-go
	} // end isnp
	fclose(outf);
	fprintf(stderr, "Wrote CMS scores to %s\n", outfilename);
	if (writeLikes == 0){fprintf(stderr, "Wrote score decomposition to %s\n", outfilename_likes);}
	free_popComp_data_multiple(&score_data);
	free_likes_data_multiple(&ihs_likes_data);
	free_likes_data_multiple(&nsl_likes_data);
	free_likes_data_multiple(&delihh_likes_data);
	free_likes_data_multiple(&xpehh_likes_data);
	free_likes_data_multiple(&fst_likes_data);				
	free_likes_data_multiple(&deldaf_likes_data);		
	return 0;
} // end main