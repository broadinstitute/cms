// 03.05.17: quick duplicate program of combine_scores with per-component likelihood contribution given 
// gcc -O0 -ggdb3 -lm -Wall -o combine_scores_local_outlikes combine_scores_outlikes_local.c cms_data.c
// ./combine_scores_local_outlikes test_out.txt test_masterlikes_params.txt testpair1.txt testpair2.txt

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
	FILE *inf=NULL, *outf=NULL;	
	char *token, *running;
	char cms_param_filename[528], paramline[528], outfilename[256];
	char ihs_master_likesfilename[256], nsl_master_likesfilename[256], delihh_master_likesfilename[256];
    char xpehh_master_likesfilename[256], fst_master_likesfilename[256], deldaf_master_likesfilename[256];
	float delihh_hitprob, delihh_missprob, delihh_prob, delihh_minprob, delihh_maxprob; 
	float nsl_hitprob, nsl_missprob, nsl_prob, nsl_minprob, nsl_maxprob; //bayes factor
	float ihs_hitprob, ihs_missprob, ihs_prob, ihs_minprob, ihs_maxprob;
	float xpehh_hitprob, xpehh_missprob, xpehh_prob, xpehh_minprob, xpehh_maxprob;
	float fst_hitprob, fst_missprob, fst_prob, fst_minprob, fst_maxprob;
	float deldaf_hitprob, deldaf_missprob, deldaf_prob, deldaf_minprob, deldaf_maxprob;
	int isnp, iComp, itoken, thisPos, likesFreqIndex, nComparisons, maxPos, minPos;
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double compLike, minDaf;
	//int ibin;  //for debug
	int proceed; //Boolean used to log whether each SNP passes filter 0T 1F
	int takeIhs, takeDelihh, takeNsl, takeXpehh, takeFst, takeDeldaf; //Bools as above
	//char takeScoreString[6];
	double prior; // = 1/nSNP for region
	int nsnps_regional;

	if (argc <= 3) {
		fprintf(stderr, "Usage: ./combine_scores_local_outlikes <savefilename> <cms_run_paramfile> <input_pair_file1> ...\n");
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

	/////////////////////
	/// DEFINE REGION ///
	/////////////////////
	nsnps_regional = 0;
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
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = score_data.daf_selpop[iComp][isnp];
		if (thisdaf < minDaf){proceed=1;} 
		//if still a go...
		if(proceed == 0){nsnps_regional++;
		}
	} //end for isnp
	prior = 1. / nsnps_regional;

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	strcpy(outfilename, argv[1]);
	//fprintf(stderr, "Preparing to write to: %s\n", outfilename);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
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
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = score_data.daf_selpop[iComp][isnp];
		if (thisdaf < minDaf){proceed=1;} 
		//if still a go...
		if(proceed == 0){
			compLike = 1;
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

			delihh_minprob = getMinProb(&delihh_likes_data, likesFreqIndex, prior);
			nsl_minprob = getMinProb(&nsl_likes_data, likesFreqIndex, prior);
			ihs_minprob = getMinProb(&ihs_likes_data, likesFreqIndex, prior);	
			fst_minprob = getMinProb(&fst_likes_data, likesFreqIndex, prior);
			deldaf_minprob = getMinProb(&deldaf_likes_data, likesFreqIndex, prior);
			xpehh_minprob = getMinProb(&xpehh_likes_data, likesFreqIndex, prior);			
			
			delihh_maxprob = getMaxProb(&delihh_likes_data, likesFreqIndex, prior);
			nsl_maxprob = getMaxProb(&nsl_likes_data, likesFreqIndex, prior);
			ihs_maxprob = getMaxProb(&ihs_likes_data, likesFreqIndex, prior);	
			fst_maxprob = getMaxProb(&fst_likes_data, likesFreqIndex, prior);
			deldaf_maxprob = getMaxProb(&deldaf_likes_data, likesFreqIndex, prior);
			xpehh_maxprob = getMaxProb(&xpehh_likes_data, likesFreqIndex, prior);			
			
			///////////////////////////////////////////////////////
			//catch pseudocounts per SG/IS CMS 1.0 implementation// make this toggleable as well?
			///////////////////////////////////////////////////////
			if (delihh_missprob < 2e-10 && delihh_hitprob > 2e-10){delihh_prob = delihh_maxprob;}
			if (delihh_hitprob < 2e-10 && delihh_missprob > 2e-10){delihh_prob = delihh_minprob;}
			else{delihh_prob = (prior*delihh_hitprob) / ((prior*delihh_hitprob) + ((1.-prior)*delihh_missprob));} 

			if (nsl_missprob < 2e-10 && nsl_hitprob > 2e-10){nsl_prob = nsl_maxprob;}
			if (nsl_hitprob < 2e-10 && nsl_missprob > 2e-10){nsl_prob = nsl_minprob;}
			else{nsl_prob = (prior*nsl_hitprob) / ((prior*nsl_hitprob) + ((1.-prior)*nsl_missprob));} 

			if (ihs_missprob < 2e-10 && ihs_hitprob > 2e-10){ihs_prob = ihs_maxprob;}
			if (ihs_hitprob < 2e-10 && ihs_missprob > 2e-10){ihs_prob = ihs_minprob;}
			else {ihs_prob = (prior*ihs_hitprob) / ((prior*ihs_hitprob) + ((1.-prior)*ihs_missprob));} 

			if (fst_missprob < 2e-10 && fst_hitprob > 2e-10){fst_prob = fst_maxprob;}
			if (fst_hitprob < 2e-10 && fst_missprob > 2e-10){fst_prob = fst_minprob;}
			else{fst_prob = (prior*fst_hitprob) / ((prior*fst_hitprob) + ((1.-prior)*fst_missprob));} 

			if (deldaf_missprob < 2e-10 && deldaf_hitprob > 2e-10){deldaf_prob = deldaf_maxprob;}
			if (deldaf_hitprob < 2e-10 && deldaf_missprob > 2e-10){deldaf_prob = deldaf_minprob;}
			else{deldaf_prob = (prior*deldaf_hitprob) / ((prior*deldaf_hitprob) + ((1.-prior)*deldaf_missprob));} 

			if (xpehh_missprob < 2e-10 && xpehh_hitprob > 2e-10){xpehh_prob = xpehh_maxprob;}
			if (xpehh_hitprob < 2e-10 && xpehh_missprob > 2e-10){xpehh_prob = xpehh_minprob;}
			else{xpehh_prob = (prior*xpehh_hitprob) / ((prior*xpehh_hitprob) + ((1.-prior)*xpehh_missprob));} 
			
			/////////////////////
			/// GET CMS SCORE ///
			/////////////////////		
			if(takeIhs == 0){compLike *= ihs_prob;}// fprintf(stderr, "ihs\t");}
			if(takeDelihh == 0){compLike *= delihh_prob;}//fprintf(stderr, "delihh\t");}
			if(takeNsl == 0){compLike *= nsl_prob;}//;fprintf(stderr, "nsl\t");}
			if(takeFst == 0){compLike *= fst_prob;}//;fprintf(stderr, "fst\t");}
			if(takeDeldaf == 0){compLike *= deldaf_prob;}//;fprintf(stderr, "deldaf\t");}
			if(takeXpehh == 0){compLike *= xpehh_prob;}//;fprintf(stderr, "xp\n");}
		
			//DEBUG 
			/*fprintf(stderr, "ihs %f\t hit %e\tmiss %e\tbf %e\n", thisihs, ihs_hitprob, ihs_missprob, ihs_prob); //debug
			fprintf(stderr, "delihh %f\t hit %e\tmiss %e\tbf %e\n", thisihh, delihh_hitprob, delihh_missprob, delihh_prob); //debug
			fprintf(stderr, "fst %f\t hit %e\tmiss %e\tbf %e\n", thisfst, fst_hitprob, fst_missprob, fst_prob); //debug
			fprintf(stderr, "deldaf %f\t hit %e\tmiss %e\tbf %e\n", thisdelDaf, deldaf_hitprob, deldaf_missprob, deldaf_prob); //debug
			fprintf(stderr, "xp %f\t hit %e\tmiss %e\tbf %e\n", thisxpehh, xpehh_hitprob, xpehh_missprob, xpehh_prob); //debug
			fprintf(stderr, "clr: %e\n", compLike);
			fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf);*/
			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, ihs_prob, delihh_prob, nsl_prob, xpehh_prob, fst_prob, deldaf_prob, compLike);
		}//end if-a-go
	} // end isnp
	fclose(outf);
	fprintf(stderr, "Wrote to %s\n", outfilename);
	free_popComp_data_multiple(&score_data);
	free_likes_data_multiple(&ihs_likes_data);
	free_likes_data_multiple(&nsl_likes_data);
	free_likes_data_multiple(&delihh_likes_data);
	free_likes_data_multiple(&xpehh_likes_data);
	free_likes_data_multiple(&fst_likes_data);				
	free_likes_data_multiple(&deldaf_likes_data);		
	return 0;
} // end main