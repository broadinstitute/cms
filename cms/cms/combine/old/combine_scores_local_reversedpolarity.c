// 	calculate within-region Composite of Multiple Signals (CMS) statistic - used alternate encoding (0der)
//	last updated 09.19.2017 	vitti@broadinstitute.org

//(cms-env3) [jvitti@holy2a18308 combine]$ gcc -c combine_scores_local_reversedpolarity.c
//(cms-env3) [jvitti@holy2a18308 combine]$ gcc -O0 -ggdb3 -lm -Wall -o combine_scores_local_reversedpolarity combine_scores_local_reversedpolarity.o cms_data.c

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
	float delihh_hitprob, delihh_missprob, nsl_hitprob, nsl_missprob, ihs_hitprob, ihs_missprob;
	float xpehh_hitprob, xpehh_missprob, fst_hitprob, fst_missprob, deldaf_hitprob, deldaf_missprob;
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double compLike, minDaf, minGenLen;
	double prior; // = 1/nSNP for region for now
	double gendist;
	double joint_score_prob_sel, joint_score_prob_neut; //useful for re-arranging CMS calculation to use prior correctly
	double cms_numerator, cms_denominator;
	int isnp, iComp, itoken, thisPos, likesFreqIndex, nComparisons, maxPos, minPos;
	int nsnps_regional;
	int istart, iend;
	int proceed; //Boolean used to log whether each SNP passes filter 0T 1F
	int takeIhs, takeDelihh, takeNsl, takeXpehh, takeFst, takeDeldaf; //Bools as above
	//int ibin;  //for debug

	if (argc <= 3) {
		fprintf(stderr, "Usage: ./combine_scores_local <savefilename> <cms_run_paramfile> <input_pair_file1> ...\n");
		exit(0);
	}
	nComparisons = argc - 3;
	
	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "\nPreparing to load component scores...\n");
	get_popComp_anyData(&score_data, nComparisons, argc, argv); 
	fprintf(stderr, "\tloaded data object with %d snps and %d population comparisons.\n", score_data.nsnps, score_data.ncomp);
	
	//for (isnp = 0; isnp < 5; isnp++ ){fprintf(stderr, "%d\t%f\t%f\n", score_data.physpos[1][isnp], score_data.ihs_normed[1][isnp], score_data.nsl_normed[1][isnp]);} // DEBUG

	////////////////////////////////////////
	// LOAD SCORE LIKELIHOODS (DEM MODEL) //
	// AND OTHER RUN PARAMETERS ////////////
	////////////////////////////////////////
	fprintf(stderr, "Preparing to load score likelihoods and composite parameters...\n");
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
	minGenLen = .5; //take no regions < .5 cM
	takeIhs = takeDelihh = takeNsl = takeXpehh = takeFst = takeDeldaf = 0; //all T by default
	//if additional lines is included, parse 
	if (fgets(paramline, line_size, inf) != NULL){
		for (running = paramline, itoken=0; (token = strsep(&running, " \t")) != NULL; itoken++){
			if (itoken == 0) {minPos = atoi(token);}
			else if (itoken == 1){maxPos = atoi(token);}
			else if (itoken == 2){minDaf = atof(token);}
			//else if (itoken == 3){writeLikes = atof(token);}	
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
	fprintf(stderr, "\tloaded parameters: minPos %d maxPos %d minDaf %f minGenLen %f\n", minPos, maxPos, minDaf, minGenLen);		
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
	/*
	if (minGenLen > 0){ //enforce minimum genetic length of region, if specified
		//get index of current bounds
		istart = 0;
		iend = 0;
		for (isnp = 0; isnp < score_data.nsnps; isnp++){
			thisPos = score_data.physpos[0][isnp];
			if (thisPos <= minPos){istart=isnp;}
			if (thisPos >= maxPos){iend=isnp; break;}
		}

		gendist = score_data.genpos[0][iend] - score_data.genpos[0][istart];
		//fprintf(stderr,"starting gendist: %f\n", gendist);
		while(gendist < minGenLen){
			istart-=1;  //move out one SNP on each end
			iend+=1;
			gendist = score_data.genpos[0][iend] - score_data.genpos[0][istart];
		}	
		
		minPos = score_data.physpos[0][istart];
		maxPos = score_data.physpos[0][iend];
		//fprintf(stderr," %d\t%d\n", minPos, maxPos);
		//fprintf(stderr, "adjusted region bounds to enforce minimum genetic length: %f\n", minGenLen);
	} // end adjust region bounds
	*/
	nsnps_regional = 0;
	for (isnp = 0; isnp < score_data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < score_data.ncomp; iComp++){
			if (score_data.physpos[iComp][isnp] != 0){break;}
		} //advance to the first comparison for which we have any data
		if (iComp >= score_data.ncomp){iComp = 0;} //catch SNPs at position 0
		thisihs = -1 * score_data.ihs_normed[iComp][isnp];
		thisihh = score_data.delihh_normed[iComp][isnp];
		thisnsl = -1 * score_data.nsl_normed[iComp][isnp];
		thisxpehh = compareXp(&score_data, isnp); //determine others by comparison. XP: take maximum.
		thisfst = compareFst_PBS(&score_data, isnp);
		thisdelDaf = -1 * comparedelDaf_outgroup_ave(&score_data, isnp);	 
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = 1 - score_data.daf_selpop[iComp][isnp];
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
	fprintf(stderr, "Preparing to write to: %s\n", outfilename);
	fprintf(stderr, "debug: score_data.nsnps = %d\n", score_data.nsnps);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(outf, "physPos\tgenPos\tpopDAF\tnormed_iHS\tnormed_deliHH\tnormed_nsl\tnormed_xp-ehh\tfst\tdelDAF\tcompLike_CMS\n");
	for (isnp = 0; isnp < score_data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < score_data.ncomp; iComp++){
			if (score_data.physpos[iComp][isnp] != 0){break;}
		} //advance to the first comparison for which we have any data
		if (iComp >= score_data.ncomp){iComp = 0;} //catch SNPs at position 0
		thisihs = -1 * score_data.ihs_normed[iComp][isnp];
		thisihh = score_data.delihh_normed[iComp][isnp];
		thisnsl = -1 * score_data.nsl_normed[iComp][isnp];
		thisxpehh = compareXp(&score_data, isnp);
		thisfst = compareFst_PBS(&score_data, isnp);	
		thisdelDaf = -1 * comparedelDaf_outgroup_ave(&score_data, isnp);	
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = 1 - score_data.daf_selpop[iComp][isnp];
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

			delihh_hitprob = getHitProb(&delihh_likes_data, likesFreqIndex, thisihh); //build in a check for pseudocounts here?
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

			///////////////////////////
			/// GET LOCAL CMS SCORE ///
			///////////////////////////
			joint_score_prob_sel = 1;
			joint_score_prob_neut = 1;

			if(takeIhs == 0 && isnan(thisihs) == 0){	joint_score_prob_sel *= ihs_hitprob; 
								joint_score_prob_neut *= ihs_missprob;	}	
			if(takeDelihh == 0 && isnan(thisihh) == 0){	joint_score_prob_sel *= delihh_hitprob; 
								joint_score_prob_neut *= delihh_missprob;	}
			if(takeNsl == 0 && isnan(thisnsl) == 0){	joint_score_prob_sel *= nsl_hitprob; 
								joint_score_prob_neut *= nsl_missprob;	}	
			if(takeFst == 0 && isnan(thisfst) == 0){	joint_score_prob_sel *= fst_hitprob; 
								joint_score_prob_neut *= fst_missprob;	}	
			if(takeDeldaf == 0 && isnan(thisdelDaf) == 0){	joint_score_prob_sel *= deldaf_hitprob; 
								joint_score_prob_neut *= deldaf_missprob;	}	
			if(takeXpehh == 0 && isnan(thisxpehh) == 0){	joint_score_prob_sel *= xpehh_hitprob; 
								joint_score_prob_neut *= xpehh_missprob;	}				
			cms_numerator = prior * joint_score_prob_sel;
			cms_denominator = cms_numerator + (1. - prior) * joint_score_prob_neut;
			compLike = cms_numerator / cms_denominator;			
			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf, compLike);

			//DEBUG 
			/*fprintf(stderr, "%d\t%d\n", cms_numerator, cms_denominator);
			fprintf(stderr, "ihs %f\t hit %e\tmiss %e\n", thisihs, ihs_hitprob, ihs_missprob); 
			fprintf(stderr, "delihh %f\t hit %e\tmiss %e\n", thisihh, delihh_hitprob, delihh_missprob);
			fprintf(stderr, "nsl %f\t hit %e\tmiss %e\n", thisnsl, nsl_hitprob, nsl_missprob);			
			fprintf(stderr, "fst %f\t hit %e\tmiss %e\n", thisfst, fst_hitprob, fst_missprob); 
			fprintf(stderr, "deldaf %f\t hit %e\tmiss %e\n", thisdelDaf, deldaf_hitprob, deldaf_missprob); 
			fprintf(stderr, "xp %f\t hit %e\tmiss %e\n", thisxpehh, xpehh_hitprob, xpehh_missprob); 
			fprintf(stderr, "cl: %e\n", compLike);
			fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf);
			*/
		}//end if-a-go
	} // end isnp
	fclose(outf);
	fprintf(stderr, "\nWrote to %s\n", outfilename);
	free_popComp_data_multiple(&score_data);
	free_likes_data_multiple(&ihs_likes_data);
	free_likes_data_multiple(&nsl_likes_data);
	free_likes_data_multiple(&delihh_likes_data);
	free_likes_data_multiple(&xpehh_likes_data);
	free_likes_data_multiple(&fst_likes_data);				
	free_likes_data_multiple(&deldaf_likes_data);		
	return 0;
} // end main