// 	this implementation calculates a Bayesian posterior for /ALL/ SNPs (gw) using 
//	a uniform prior (experimental)	10.31.2017

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
	likes_data_multiple H12_likes_data, iSAFE_likes_data;
	FILE *inf=NULL, *outf=NULL, *outf2=NULL;	
	char *token, *running;
	char cms_param_filename[528], paramline[528], outfilename[256], outfilename_likes[256];
	char ihs_master_likesfilename[256], nsl_master_likesfilename[256], delihh_master_likesfilename[256];
    char xpehh_master_likesfilename[256], fst_master_likesfilename[256], deldaf_master_likesfilename[256];
    char H12_master_likesfilename[256], iSAFE_master_likesfilename[256];
	float delihh_hitprob, delihh_missprob, delihh_bf, delihh_minbf, delihh_maxbf; 
	float nsl_hitprob, nsl_missprob, nsl_bf, nsl_minbf, nsl_maxbf; //bayes factor
	float ihs_hitprob, ihs_missprob, ihs_bf, ihs_minbf, ihs_maxbf;
	float xpehh_hitprob, xpehh_missprob, xpehh_bf, xpehh_minbf, xpehh_maxbf;
	float fst_hitprob, fst_missprob, fst_bf, fst_minbf, fst_maxbf;
	float deldaf_hitprob, deldaf_missprob, deldaf_bf, deldaf_minbf, deldaf_maxbf;
	float H12_hitprob, H12_missprob, H12_bf, H12_minbf, H12_maxbf;
	float iSAFE_hitprob, iSAFE_missprob, iSAFE_bf, iSAFE_minbf, iSAFE_maxbf;
	int isnp, iComp, itoken, thisPos, likesFreqIndex, nComparisons, maxPos, minPos;
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double thisH12, thisH2H1, thisiSAFE;
	double compLike, minDaf;
	double joint_score_prob_sel, joint_score_prob_neut; //useful for re-arranging CMS calculation to use prior correctly
	double cms_numerator, cms_denominator;

	int proceed; //Boolean used to log whether each SNP passes filter 0T 1F
	int takeIhs, takeDelihh, takeNsl, takeXpehh, takeFst, takeDeldaf, writeLikes; //Bools as above
	int takeH12, takeiSAFE;
	//int ibin;  //for debug
	double prior; // uniform for naive CMS_gw_posterior
	
	if (argc <= 3) {
		fprintf(stderr, "Usage: ./combine_gw_posterior <savefilename> <cms_run_paramfile> <input_pair_file1> ...\n");
		exit(0);
	}
	nComparisons = argc - 3;
	
	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "\nPreparing to load component scores...\n");
	get_popComp_anyData(&score_data, nComparisons, argc, argv); 
	fprintf(stderr, "\tloaded data object with %d snps and %d population comparisons.\n", score_data.nsnps, score_data.ncomp);
	//for (isnp = 0; isnp < score_data.nsnps; isnp++ ){fprintf(stderr, "%f\t", score_data.ihs_normed[1][isnp]);} // DEBUG

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
	fgets(H12_master_likesfilename, line_size, inf);
	strtok(H12_master_likesfilename, "\n");
	fgets(iSAFE_master_likesfilename, line_size, inf);
	strtok(iSAFE_master_likesfilename, "\n");
	fgets(xpehh_master_likesfilename, line_size, inf);
	strtok(xpehh_master_likesfilename, "\n");
	fgets(fst_master_likesfilename, line_size, inf);
	strtok(fst_master_likesfilename, "\n");
	fgets(deldaf_master_likesfilename, line_size, inf);
	strtok(deldaf_master_likesfilename, "\n");

	//STANDARDIZE ORDERING
	//LOSE LIKESFREQS

	//set defaults
	minPos = -1;		
	maxPos = 2147483647;
	minDaf = 0;
	takeIhs = takeDelihh = takeNsl = takeXpehh = takeFst = takeDeldaf = 0; //all T by default
	takeH12 = takeiSAFE = 0;
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
			else if (itoken == 1){takeNsl = atoi(token);}
			else if (itoken == 2){takeDelihh = atoi(token);}			
			else if (itoken == 3){takeH12 = atoi(token);}
			else if (itoken == 4){takeiSAFE = atoi(token);}			
			else if (itoken == 5){takeXpehh = atoi(token);}	
			else if (itoken == 6){takeFst = atoi(token);}			
			else if (itoken == 7){takeDeldaf = atoi(token);}	


		} // end for running
	}  //end if fgets paramline
	fclose(inf);
	fprintf(stderr, "\tloaded parameters: minPos %d maxPos %d minDaf %f\n", minPos, maxPos, minDaf);		
	get_likes_data_multiple(&ihs_likes_data, ihs_master_likesfilename); 
	get_likes_data_multiple(&nsl_likes_data, nsl_master_likesfilename); 
	get_likes_data_multiple(&delihh_likes_data, delihh_master_likesfilename); 
	get_likes_data_multiple(&H12_likes_data, H12_master_likesfilename); 
	get_likes_data_multiple(&iSAFE_likes_data, iSAFE_master_likesfilename); 
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
	fprintf(outf, "physPos\tgenPos\tpopDAF\tnormed_iHS\tnormed_nsl\tnormed_delihh\tH12\tiSAFE\tnormed_xp-ehh\tfst\tdelDAF\tcompLike_CMS\n");
	/*
	if (writeLikes == 0){
		strcpy(outfilename_likes, argv[1]);
		strcat(outfilename_likes, ".likes");
		outf2 = fopen(outfilename_likes, "w");
		assert(outf2 != NULL);
		fprintf(outf2, "physPos\tgenPos\tpopDAF\tlike_iHS\tlike_nsl\tlike_deliHH\tlike_H12\tlike_iSAFE\tlike_xp-ehh\tlikefst\tlikedelDAF\tcompLike_CMS\n");	
	} //end if write likes
	*/
	for (isnp = 0; isnp < score_data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < score_data.ncomp; iComp++){
			if (score_data.physpos[iComp][isnp] != 0){break;}
		} //advance to the first comparison for which we have any data
		if (iComp >= score_data.ncomp){iComp = 0;} //catch SNPs at position 0

		thisihs = score_data.ihs_normed[iComp][isnp];
		thisnsl = score_data.nsl_normed[iComp][isnp];
		thisihh = score_data.delihh_normed[iComp][isnp];
		thisH12 = score_data.H12[iComp][isnp];
		thisH2H1 = score_data.H2H1[iComp][isnp];
		thisiSAFE = score_data.iSAFE[iComp][isnp];
		thisxpehh = compareXp(&score_data, isnp);
		thisfst = compareFst_PBS(&score_data, isnp);
		thisdelDaf = comparedelDaf_outgroup_ave(&score_data, isnp);
		
		proceed = 0;
		//check position
		thisPos = score_data.physpos[iComp][isnp];
		if (thisPos < minPos){proceed=1;}
		if (thisPos > maxPos){proceed=1;}
		//check daf
		thisdaf = score_data.daf_selpop[iComp][isnp];
		if (thisdaf < minDaf){proceed=1;} 
		compLike = 0;
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
			H12_hitprob = getHitProb(&H12_likes_data, likesFreqIndex, thisH12);
			iSAFE_hitprob = getHitProb(&iSAFE_likes_data, likesFreqIndex, thisiSAFE);

			delihh_missprob = getMissProb(&delihh_likes_data, thisihh);
			nsl_missprob = getMissProb(&nsl_likes_data, thisnsl);			
			ihs_missprob = getMissProb(&ihs_likes_data, thisihs);
			fst_missprob = getMissProb(&fst_likes_data, thisfst);
			deldaf_missprob = getMissProb(&deldaf_likes_data, thisdelDaf);
			xpehh_missprob = getMissProb(&xpehh_likes_data, thisxpehh);
			H12_missprob = getMissProb(&H12_likes_data, thisH12);
			iSAFE_missprob = getMissProb(&iSAFE_likes_data, thisiSAFE);

			/*
			delihh_minbf = getMinBf(&delihh_likes_data, likesFreqIndex);
			nsl_minbf = getMinBf(&nsl_likes_data, likesFreqIndex);
			ihs_minbf = getMinBf(&ihs_likes_data, likesFreqIndex);	
			fst_minbf = getMinBf(&fst_likes_data, likesFreqIndex);
			deldaf_minbf = getMinBf(&deldaf_likes_data, likesFreqIndex);
			xpehh_minbf = getMinBf(&xpehh_likes_data, likesFreqIndex);		
			H12_minbf = getMinBf(&H12_likes_data, likesFreqIndex);
			iSAFE_minbf = getMinBf(&iSAFE_likes_data, likesFreqIndex);			
			
			delihh_maxbf = getMaxBf(&delihh_likes_data, likesFreqIndex);
			nsl_maxbf = getMaxBf(&nsl_likes_data, likesFreqIndex);
			ihs_maxbf = getMaxBf(&ihs_likes_data, likesFreqIndex);	
			fst_maxbf = getMaxBf(&fst_likes_data, likesFreqIndex);
			deldaf_maxbf = getMaxBf(&deldaf_likes_data, likesFreqIndex);
			xpehh_maxbf = getMaxBf(&xpehh_likes_data, likesFreqIndex);			
			H12_maxbf = getMaxBf(&H12_likes_data, likesFreqIndex);
			iSAFE_maxbf = getMaxBf(&iSAFE_likes_data, likesFreqIndex);	

			delihh_bf = 1;
			if (delihh_missprob > 2e-10 && delihh_hitprob > 2e-10){delihh_bf = delihh_hitprob / delihh_missprob;}
			if (delihh_missprob < 2e-10 && delihh_hitprob > 2e-10){delihh_bf = delihh_maxbf;}
			if (delihh_hitprob < 2e-10 && delihh_missprob > 2e-10){delihh_bf = delihh_minbf;}

			nsl_bf = 1;
			if (nsl_missprob > 2e-10 && nsl_hitprob > 2e-10){nsl_bf = nsl_hitprob / nsl_missprob;}
			if (nsl_missprob < 2e-10 && nsl_hitprob > 2e-10){nsl_bf = nsl_maxbf;}
			if (nsl_hitprob < 2e-10 && nsl_missprob > 2e-10){nsl_bf = nsl_minbf;}

			ihs_bf = 1;
			if (ihs_missprob > 2e-10 && ihs_hitprob > 2e-10){ihs_bf = ihs_hitprob / ihs_missprob;}
			if (ihs_missprob < 2e-10 && ihs_hitprob > 2e-10){ihs_bf = ihs_maxbf;}
			if (ihs_hitprob < 2e-10 && ihs_missprob > 2e-10){ihs_bf = ihs_minbf;}

			fst_bf = 1;
			if (fst_missprob > 2e-10 && fst_hitprob > 2e-10){fst_bf = fst_hitprob / fst_missprob;}
			if (fst_missprob < 2e-10 && fst_hitprob > 2e-10){fst_bf = fst_maxbf;}
			if (fst_hitprob < 2e-10 && fst_missprob > 2e-10){fst_bf = fst_minbf;}

			deldaf_bf = 1;
			if (deldaf_missprob > 2e-10 && deldaf_hitprob > 2e-10){deldaf_bf = deldaf_hitprob / deldaf_missprob;}
			if (deldaf_missprob < 2e-10 && deldaf_hitprob > 2e-10){deldaf_bf = deldaf_maxbf;}
			if (deldaf_hitprob < 2e-10 && deldaf_missprob > 2e-10){deldaf_bf = deldaf_minbf;}
			
			xpehh_bf = 1;
			if (xpehh_missprob > 2e-10 && xpehh_hitprob > 2e-10){xpehh_bf = xpehh_hitprob / xpehh_missprob;}			
			if (xpehh_missprob < 2e-10 && xpehh_hitprob > 2e-10){xpehh_bf = xpehh_maxbf;}
			if (xpehh_hitprob < 2e-10 && xpehh_missprob > 2e-10){xpehh_bf = xpehh_minbf;}

			H12_bf = 1;
			if (H12_missprob > 2e-10 && H12_hitprob > 2e-10){H12_bf = H12_hitprob / H12_missprob;}
			if (H12_missprob < 2e-10 && H12_hitprob > 2e-10){H12_bf = H12_maxbf;}
			if (H12_hitprob < 2e-10 && H12_missprob > 2e-10){H12_bf = H12_minbf;}
			
			iSAFE_bf = 1;
			if (iSAFE_missprob > 2e-10 && iSAFE_hitprob > 2e-10){iSAFE_bf = iSAFE_hitprob / iSAFE_missprob;}			
			if (iSAFE_missprob < 2e-10 && iSAFE_hitprob > 2e-10){iSAFE_bf = iSAFE_maxbf;}
			if (iSAFE_hitprob < 2e-10 && iSAFE_missprob > 2e-10){iSAFE_bf = iSAFE_minbf;}
			*/


					
			//////////////////////////////////////
			/// GET CMS SCORE - GW / POSTERIOR ///
			//////////////////////////////////////
			prior = .5;
			joint_score_prob_sel = 1;
			joint_score_prob_neut = 1;

			if(takeIhs == 0 && isnan(thisihs) == 0){	joint_score_prob_sel *= ihs_hitprob; 
								joint_score_prob_neut *= ihs_missprob;	}	
			if(takeDelihh == 0 && isnan(thisihh) == 0){	joint_score_prob_sel *= delihh_hitprob; 
								joint_score_prob_neut *= delihh_missprob;	}
			if(takeNsl == 0 && isnan(thisnsl) == 0){	joint_score_prob_sel *= nsl_hitprob; 
								joint_score_prob_neut *= nsl_missprob;	}	
			if(takeiSAFE == 0 && isnan(thisiSAFE) == 0){	joint_score_prob_sel *= iSAFE_hitprob; 
								joint_score_prob_neut *= iSAFE_missprob;	}	
			if(takeH12 == 0 && isnan(thisH12) == 0){	joint_score_prob_sel *= H12_hitprob; 
								joint_score_prob_neut *= H12_missprob;	}	
			if(takeFst == 0 && isnan(thisfst) == 0){	joint_score_prob_sel *= fst_hitprob; 
								joint_score_prob_neut *= fst_missprob;	}	
			if(takeDeldaf == 0 && isnan(thisdelDaf) == 0){	joint_score_prob_sel *= deldaf_hitprob; 
								joint_score_prob_neut *= deldaf_missprob;	}	
			if(takeXpehh == 0 && isnan(thisxpehh) == 0){	joint_score_prob_sel *= xpehh_hitprob; 
								joint_score_prob_neut *= xpehh_missprob;	}		

			cms_numerator = prior * joint_score_prob_sel;
			cms_denominator = cms_numerator + (1. - prior) * joint_score_prob_neut;

			compLike = cms_numerator / cms_denominator;			
			
			//fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf, compLike);
			

			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, thisihs, thisnsl,  thisihh, thisH12, thisiSAFE, thisxpehh, thisfst, thisdelDaf, compLike);
			
			/*

			/////////////////////
			/// GET CMS SCORE ///
			/////////////////////		
			if(takeIhs == 0){compLike *= ihs_bf;}			//;fprintf(stderr, "ihs\t");}
			if(takeDelihh == 0){compLike *= delihh_bf;}	//;fprintf(stderr, "delihh\t");}
			if(takeNsl == 0){compLike *= nsl_bf;}			//;fprintf(stderr, "nsl\t");}
			if(takeFst == 0){compLike *= fst_bf;}			//;fprintf(stderr, "fst\t");}
			if(takeDeldaf == 0){compLike *= deldaf_bf;}	//;fprintf(stderr, "deldaf\t");}
			if(takeXpehh == 0){compLike *= xpehh_bf;}		//;fprintf(stderr, "xp\n");}
			if(takeH12 == 0){compLike *= H12_bf;}			
			if(takeiSAFE == 0){compLike *= iSAFE_bf;}		


			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, thisihs, thisnsl,  thisihh, thisH12, thisiSAFE, thisxpehh, thisfst, thisdelDaf, compLike);
			if (writeLikes == 0){fprintf(outf2, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", score_data.physpos[iComp][isnp], score_data.genpos[iComp][isnp], thisdaf, ihs_bf, nsl_bf, delihh_bf, H12_bf, iSAFE_bf, xpehh_bf, fst_bf, deldaf_bf, compLike);} //end if write likes
			*/

			//DEBUG 
			/*
			fprintf(stderr, "ihs %f\t hit %e\tmiss %e\tbf %e\n", thisihs, ihs_hitprob, ihs_missprob, ihs_bf); //debug
			fprintf(stderr, "delihh %f\t hit %e\tmiss %e\tbf %e\n", thisihh, delihh_hitprob, delihh_missprob, delihh_bf); //debug
			fprintf(stderr, "nsl %f\t hit %e\tmiss %e\tbf %e\n", thisnsl, nsl_hitprob, nsl_missprob, nsl_bf); //debug
			fprintf(stderr, "H12 %f\t hit %e\tmiss %e\tbf %e\n", thisH12, H12_hitprob, H12_missprob, H12_bf); //debug
			fprintf(stderr, "iSAFE %f\t hit %e\tmiss %e\tbf %e\n", thisiSAFE, iSAFE_hitprob, iSAFE_missprob, iSAFE_bf); //debug
			fprintf(stderr, "fst %f\t hit %e\tmiss %e\tbf %e\n", thisfst, fst_hitprob, fst_missprob, fst_bf); //debug
			fprintf(stderr, "xp %f\t hit %e\tmiss %e\tbf %e\n", thisxpehh, xpehh_hitprob, xpehh_missprob, xpehh_bf); //debug
			fprintf(stderr, "clr: %e\n", compLike);
			fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf);
			*/
		}//end if-a-go
	} // end isnp
	fclose(outf);
	fprintf(stderr, "\nWrote CMS scores to %s\n", outfilename);
	if (writeLikes == 0){fprintf(stderr, "Wrote score decomposition to %s\n", outfilename_likes);}
	free_popComp_data_multiple(&score_data);
	free_likes_data_multiple(&ihs_likes_data);
	free_likes_data_multiple(&nsl_likes_data);
	free_likes_data_multiple(&delihh_likes_data);
	free_likes_data_multiple(&xpehh_likes_data);
	free_likes_data_multiple(&fst_likes_data);				
	free_likes_data_multiple(&deldaf_likes_data);
	free_likes_data_multiple(&H12_likes_data);				
	free_likes_data_multiple(&iSAFE_likes_data);

	return 0;
} // end main