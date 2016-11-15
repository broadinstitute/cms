// methods for running CMS with a putative selPop and 2+ outgroups. 
// last updated 11.15.16 	vitti@broadinstitute.org

#define EXTRAARGS 32
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
//#include <zlib.h>
#include "cms_data.h"
#include "pop_comparison.h"

float getMinBf(likes_data* data_hit, likes_data* data_miss){
	int ibin;
	float thisBf;
	float minBf = 1;

	for (ibin = 0; ibin < data_hit->nbins; ibin++){
		if(data_hit->probs[ibin] != 1e-10 && data_miss->probs[ibin] != 1e-10){
			thisBf = data_hit->probs[ibin] / data_miss->probs[ibin];
			if (thisBf < minBf){minBf = thisBf;}
		}
	}//end ibin
	
	return minBf;
}//end function

float getProb(likes_data* data, double value){
	int ibin;
	for (ibin = 0; ibin < data->nbins; ibin++){
		if (value >= data->start_bin[ibin] && value <= data->end_bin[ibin]){return data->probs[ibin];}
	}
	if (value < data->start_bin[0]){return data->probs[0];}
	if (value > data->end_bin[data->nbins - 1]){return data->probs[data->nbins - 1];}

	return 0;
} //end function
float compareXp(popComp_data_multiple* data, int isnp){//currently: takes max val
	double xp;
	int iComp;
	xp = -1e10;
	for (iComp = 0; iComp < data->ncomp; iComp++){
		if (data->xp_normed[iComp][isnp] > xp){xp = data->xp_normed[iComp][isnp];}
	}
	return xp;
} //end function
float compareFst(popComp_data_multiple* data, int isnp){ //currently: takes average
//takes LSBL? or PBS? I would need to include outgroup-pairs-Fst. 
//previously: takes max val
	double fst;
	double ave;
	int iComp;
	int theseComp=0;
	fst =	0.;
	for (iComp = 0; iComp < data->ncomp; iComp++){
		if (data->fst[iComp][isnp] != 0){
		fst += data->fst[iComp][isnp];
		theseComp ++;}
	}
	ave = fst / (double)theseComp;
	return ave;
} //end function
float comparedelDaf(popComp_data_multiple* data, int isnp){//currently: takes average
//previously: takes max val
	double deldaf;
	int iComp;
	double ave;
	int theseComp=0;
	deldaf = 0;//-100;
	for (iComp = 0; iComp < data->ncomp; iComp++){
		if (data->delDAF[iComp][isnp] !=0){
		deldaf += data->delDAF[iComp][isnp];
		theseComp++;}
		//if (data->delDAF[iComp][isnp] > deldaf){deldaf = data->delDAF[iComp][isnp];}

	}
	ave = deldaf / (double)theseComp;
	return ave;
} //end function

void get_popComp_data_multiple(popComp_data_multiple* data, int argc, char *argv[]){
	const int line_size = 15000000; 
	popComp_data data_sing;
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilename[512];
	int	isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique;
	int *allSnps, *allUniqueSnps;
	int numLikesFiles;

	//////////////////
	/// INITIALIZE ///
	//////////////////
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 
	data->nsnps = 0;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

	//////////////////////////
	/// COLLATE LOCI: LOAD ///
	/////////////////////////

	nComparisons = argc - EXTRAARGS;  
	numLikesFiles = EXTRAARGS;
	//fprintf(stderr,"\n\n%d\n\n", nComparisons);
	totNsnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+(numLikesFiles)]);
		inf = fopen(infilename, "r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header
		while (fgets(newLine, line_size, inf) != NULL){totNsnp++;}
		fclose(inf);
	} // end iComp
	//fprintf(stderr, "Found %d loci \n", totNsnp);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+(numLikesFiles)]);
		inf = fopen(infilename,"r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header

		while (fgets(newLine, line_size, inf) != NULL){
			for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 1) {
				 thisPhysPos = atoi(token);
				 allSnps[isnp] =thisPhysPos;
				 break;
				 } // end if
			} // end for 
			isnp +=1;	 
		} // end while
		fclose(inf);
	}//end icomp

	//////////////////////////
	/// COLLATE LOCI: SORT ///
	/////////////////////////

	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	fprintf(stderr, "Sorted SNPs from pos %d to %d\n", allSnps[0], allSnps[totNsnp-1]);
	nunique = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{nunique++;}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison...\n", nunique);

	allUniqueSnps	= malloc(nunique * sizeof(int));
	jsnp = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
	} // end for isnp

	///////////////////////
	/// ALLOCATE MEMORY ///
	//////////////////////

	fprintf(stderr, "Allocating memory...\n");
	data->nsnps = nunique;
	data->ncomp =nComparisons;
	data->physpos = malloc(nComparisons * sizeof(int*));
	data->genpos = malloc(nComparisons * sizeof(double*));
	data->daf_selpop = malloc(nComparisons * sizeof(double*));
	data->delDAF = malloc(nComparisons * sizeof(double*));
	data->fst = malloc(nComparisons * sizeof(double*));
	data->xp_normed = malloc(nComparisons * sizeof(double*));
	data->ihs_normed = malloc(nComparisons * sizeof(double*));
	data->delihh_normed = malloc(nComparisons * sizeof(double*));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);

	for (iComp = 0; iComp < nComparisons; iComp++){
		data->physpos[iComp] = calloc(nunique, sizeof(int));
		data->genpos[iComp] = calloc(nunique, sizeof(double));
		data->daf_selpop[iComp] = calloc(nunique, sizeof(double));
		data->delDAF[iComp] = calloc(nunique, sizeof(double));
		data->fst[iComp] = calloc(nunique, sizeof(double));
		data->xp_normed[iComp] = calloc(nunique, sizeof(double));		
		data->ihs_normed[iComp] = calloc(nunique, sizeof(double));
		data->delihh_normed[iComp] = calloc(nunique, sizeof(double));		
		assert(data->physpos[iComp] != NULL);
		assert(data->genpos[iComp] != NULL);
		assert(data->daf_selpop[iComp] != NULL);
		assert(data->delDAF[iComp] != NULL);
		assert(data->fst[iComp] != NULL);
		assert(data->xp_normed[iComp] != NULL);
		assert(data->ihs_normed[iComp] != NULL);
		assert(data->delihh_normed[iComp] != NULL);
	} // end for icomp

	/////////////////////////////////////////////
	// LOAD ALL COMPARISONS TO ONE DATA OBJECT //
	/////////////////////////////////////////////

	fprintf(stderr, "Loading all component scores...\n");
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+numLikesFiles]);
		//fprintf(stderr, infilename);
		get_popComp_data(&data_sing, infilename);
		jsnp = 0; //isnp iterates (0, nunique) over allUnique Snps; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
			//fprintf(stderr, "%d\t%d\t%d\t%d\n", isnp, jsnp, allUniqueSnps[isnp], data_sing.physpos[jsnp]);
			if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
				data->physpos[iComp][isnp] = data_sing.physpos[jsnp];	
				data->genpos[iComp][isnp] = data_sing.genpos[jsnp];	 
				data->daf_selpop[iComp][isnp] = data_sing.daf_selpop[jsnp];	 
				data->delDAF[iComp][isnp] = data_sing.delDAF[jsnp];	
				data->fst[iComp][isnp] = data_sing.fst[jsnp];	 
				data->xp_normed[iComp][isnp] = data_sing.xp_normed[jsnp];							 
				data->ihs_normed[iComp][isnp] = data_sing.ihs_normed[jsnp];	 
				data->delihh_normed[iComp][isnp] = data_sing.delihh_normed[jsnp];			 
				jsnp++; //assert(jsnp<=data_sing.nsnps);
				if (jsnp >= data_sing.nsnps){break;}
			}
			else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; if (jsnp >= data_sing.nsnps){break;}}//assert(jsnp<=data_sing.nsnps);}
			//else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
		}// end for isnp loop

		free_popComp_data(&data_sing); 
	} // end for icomp
} //end method
void get_popComp_data_multiple_region(popComp_data_multiple* data, int argc, char *argv[]){
	const int line_size = 15000000; 
	popComp_data data_sing;
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilename[512];
	int	isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique;
	int *allSnps, *allUniqueSnps;
	int numLikesFiles;
	int startPos, endPos;

	//////////////////
	/// INITIALIZE ///
	//////////////////
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 
	data->nsnps = 0;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

	startPos = atoi(argv[1]);
	endPos = atoi(argv[2]);

	//////////////////////////
	/// COLLATE LOCI: LOAD ///
	/////////////////////////

	nComparisons = argc - (EXTRAARGS+2); 
	numLikesFiles = (EXTRAARGS+2);
	totNsnp = 0;
	//fprintf(stderr, "numcomparisons: %d\n", nComparisons);
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+(numLikesFiles)]);
		//fprintf(stderr, infilename);
		inf = fopen(infilename, "r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header
		while (fgets(newLine, line_size, inf) != NULL){totNsnp++;}
		fclose(inf);
	} // end iComp
	//fprintf(stderr, "Found %d loci \n", totNsnp);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+(numLikesFiles)]);
		inf = fopen(infilename,"r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header

		while (fgets(newLine, line_size, inf) != NULL){
			for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 1) {
				 thisPhysPos = atoi(token);
				 allSnps[isnp] =thisPhysPos;
				 break;
				 } // end if
			} // end for 
			isnp +=1;	 
		} // end while
		fclose(inf);
	}//end icomp

	//////////////////////////
	/// COLLATE LOCI: SORT ///
	/////////////////////////

	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	fprintf(stderr, "Sorted SNPs from pos %d to %d\n", allSnps[0], allSnps[totNsnp-1]);
	nunique = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] >= startPos && allSnps[isnp] <= endPos){
			if (allSnps[isnp] == allSnps[isnp+1]){continue;}
			else{nunique++;}
		}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison...\n", nunique);

	allUniqueSnps	= malloc(nunique * sizeof(int));
	jsnp = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] >= startPos && allSnps[isnp] <= endPos){
			if (allSnps[isnp] == allSnps[isnp+1]){continue;}
			else{allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
		}
	} // end for isnp

	///////////////////////
	/// ALLOCATE MEMORY ///
	//////////////////////

	fprintf(stderr, "Allocating memory...\n");
	data->nsnps = nunique;
	data->ncomp =nComparisons;
	data->physpos = malloc(nComparisons * sizeof(int*));
	data->genpos = malloc(nComparisons * sizeof(double*));
	data->daf_selpop = malloc(nComparisons * sizeof(double*));
	data->delDAF = malloc(nComparisons * sizeof(double*));
	data->fst = malloc(nComparisons * sizeof(double*));
	data->xp_normed = malloc(nComparisons * sizeof(double*));
	data->ihs_normed = malloc(nComparisons * sizeof(double*));
	data->delihh_normed = malloc(nComparisons * sizeof(double*));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);

	for (iComp = 0; iComp < nComparisons; iComp++){
		data->physpos[iComp] = calloc(nunique, sizeof(int));
		data->genpos[iComp] = calloc(nunique, sizeof(double));
		data->daf_selpop[iComp] = calloc(nunique, sizeof(double));
		data->delDAF[iComp] = calloc(nunique, sizeof(double));
		data->fst[iComp] = calloc(nunique, sizeof(double));
		data->xp_normed[iComp] = calloc(nunique, sizeof(double));		
		data->ihs_normed[iComp] = calloc(nunique, sizeof(double));
		data->delihh_normed[iComp] = calloc(nunique, sizeof(double));		
		assert(data->physpos[iComp] != NULL);
		assert(data->genpos[iComp] != NULL);
		assert(data->daf_selpop[iComp] != NULL);
		assert(data->delDAF[iComp] != NULL);
		assert(data->fst[iComp] != NULL);
		assert(data->xp_normed[iComp] != NULL);
		assert(data->ihs_normed[iComp] != NULL);
		assert(data->delihh_normed[iComp] != NULL);
	} // end for icomp

	/////////////////////////////////////////////
	// LOAD ALL COMPARISONS TO ONE DATA OBJECT //
	/////////////////////////////////////////////

	fprintf(stderr, "Loading all component scores...\n");
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp+numLikesFiles]);
		get_popComp_data(&data_sing, infilename);
		//fprintf(stderr, "\n\tloaded data for one population comparison; nSnps: %d\n", data_sing.nsnps);
		//jsnp = 0; //isnp iterates (0, nunique) over allUnique Snps; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
			for (jsnp = 0; jsnp < data_sing.nsnps; jsnp++){
				if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
					data->physpos[iComp][isnp] = data_sing.physpos[jsnp];	
					data->genpos[iComp][isnp] = data_sing.genpos[jsnp];	 
					data->daf_selpop[iComp][isnp] = data_sing.daf_selpop[jsnp];	 
					data->delDAF[iComp][isnp] = data_sing.delDAF[jsnp];	
					data->fst[iComp][isnp] = data_sing.fst[jsnp];	 
					data->xp_normed[iComp][isnp] = data_sing.xp_normed[jsnp];							 
					data->ihs_normed[iComp][isnp] = data_sing.ihs_normed[jsnp];	 
					data->delihh_normed[iComp][isnp] = data_sing.delihh_normed[jsnp];			 
					break;
					//jsnp++; //assert(jsnp<=data_sing.nsnps);
					//if (jsnp >= data_sing.nsnps){break;}
				}
				//else {pass;}
				//else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){jsnp++; if (jsnp >= data_sing.nsnps){break;}}//assert(jsnp<=data_sing.nsnps);}
				//else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
				} //end for jsnp loop

		}// end for isnp loop

		free_popComp_data(&data_sing); 
	} // end for icomp
} //end method
void free_popComp_data_multiple(popComp_data_multiple* data){
	int iComp;
	if (data == NULL) {return;}
	for (iComp = 0; iComp < data->ncomp; iComp++){
		free(data->physpos[iComp]);
		free(data->genpos[iComp]);
		free(data->daf_selpop[iComp]);
		free(data->delDAF[iComp]);
		free(data->fst[iComp]);
		free(data->xp_normed[iComp]);
		free(data->ihs_normed[iComp]);
		free(data->delihh_normed[iComp]);
	}
	free(data->physpos);
	free(data->genpos);
	free(data->daf_selpop);
	free(data->delDAF);
	free(data->fst);
	free(data->xp_normed);
	free(data->ihs_normed);
	free(data->delihh_normed);
	data->nsnps = 0;
	data->ncomp = 0;
} //end method
