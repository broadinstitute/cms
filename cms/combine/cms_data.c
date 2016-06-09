// {{COMPONENT SCORES + SCORE DISTRIBUTIONS FROM SIMULATED DATA --> COMPOSITE SCORES}}
// functions for handling cms component(+composite) score datastructures
// last updated: 06.04.15 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>
#include "cms_data.h"

int intcmp(const void *v1, const void *v2);
int intcmp(const void *v1, const void *v2) {return (*(int *)v1 - *(int *)v2);}

/**********************/
/***COMPONENT SCORES***/
/**********************/

void get_fst_deldaf_data(fst_deldaf_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
		
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->fst = NULL;
	data->deldaf = NULL;

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: "); fprintf(stderr, filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->fst = malloc(data->nsnps * sizeof(double*)); assert(data->fst != NULL);
	data->deldaf = malloc(data->nsnps * sizeof(double*)); assert(data->deldaf != NULL);

	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	fgets(newLine, line_size, inf); // strip header
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 0) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 1) {
				data->fst[isnp] = atof(token);
			}
			else if (itoken == 2) {
				data->deldaf[isnp] = atof(token);
			}			
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method

void free_fst_deldaf_data(fst_deldaf_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->fst);
	free(data->deldaf);
	data->nsnps = 0;
} //end method

void get_delihh_data(delihh_data* data, char filename[]) {
	const int line_size = 15000000;	
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
		
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->ihs_unnormed = NULL;
	data->delihh_unnormed = NULL;
	data->delihh_normed = NULL;
	data->lastcol = NULL;

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: "); fprintf(stderr, filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihs_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_unnormed != NULL);
	data->delihh_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->delihh_unnormed != NULL);	
	data->delihh_normed = malloc(data->nsnps * sizeof(double*)); assert(data->delihh_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);
	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 1) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 2) {
				data->freq1[isnp] = atof(token);
			}
			else if (itoken == 3) {
				data->ihs_unnormed[isnp] = atof(token);
			}			
			else if (itoken == 4) {
				data->delihh_unnormed[isnp] = atof(token);
			}
			//5 is duplicate!!!!
			else if (itoken == 6) {
				data->delihh_normed[isnp] = atof(token);
			}		 
			else if (itoken == 7) {
				data->lastcol[isnp] = atoi(token);
			}						
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method

void free_delihh_data(delihh_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->freq1);
	free(data->ihs_unnormed);
	free(data->delihh_unnormed);
	free(data->delihh_normed);	
	free(data->lastcol);
	data->nsnps = 0;
} //end method

void get_xpehh_data(xpehh_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int	isnp, itoken;
		
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->genpos = NULL;
	data->freq1 = NULL;
	data->ihh1 = NULL; 
	data->freq2 = NULL;
	data->ihh2 = NULL;
	data->xpehh_unnormed = NULL; 
	data->xpehh_normed = NULL;
	data->lastcol = NULL;

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: "); fprintf(stderr, filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->genpos = malloc(data->nsnps * sizeof(double*)); assert(data->genpos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihh1 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh1 != NULL);
	data->freq2 = malloc(data->nsnps * sizeof(double*)); assert(data->freq2 != NULL);
	data->ihh2 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh2 != NULL);
	data->xpehh_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->xpehh_unnormed != NULL);
	data->xpehh_normed = malloc(data->nsnps * sizeof(double*)); assert(data->xpehh_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);

	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	fgets(newLine, line_size, inf); // strip header
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 1) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 2) {
				data->genpos[isnp] = atof(token);
			}
			else if (itoken == 3) {
				data->freq1[isnp] = atof(token);
			}		
			else if (itoken == 4) {
				data->ihh1[isnp] = atof(token);
			}
			else if (itoken == 5) {
				data->freq2[isnp] = atof(token);
			}		
			else if (itoken == 6) {
				data->ihh2[isnp] = atof(token);
			}
			else if (itoken == 7) {
				data->xpehh_unnormed[isnp] = atof(token);
			}		
			else if (itoken == 8) {
				data->xpehh_normed[isnp] = atof(token);
			}
			else if (itoken == 9) {
				data->lastcol[isnp] = atoi(token);
			} 
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method

void free_xpehh_data(xpehh_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->genpos);
	free(data->freq1);
	free(data->ihh1);
	free(data->freq2); 
	free(data->ihh2); 
	free(data->xpehh_unnormed);
	free(data->xpehh_normed);
	free(data->lastcol);
	data->nsnps = 0;
} //end method

void get_ihs_data(ihs_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->ihh0 = NULL;
	data->ihh1 = NULL;
	data->ihs_unnormed = NULL;
	data->ihs_normed = NULL;
	data->lastcol = NULL; //not sure what information this field contains, selscan documentation is sparse. 0/1

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: "); fprintf(stderr, filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihh0 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh0 != NULL);
	data->ihh1 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh1 != NULL);
	data->ihs_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_unnormed != NULL);
	data->ihs_normed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);

	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 1) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 2) {
				data->freq1[isnp] = atof(token);
			}
			else if (itoken == 3) {
				data->ihh0[isnp] = atof(token);
			}			
			else if (itoken == 4) {
				data->ihh1[isnp] = atof(token);
			}
			else if (itoken == 5) {
				data->ihs_unnormed[isnp] = atof(token);
			}			
			else if (itoken == 6) {
				data->ihs_normed[isnp] = atof(token);
			}
			else if (itoken == 7) {
				data->lastcol[isnp] = atoi(token);
			}
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method

void free_ihs_data(ihs_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->freq1);
	free(data->ihh0);
	free(data->ihh1);
	free(data->ihs_unnormed);
	free(data->ihs_normed);
	free(data->lastcol);	
	data->nsnps = 0;
} //end method

/*************************/
/***SCORE DISTRIBUTIONS***/
/*************************/

void get_likes_data(likes_data* data, char filename[]){
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int ibin, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nbins = 0;
	data->start_bin = NULL; 
	data->end_bin = NULL;
	data->probs = NULL;
	fprintf(stderr, "\tloading likes tables from ");
	fprintf(stderr, filename);
	fprintf(stderr, "\n");

	inf = fopen(filename, "r");
	assert(inf != NULL);
	fgets(newLine, line_size, inf); //strip
	while (fgets(newLine, line_size, inf) != NULL){
		data->nbins++;
	}

	data->start_bin = malloc(data->nbins * sizeof(double));
	data->end_bin = malloc(data->nbins * sizeof(double));
	data->probs = malloc(data->nbins * sizeof(double));
	fclose(inf);

	inf = fopen(filename, "r");
	assert(inf != NULL);
	ibin = 0;
	while (fgets(newLine, line_size, inf) != NULL){
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
			//fprintf(stderr, token);
		 // fprintf(stderr, "\n%d\n", itoken);
			if (itoken == 0) {
				data->start_bin[ibin] = atof(token);
			}
			else if (itoken ==1){
				 data->end_bin[ibin] = atof(token);
			}
			else if (itoken ==2){
				 data->probs[ibin] = atof(token);
				 break;
			} 
		} // end for loop
		ibin++;
	} // end while loop

	fclose(inf);
	free(newLine);
} //end method

void free_likes_data(likes_data* data) {
	if (data == NULL) {return;}
	data->nbins = 0;
	free(data->start_bin);
	free(data->end_bin);
	free(data->probs);
} //end method

/**********************/
/***POP COMPARISONS****/
/**********************/

void get_popComp_data_multiple_region(popComp_data_multiple* data, int argc, char *argv[]){
  ///collates snps for an arbitrary number of populations compared to a putative selpop,
  ///and returns component scores
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilebase[256], infilesuffix[256], infilename[512];
	int isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique; //maxNsnp, i, 
	int startBp, endBp;
	popComp_data data_sing;
	int *allSnps, *allUniqueSnps;
	int iLine;

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

	startBp = atoi(argv[1]);
	endBp = atoi(argv[2]);

	sprintf(infilebase, "/idi/sabeti-scratch/jvitti/synth/cms_pairs/");
	strcat(infilebase, argv[4]); //selpop
	strcat(infilebase, "_");
	sprintf(infilesuffix, "_chr");
	strcat(infilesuffix, argv[3]); //altpop, this gets iterated
	strcat(infilesuffix, "_062215_2pm_nomultialleles_strictMask.cms");



  //////////////////////////
  /// COLLATE LOCI: LOAD ///
  /////////////////////////

	nComparisons = argc - 5;
	totNsnp = 0;
	for (iComp = 1; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 5]);
		strcat(infilename, infilesuffix);
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
	for (iComp = 1; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 5]);
		strcat(infilename, infilesuffix);

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
	for (isnp = 0; isnp < totNsnp-1; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		//else{nunique++;}
		else{
			if (allSnps[isnp] >= startBp && allSnps[isnp] <= endBp){nunique++; iLine = isnp;}
		}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison within region...\n", nunique);

  allUniqueSnps  = malloc(nunique * sizeof(int));
  jsnp = 0;
  for (isnp = 0; isnp < totNsnp-1; isnp++){
    //fprintf(stderr, "%d\t", allSnps[isnp]);
    if (allSnps[isnp] == allSnps[isnp+1]){continue;}
    else{
    	if (allSnps[isnp] >= startBp && allSnps[isnp] <= endBp)
    		{ //fprintf(stderr, " %d ", allSnps[isnp]);
    		allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
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
	for (iComp = 1; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 5]);
		strcat(infilename, infilesuffix);
		//fprintf(stderr, "CHECK THIS ");
		//fprintf(stderr, infilename);
		//fprintf(stderr, "\n");
		get_popComp_data_region(&data_sing, infilename, startBp, endBp);

		fprintf(stderr, "Region starts with SNP at pos: %d\n", data_sing.physpos[0]);

            //isnp iterates (0, nunique) over allUnique Snps 
    jsnp = 0; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
      //fprintf(stderr, "%d\t%d\t%d\t%d\n", isnp, jsnp, allUniqueSnps[isnp], data_sing.physpos[jsnp]);

      if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
        data->physpos[iComp-1][isnp] = data_sing.physpos[jsnp];  
        data->genpos[iComp-1][isnp] = data_sing.genpos[jsnp];   
        data->daf_selpop[iComp-1][isnp] = data_sing.daf_selpop[jsnp];   
        data->delDAF[iComp-1][isnp] = data_sing.delDAF[jsnp];  
        data->fst[iComp-1][isnp] = data_sing.fst[jsnp];   
        data->xp_normed[iComp-1][isnp] = data_sing.xp_normed[jsnp];               
        data->ihs_normed[iComp-1][isnp] = data_sing.ihs_normed[jsnp];   
        data->delihh_normed[iComp-1][isnp] = data_sing.delihh_normed[jsnp];       
        jsnp++; assert(jsnp<=data_sing.nsnps);
      }
      else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; assert(jsnp<=data_sing.nsnps);}
      //else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
    }// end for isnp loop

		free_popComp_data(&data_sing); 
	} // end for icomp
} //end method

void get_popComp_data_region(popComp_data* data, char filename[], int startBp, int endBp){
	// assumes just a single pop comparison
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running, *this_locus_id;
	int isnp, itoken;
	int toTake; //Boolean
	int thisPhysPos;
	float this_genPos, this_daf, this_deldaf, this_xp, this_ihs, this_delihh, this_fst;


	newLine = malloc((line_size+1) * sizeof(char));
	this_locus_id = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->locus_id = NULL;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

	fprintf(stderr, "\tloading from ");
	fprintf(stderr, filename);

	inf = fopen(filename, "r");
	assert(inf != NULL);
	fgets(newLine, line_size, inf); // header
	while (fgets(newLine, line_size, inf) != NULL){
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 1){
					 thisPhysPos = atoi(token);
						if (thisPhysPos >= startBp && thisPhysPos <= endBp)
							{data->nsnps++;}

				} // end if itoken == 1
		} // end for running
	} // end while fgets new line
	fclose(inf);
	fprintf(stderr, "\tnsnp in region: %d\n", data->nsnps);

	data->locus_id = malloc(data->nsnps * sizeof(char*));
	for (isnp = 0; isnp < data->nsnps; isnp++){
		data->locus_id[isnp] = malloc(256*sizeof(char));
		assert(data->locus_id[isnp] != NULL);
	} // end for isnp
	data->physpos = malloc(data->nsnps * sizeof(int));
	data->genpos = malloc(data->nsnps * sizeof(double));
	data->daf_selpop = malloc(data->nsnps * sizeof(double));
	data->delDAF = malloc(data->nsnps * sizeof(double));
	data->fst = malloc(data->nsnps * sizeof(double));
	data->xp_normed = malloc(data->nsnps * sizeof(double));
	data->ihs_normed = malloc(data->nsnps * sizeof(double));
	data->delihh_normed = malloc(data->nsnps * sizeof(double));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);

	inf = fopen(filename, "r");
	fgets(newLine, line_size, inf); // header
	assert(inf != NULL);
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL){
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				toTake = 0;
				if (itoken == 0){strcpy(this_locus_id, token);}//strcpy(data->locus_id[isnp], token);}
				else if (itoken == 1){
					thisPhysPos = atoi(token);
					if (thisPhysPos >= startBp && thisPhysPos <= endBp){toTake = 1;}//fprintf(stderr, "zoopy!\n");}
				} // end else if
				else if (itoken == 2){this_genPos = atof(token);} 
				else if (itoken == 3){this_daf = atof(token);} 
				else if (itoken == 4){this_deldaf = atof(token);} 
				else if (itoken == 5){this_fst = atof(token);} 
				else if (itoken == 6){this_xp = atof(token);} 
				else if (itoken == 7){this_ihs = atof(token);} 
				else if (itoken == 8){this_delihh = atof(token);}
					if (toTake == 1){
						//fprintf(stderr, "whoopy!\n"); 
						strcpy(data->locus_id[isnp], this_locus_id);
						data->physpos[isnp] = thisPhysPos;
						data->genpos[isnp] = this_genPos;
						data->daf_selpop[isnp] = this_daf;
						data->delDAF[isnp] = this_deldaf;
						data->fst[isnp] = this_fst;
						data->xp_normed[isnp] = this_xp;
						//fprintf(stderr, "%f\t", this_xp);
						data->ihs_normed[isnp] = this_ihs;
						data->delihh_normed[isnp] = this_delihh;
						isnp+=1;
					} //end if totake





		


			}  // end for running


	} // end while fgets
		

		fclose(inf);
		free(newLine);
} //end method
void free_popComp_data_multiple_region(popComp_data_multiple* data){
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

void get_popComp_data_multiple(popComp_data_multiple* data, int argc, char *argv[]){
  ///collates snps for an arbitrary number of populations compared to a putative selpop,
  ///and returns component scores
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilebase[256], infilesuffix[256], infilename[512];
	int  isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique; //i,maxNsnp
	popComp_data data_sing;
	int *allSnps, *allUniqueSnps;

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

	sprintf(infilebase, "/idi/sabeti-scratch/jvitti/synth/cms_pairs/");
	strcat(infilebase, argv[2]); //selpop
	strcat(infilebase, "_");
	sprintf(infilesuffix, "_chr");
	strcat(infilesuffix, argv[1]);
	strcat(infilesuffix, "_062215_2pm_nomultialleles_strictMask.cms");

  //////////////////////////
  /// COLLATE LOCI: LOAD ///
  /////////////////////////

	nComparisons = argc - 4;
	totNsnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 6]);
		strcat(infilename, infilesuffix);
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
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 6]);
		strcat(infilename, infilesuffix);
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
	for (isnp = 0; isnp <= totNsnp-1; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{nunique++;}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison...\n", nunique);

  allUniqueSnps  = malloc(nunique * sizeof(int));
  jsnp = 0;
  for (isnp = 0; isnp <= totNsnp-1; isnp++){
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
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 6]);
		strcat(infilename, infilesuffix);
		get_popComp_data(&data_sing, infilename);

            //isnp iterates (0, nunique) over allUnique Snps 
    jsnp = 0; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
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
        jsnp++; assert(jsnp<=data_sing.nsnps);
      }
      else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; assert(jsnp<=data_sing.nsnps);}
      //else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
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

void get_popComp_data(popComp_data* data, char filename[]){
	// assumes just a single pop comparison
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->locus_id = NULL;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

	fprintf(stderr, "\tloading from ");
	fprintf(stderr, filename);

	inf = fopen(filename, "r");
	assert(inf != NULL);
	fgets(newLine, line_size, inf); // header
	while (fgets(newLine, line_size, inf) != NULL){
		data->nsnps++;
	}
	fclose(inf);
	fprintf(stderr, "\tnsnp: %d\n", data->nsnps);

	data->locus_id = malloc(data->nsnps * sizeof(char*));
	for (isnp = 0; isnp < data->nsnps; isnp++){
		data->locus_id[isnp] = malloc(256*sizeof(char));
		assert(data->locus_id[isnp] != NULL);
	}
	data->physpos = malloc(data->nsnps * sizeof(int));
	data->genpos = malloc(data->nsnps * sizeof(double));
	data->daf_selpop = malloc(data->nsnps * sizeof(double));
	data->delDAF = malloc(data->nsnps * sizeof(double));
	data->fst = malloc(data->nsnps * sizeof(double));
	data->xp_normed = malloc(data->nsnps * sizeof(double));
	data->ihs_normed = malloc(data->nsnps * sizeof(double));
	data->delihh_normed = malloc(data->nsnps * sizeof(double));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);

	inf = fopen(filename, "r");
	fgets(newLine, line_size, inf); // header
	assert(inf != NULL);
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL){
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 0) {
					strcpy(data->locus_id[isnp], token); // must check that this works properly
				}
				else if (itoken == 1){
					 data->physpos[isnp] = atoi(token);
				}
				else if (itoken == 2){
					 data->genpos[isnp] = atof(token);
				} 
				else if (itoken == 3){
					 data->daf_selpop[isnp] = atof(token);
				} 
				else if (itoken == 4){
					 data->delDAF[isnp] = atof(token);
				} 
				else if (itoken == 5){
					 data->fst[isnp] = atof(token);
				} 
				else if (itoken == 6){
					 data->xp_normed[isnp] = atof(token);
				} 
				else if (itoken == 7){
					 data->ihs_normed[isnp] = atof(token);
				} 
				else if (itoken == 8){
					 data->delihh_normed[isnp] = atof(token);
					 break;
				} 
			} // end for loop
			isnp++;
		} // end while loop
		fclose(inf);
		free(newLine);
} //end method

void free_popComp_data(popComp_data* data){
	int isnp;
	if (data == NULL) {return;}
	for (isnp = 0; isnp < data->nsnps; isnp++){
		free(data->locus_id[isnp]);
	}
	free(data->locus_id);
	free(data->physpos);
	free(data->genpos);
	free(data->daf_selpop);
	free(data->delDAF);
	free(data->fst);
	free(data->xp_normed);
	free(data->ihs_normed);
	free(data->delihh_normed);
	data->nsnps = 0;
} //end method