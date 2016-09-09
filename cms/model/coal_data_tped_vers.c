// last updated 09.08.16	vitti@broadinstitute.org 

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>
#include "coal_data_tped_vers.h"

/************************/
/***HELPER FUNCTIONS*****/
/************************/

double getGenDist(coal_data* data, int pos_i, int pos_j){
		int pointer_i = (data->nRecom - 1);
		int pointer_j = (data->nRecom - 1);
		int dist, dist_i, dist_j, interval, ibin, jbin;
		double totaldist;
		
		//find correct bin for each pointer
		for (ibin = 0; ibin < data->nRecom; ibin++){
				if (data->physPos[ibin] > pos_i){
						pointer_i = ( ibin - 1 );
						break;
				}
		} // end ibin
		for (jbin = 0; jbin < data->nRecom; jbin++){
				if (data->physPos[jbin] > pos_j){
						pointer_j = ( jbin - 1 );
						break;
				}
		} // end jbin
		
		if (pointer_i == pointer_j){
				dist = (double)(pos_j - pos_i);
				return (dist*data->genPos[pointer_i]);
		} // end if pointer_i
		else{
				dist_i = (data->physPos[pointer_i+1] - pos_i);
				dist_j = (pos_j - data->physPos[pointer_j]);
				totaldist = (double)((dist_i*data->genPos[pointer_i]) + (dist_j*data->genPos[pointer_j])); //bookends
				if ((pointer_i - pointer_j == 1) || (pointer_j - pointer_i == 1)){return totaldist;}
				else
				{pointer_i++;
						while(pointer_i < pointer_j){
						interval = (double)(data->physPos[pointer_i + 1] - data->physPos[pointer_i]);
						totaldist += (interval * data->genPos[pointer_i]);
								pointer_i++;}
				}
				return totaldist;
		} // end else
		return 0.;
} // end function
int getIndexOfItem(int *values, int numVals, int itemToFind){
	int i;
	for (i = 0; i < numVals; i++){
		if (values[i] == itemToFind){
			return i;
		}
	}
	return -1;
} // end function
int getLowerIndexOfItem(int *values, int numVals, int itemToFind){
	/*returns the index of just below the value to be sought*/
	int i;
	for (i = 0; i < numVals; i++){
		if (values[i] == itemToFind){
			return i;
		}
		if (values[i] > itemToFind){
			if (i == 0){return 0;}
			else{return i-1;}
		}
	}
	return numVals-1; //last index
} // end function
int getUpperIndexOfItem(int *values, int numVals, int itemToFind){
	int i;
	for (i = 0; i < numVals; i++){
		if (values[i] == itemToFind){
			return i;
		}
		if (values[i] > itemToFind){
			return i;
		}
	}
	return numVals-1; 
} // end function

/*********************************************/
/***FUNCTIONS FOR HANDLING DATASTRUCTURES*****/
/*********************************************/

void get_coal_data_tped_vers(coal_data* data, char tpedfilename[], char recomfilename[]) {
	const int line_size = 999999999; // upper limit
	const int numRecomLines = 500000;
	//char cmd[600];
	char *newLine, *token, *running;
	int isamp, isnp, itoken, iRecom;
  double genrate;
  FILE *inf=NULL;

	/**************************
	INITIALIZATION, ALLOCATION,
	COUNTING nsample nsnp nrecom
	**************************/
		
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsample = 0;
	data->nsnp = 0;
	data->nRecom = 0;
	data->samp_id = NULL;
	data->pos = NULL;
	data->snp_id = NULL;
	data->chrom = NULL;
	data->genotype = NULL;
	data->anc_base = NULL; 
	data->nallele = NULL;
	data->genPos = NULL;
	data->physPos = NULL;
	data->genloc = NULL;

	//fprintf(stderr, "Getting information from file: %s\n", tpedfilename);

	// Count number of SNPs in tped
	inf = fopen(tpedfilename, "r");

	//handle zipped
	//sprintf(cmd, "gunzip -c %s", tpedfilename);
	//inf = popen(cmd, "r");

	if (inf == NULL) {fprintf(stderr, "Missing TPED file: %s\n", tpedfilename);}
	assert (inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
		assert(strlen(newLine) < line_size);
		data->nsnp++;
	}

	//get number of samples from last line
	int ik = 0;
	char *pch=strchr(newLine,' ');
	while (pch!=NULL) {
		ik++;
		pch=strchr(pch+1,' ');
	}
	data->nsample = ik-3;

	fclose(inf);
	fprintf(stderr, "nSNP: %d\n", data->nsnp);
	fprintf(stderr, "nSample: %d\n", data->nsample);
		
	// Count number of lines in recombination file
  inf = fopen(recomfilename, "r"); assert(inf != NULL);
	fgets(newLine, line_size, inf); //header
	if (inf == NULL) {fprintf(stderr, "Missing recombination file: %s\n", recomfilename);}
	while (fgets(newLine, line_size, inf) != NULL) {
		assert(strlen(newLine) < line_size);
		data->nRecom++;
	}
	fclose(inf);

	// Allocate memory; initialize
	data->samp_id = malloc(data->nsample * sizeof(char*)); assert(data->samp_id != NULL);
	data->genotype = malloc(data->nsample * sizeof(int*)); assert(data->genotype != NULL);
	for (isamp = 0; isamp < data->nsample; isamp++) {
		data->samp_id[isamp] = malloc(64 * sizeof(char)); assert(data->samp_id[isamp] != NULL);
		data->genotype[isamp] = malloc(data->nsnp * sizeof(int)); assert(data->genotype[isamp] !=NULL);
	}
	data->snp_id = malloc(data->nsnp * sizeof(char*)); assert(data->snp_id != NULL);
	for (isnp = 0; isnp < data->nsnp; isnp++) {
		data->snp_id[isnp] = malloc(64 * sizeof(char)); assert(data->snp_id[isnp] != NULL);
	}
	data->pos = malloc(data->nsnp * sizeof(long int)); assert(data->pos != NULL);
	data->chrom = malloc(data->nsnp * sizeof(int)); assert(data->chrom != NULL);
	data->anc_base = malloc(data->nsnp * sizeof(int)); assert(data->anc_base != NULL);
	data->snp_base[0] = malloc(data->nsnp * sizeof(char*)); assert(data->snp_base[0] != NULL);
	data->snp_base[1] = malloc(data->nsnp * sizeof(char*)); assert(data->snp_base[1] != NULL);
	for (isnp = 0; isnp < data->nsnp; isnp++) {
		data->snp_base[0][isnp] = malloc(3 * sizeof(char)); assert(data->snp_base[0][isnp] != NULL);
		data->snp_base[1][isnp] = malloc(3 * sizeof(char)); assert(data->snp_base[1][isnp] != NULL);
	}
	data->genloc = malloc(data->nsnp * sizeof(double)); assert(data->genloc != NULL);
	data->snp_base[2] = calloc(data->nsnp, sizeof(char*)); // STRIKE?
	data->snp_base[3] = calloc(data->nsnp, sizeof(char*)); // STRIKE?
	data->nallele = malloc(data->nsnp * sizeof(int)); assert(data->nallele != NULL);
	data->genPos = malloc(numRecomLines * sizeof(double)); assert(data->genPos !=NULL);
	data->physPos = malloc(numRecomLines * sizeof(int)); assert(data->physPos != NULL);

	/*******************
	GET DATA FROM TPED
	*******************/
	//handle zipped
	//sprintf(cmd, "gunzip -c %s", tpedfilename);
	///inf = popen(cmd, "r");
	inf = fopen(tpedfilename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing TPED file: %s\n", tpedfilename);}
	assert(inf != NULL);

	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 0) {
				//data->chrom[isnp] = chromosome;
				data->nallele[isnp] = 2; //only biallelics in tped
			}	
			else if (itoken == 1) {
					strcpy(data->snp_id[isnp], token);
			}
			else if (itoken == 2){
				data->genloc[isnp] = atof(token);
			}

			else if (itoken == 3) {
				data->pos[isnp] = atoi(token); 
			}
			else if (itoken >= 4) {
				isamp = itoken - 4; //RIGHT?
				data->genotype[isamp][isnp] = atoi(token);
			}
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	fclose(inf);

	/**************************
	GET DATA FROM RECOMB FILE
	**************************/
	inf = fopen(recomfilename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing recombination file: %s\n", recomfilename);}
	assert(inf != NULL);
	fgets(newLine, line_size, inf); //header
	iRecom = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
			if (itoken == 1) {
				data->physPos[iRecom] = atoi(token);
			}
			else if (itoken == 2) {
					genrate = atof(token);
				data->genPos[iRecom] = (genrate / 1000000.); //per Mb rate
				iRecom++;
				break;
			}
		}
	} // end while fgets 
	fclose(inf);
	free(newLine);
} // end function


void free_coal_data(coal_data* data) {
	int isamp, isnp;
	if (data == NULL) {return;}
	if (data->samp_id == NULL) {return;}
	for (isamp = 0; isamp < data->nsample; isamp++) {
		free(data->genotype[isamp]);
		free(data->samp_id[isamp]);
	}
	for (isnp = 0; isnp < data->nsnp; isnp++) {
		free(data->snp_id[isnp]);
	}
	free(data->pos);
	free(data->chrom);
	free(data->snp_id);
	free(data->anc_base);
	free(data->snp_base[0]);
	free(data->snp_base[1]);
	free(data->snp_base[2]);
	free(data->snp_base[3]);
	free(data->nallele);
	free(data->samp_id);
	free(data->genotype);
	free(data->genPos);
	free(data->physPos);
	free(data->genloc);
	data->nsnp = 0;
	data->nsample = 0;
	data->nRecom = 0;
} // end function