// functions for handling cms component(+composite) score datastructures
// 12.21.2017	agnostic to un/zipped 		

/// ADD LIKES 

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include "cms_data_vers.h"

/*********************/
/***HELPER METHODS****/ 
/*********************/
int intcmp(const void *v1, const void *v2) {return (*(int *)v1 - *(int *)v2);}
int count_unique_from_sorted(int arr[], int len) {
	// given an array of sorted values; returns n_transitions + 1
    int j;
    int count = 0;
    int lastPos = -1;

    for (j = 0; j < len; j++){
    	if (arr[j] != lastPos){count +=1;}
    	lastPos = arr[j];
    }
    count++;
    return count;
} //end method

/**********************/
/***COMPONENT SCORES***/
/**********************/
void get_freqs_data(freqs_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	int thisPos;
	float thisGenPos, thisPopDaf, thisDelDaf, thisFst;
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->fst = NULL;
	data->deldaf = NULL;
	data->popdaf = NULL;
	data->genpos = NULL;

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	/////////////////
	/// GET NSNPS ///
	/////////////////	
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		gzgets(zinf, newLine, line_size); // strip header
		while (gzgets(zinf, newLine, line_size) != NULL) {
				assert(strlen(newLine) < line_size);
				for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
						if (itoken == 0) {thisPos = atoi(token);}	
						if (itoken == 1) {break;}	
				} //end for running
				if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
				if (thisPos > maxPos){break;}	
		} // end while gets line
		gzclose(zinf);
	} // end zipped
	if (zipped  == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // strip header
		while (fgets(newLine, line_size, inf) != NULL) {
				assert(strlen(newLine) < line_size);
				for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
						if (itoken == 0) {thisPos = atoi(token);}	
						if (itoken == 1) {break;}	
				} //end for running
				if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
				if (thisPos > maxPos){break;}	
		} // end while gets line
		fclose(inf);
	} //end not-zipped

	///////////////////////
	/// ALLOCATE MEMORY ///
	///////////////////////
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->fst = malloc(data->nsnps * sizeof(double*)); assert(data->fst != NULL);
	data->deldaf = malloc(data->nsnps * sizeof(double*)); assert(data->deldaf != NULL);
	data->genpos = malloc(data->nsnps * sizeof(double*)); assert(data->genpos != NULL);
	data->popdaf = malloc(data->nsnps * sizeof(double*)); assert(data->popdaf != NULL);

	////////////////////////
	// GET DATA FROM FILE //
	////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		gzgets(zinf, newLine, line_size); // strip header
		isnp = 0;
		thisGenPos = thisPopDaf = thisDelDaf = thisFst = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}	
				else if (itoken == 1) {thisGenPos = atof(token);}
				else if (itoken == 2) {thisPopDaf = atof(token);}	
				else if (itoken == 3) {thisDelDaf = atof(token);}
				else if (itoken == 4) {thisFst = atof(token);}			
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->genpos[isnp] = thisGenPos;
				data->popdaf[isnp] = thisPopDaf;		
				data->deldaf[isnp] = thisDelDaf;
				data->fst[isnp] = thisFst;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))
		gzclose(zinf);
	} // end zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		fgets(newLine, line_size, inf); // strip header
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}	
				else if (itoken == 1) {thisGenPos = atof(token);}
				else if (itoken == 2) {thisPopDaf = atof(token);}	
				else if (itoken == 3) {thisDelDaf = atof(token);}
				else if (itoken == 4) {thisFst = atof(token);}			
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->genpos[isnp] = thisGenPos;
				data->popdaf[isnp] = thisPopDaf;		
				data->deldaf[isnp] = thisDelDaf;
				data->fst[isnp] = thisFst;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))
	fclose(inf);
	} // end unzipped
	free(newLine);
} //end method
void free_freqs_data(freqs_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->fst);
	free(data->deldaf);
	free(data->genpos);
	free(data->popdaf);	
	data->nsnps = 0;
} //end method
void get_delihh_data(delihh_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000;	
	FILE *inf=NULL;
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	int thisPos, thisLastcol;
	float thisFreq, thisIhs, thisUnnormedDelihh, thisNormedDelihh;	
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->ihs_unnormed = NULL;
	data->delihh_unnormed = NULL;
	data->delihh_normed = NULL;
	data->lastcol = NULL;

	/////////////////
	/// GET NSNPS ///
	/////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				if (itoken == 2) {break;}	
			} //end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
			if (thisPos > maxPos){break;}	
		} // end while
		gzclose(zinf);
	} // end zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				if (itoken == 2) {break;}	
			} //end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
			if (thisPos > maxPos){break;}	
		} // end while
		fclose(inf);
	} //end unzipped

	///////////////////////
	/// ALLOCATE MEMORY ///
	///////////////////////
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihs_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_unnormed != NULL);
	data->delihh_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->delihh_unnormed != NULL);	
	data->delihh_normed = malloc(data->nsnps * sizeof(double*)); assert(data->delihh_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);
	
	///////////////////////////
	/// LOAD DATA FROM FILE ///
	///////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisIhs = atof(token);}
				else if (itoken == 4) {thisUnnormedDelihh = atof(token);}
				//5 is duplicate!!!!
				else if (itoken == 6) {thisNormedDelihh = atof(token);}
				else if (itoken == 7) {thisLastcol = atoi(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->ihs_unnormed[isnp] = thisIhs;
				data->delihh_unnormed[isnp] = thisUnnormedDelihh;
				data->delihh_normed[isnp] = thisNormedDelihh;
				data->lastcol[isnp] = thisLastcol;			
				isnp++;
			} // end if in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))	
		gzclose(zinf);
	} // end zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisIhs = atof(token);}
				else if (itoken == 4) {thisUnnormedDelihh = atof(token);}
				//5 is duplicate!!!!
				else if (itoken == 6) {thisNormedDelihh = atof(token);}
				else if (itoken == 7) {thisLastcol = atoi(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->ihs_unnormed[isnp] = thisIhs;
				data->delihh_unnormed[isnp] = thisUnnormedDelihh;
				data->delihh_normed[isnp] = thisNormedDelihh;
				data->lastcol[isnp] = thisLastcol;			
				isnp++;
			} // end if in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))	
		fclose(inf);
	} // end zipped
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
void get_xpehh_data(xpehh_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int	isnp, itoken;
	int thisPos, thisLastCol;
	float thisGenPos, thisFreq1, thisIhh1, thisFreq2, thisIhh2, thisXpUn, thisXpNormed;	
	int zipped; // Boolean 0T 1F

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

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	//////////////////////
	// ALLOCATE MEMORY ///
	//////////////////////	
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		gzgets(zinf, newLine, line_size); //strip header
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			//fprintf(stderr, newLine);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
					if (itoken == 1) {thisPos = atoi(token);}	
					if (itoken == 2) {break;}	
			} //end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
			if (thisPos > maxPos){break;}	
		}//end while gets
		gzclose(zinf);
	} // end zipped
	if (zipped == 1){
		//zinf = open(filename, "r");
		//if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		inf = fopen(filename, "r");
		fgets(newLine, line_size, inf);
		assert(inf != NULL);
		//gzgets(zinf, newLine, line_size); //strip header
		while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			//fprintf(stderr, newLine);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
					if (itoken == 1) {thisPos = atoi(token);}	
					if (itoken == 2) {break;}	
			} //end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}	
			if (thisPos > maxPos){break;}	
		}//end while gets
		//gzclose(zinf);
		fclose(inf);
	} // end zipped

	//////////////////////
	// ALLOCATE MEMORY ///
	//////////////////////	
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->genpos = malloc(data->nsnps * sizeof(double*)); assert(data->genpos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihh1 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh1 != NULL);
	data->freq2 = malloc(data->nsnps * sizeof(double*)); assert(data->freq2 != NULL);
	data->ihh2 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh2 != NULL);
	data->xpehh_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->xpehh_unnormed != NULL);
	data->xpehh_normed = malloc(data->nsnps * sizeof(double*)); assert(data->xpehh_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);

	//////////////////////////
	// LOAD DATA FROM FILE ///
	//////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		gzgets(zinf, newLine, line_size); // strip header
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				else if (itoken == 2) {thisGenPos = atof(token);}
				else if (itoken == 3) {thisFreq1 = atof(token);}		
				else if (itoken == 4) {thisIhh1 = atof(token);}
				else if (itoken == 5) {thisFreq2 = atof(token);}
				else if (itoken == 6) {thisIhh2 = atof(token);}
				else if (itoken == 7) {thisXpUn = atof(token);}
				else if (itoken == 8) {thisXpNormed =  atof(token);}
				else if (itoken == 9) {thisLastCol = atoi(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->genpos[isnp] = thisGenPos;
				data->freq1[isnp] = thisFreq1;
				data->ihh1[isnp] = thisIhh1;
				data->freq2[isnp] = thisFreq2;
				data->ihh2[isnp] = thisIhh2;
				data->xpehh_unnormed[isnp] = thisXpUn; 
				data->xpehh_normed[isnp] = thisXpNormed;
				data->lastcol[isnp] = thisLastCol;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))
		gzclose(zinf);
	} // end zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		fgets(newLine, line_size, inf); // strip header
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				else if (itoken == 2) {thisGenPos = atof(token);}
				else if (itoken == 3) {thisFreq1 = atof(token);}		
				else if (itoken == 4) {thisIhh1 = atof(token);}
				else if (itoken == 5) {thisFreq2 = atof(token);}
				else if (itoken == 6) {thisIhh2 = atof(token);}
				else if (itoken == 7) {thisXpUn = atof(token);}
				else if (itoken == 8) {thisXpNormed =  atof(token);}
				else if (itoken == 9) {thisLastCol = atoi(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->genpos[isnp] = thisGenPos;
				data->freq1[isnp] = thisFreq1;
				data->ihh1[isnp] = thisIhh1;
				data->freq2[isnp] = thisFreq2;
				data->ihh2[isnp] = thisIhh2;
				data->xpehh_unnormed[isnp] = thisXpUn; 
				data->xpehh_normed[isnp] = thisXpNormed;
				data->lastcol[isnp] = thisLastCol;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))	
		fclose(inf);
	} // end if unzipped
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
void get_ihs_data(ihs_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	int unnormed_index, normed_index; // selscan output may or may not include ihh decomposition; account for both
	int thisPos; 
	float thisFreq, thisIhh0, thisIhh1, thisUnnormed, thisNormed;
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->ihh0 = NULL;
	data->ihh1 = NULL;
	data->ihs_unnormed = NULL;
	data->ihs_normed = NULL;
	data->lastcol = NULL; 

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	/////////////////
	// COUNT SNPS ///
	/////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				if (itoken == 2) {break;}
			} // end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos){break;}
		} // end while loop
		gzclose(zinf);
	} // end if zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				if (itoken == 2) {break;}
			} // end for running
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos){break;}
		} // end while loop
		fclose(inf);
	} // end if unzipped

	// determine formatting based on line (selscan output may or may not include ihh decomposition)
	// https://github.com/szpiech/selscan
	unnormed_index = 0; //5 9
	normed_index = 0; //6 10
	int dummy_var = 0;
	//Take the first line and check its length
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		gzgets(zinf, newLine, line_size);
		gzclose(zinf);
	} // end if zipped	
	if (zipped == 1){
		inf = fopen(filename, "r");
		fgets(newLine, line_size, inf);
		fclose(inf);
	} // end if unzipped


	//fprintf(stderr, newLine);
	for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++)
		{dummy_var ++;}
	//fprintf(stderr, "dummy_var: %d\n", dummy_var);
	if (dummy_var <= 8){unnormed_index = 5; normed_index = 6;}
	else {unnormed_index = 9; normed_index = 10;}

	/////////////////////
	// ALLOCATE MEMORY //
	/////////////////////	
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->ihh0 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh0 != NULL);
	data->ihh1 = malloc(data->nsnps * sizeof(double*)); assert(data->ihh1 != NULL);
	data->ihs_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_unnormed != NULL);
	data->ihs_normed = malloc(data->nsnps * sizeof(double*)); assert(data->ihs_normed != NULL);
	data->lastcol = malloc(data->nsnps * sizeof(int*)); assert(data->lastcol != NULL);

	////////////////////////
	// GET DATA FROM FILE //
	////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	 
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisIhh0= atof(token);}			
				else if (itoken == 4) {thisIhh1 = atof(token);}
				else if (itoken == unnormed_index){thisUnnormed = atof(token);}			
				else if (itoken == normed_index){thisNormed = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->ihh0[isnp] = thisIhh0;
				data->ihh1[isnp] = thisIhh1;
				data->ihs_unnormed[isnp] = thisUnnormed;
				data->ihs_normed[isnp] = thisNormed;
				isnp++;
			}// end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))		
		gzclose(zinf);
	} // end if zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	 
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisIhh0= atof(token);}			
				else if (itoken == 4) {thisIhh1 = atof(token);}
				else if (itoken == unnormed_index){thisUnnormed = atof(token);}			
				else if (itoken == normed_index){thisNormed = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->ihh0[isnp] = thisIhh0;
				data->ihh1[isnp] = thisIhh1;
				data->ihs_unnormed[isnp] = thisUnnormed;
				data->ihs_normed[isnp] = thisNormed;
				isnp++;
			}// end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))		
		fclose(inf);
	} // end if unzipped
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
void get_nsl_data(nsl_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;	
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	int thisPos;
	float thisFreq, thisSl0, thisSl1, thisUnnormed, thisNormed;
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->sl0 = NULL;
	data->sl1 = NULL;
	data->nsl_unnormed = NULL;
	data->nsl_normed = NULL;

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	///////////////
	// GET NSNPS //
	///////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				if (itoken == 2) {break;}
			} // end for-running	
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos){break;}
		} // end while loop
		gzclose(zinf);		
	} // end zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}
				if (itoken == 2) {break;}
			} // end for-running	
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos){break;}
		} // end while loop
		fclose(inf);
	} // end unzipped

	/////////////////////
	// ALLOCATE MEMORY //
	/////////////////////	
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->sl0 = malloc(data->nsnps * sizeof(double*)); assert(data->sl0 != NULL);
	data->sl1 = malloc(data->nsnps * sizeof(double*)); assert(data->sl1 != NULL);
	data->nsl_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->nsl_unnormed != NULL);
	data->nsl_normed = malloc(data->nsnps * sizeof(double*)); assert(data->nsl_normed != NULL);

	////////////////////////
	// GET DATA FROM FILE //
	////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisSl0 = atof(token);}			
				else if (itoken == 4) {thisSl1 = atof(token);}
				else if (itoken == 5) {thisUnnormed = atof(token);}			
				else if (itoken == 6) {thisNormed = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->sl0[isnp] = thisSl0;
				data->sl1[isnp] = thisSl1;
				data->nsl_unnormed[isnp] = thisUnnormed;
				data->nsl_normed[isnp] = thisNormed;
				isnp++;
			} // end if-in-region
			if (thisPos >= maxPos){break;}
		} //END while(fgets(newLine))
		gzclose(zinf);		
	}// end if zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 1) {thisPos = atoi(token);}	
				else if (itoken == 2) {thisFreq = atof(token);}
				else if (itoken == 3) {thisSl0 = atof(token);}			
				else if (itoken == 4) {thisSl1 = atof(token);}
				else if (itoken == 5) {thisUnnormed = atof(token);}			
				else if (itoken == 6) {thisNormed = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->freq1[isnp] = thisFreq;
				data->sl0[isnp] = thisSl0;
				data->sl1[isnp] = thisSl1;
				data->nsl_unnormed[isnp] = thisUnnormed;
				data->nsl_normed[isnp] = thisNormed;
				isnp++;
			} // end if-in-region
			if (thisPos >= maxPos){break;}
		} //END while(fgets(newLine))
		fclose(inf);
	} // end if unzipped
	free(newLine);
} //end method
void free_nsl_data(nsl_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->freq1);
	free(data->sl0);
	free(data->sl1);
	free(data->nsl_unnormed);
	free(data->nsl_normed);
	data->nsnps = 0;
} //end method
void get_H12_data(H12_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;	
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	float thisH12, thisH2H1;
	int thisPos;
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->H12_value = NULL;
	data->H2H1_value = NULL;

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	/////////////////
	// COUNT NSNPS //
	/////////////////
	if (zipped == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}	
				if (itoken == 1) {break;}
			}// end for line
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos) {break;}
		} // end while
		fclose(inf);
	} // end unzipped
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}	
				if (itoken == 1) {break;}
			}// end for line
			if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			if (thisPos > maxPos) {break;}
		} // end while
		gzclose(zinf);
	} // end zipped

	/////////////////////
	// ALLOCATE MEMORY //
	/////////////////////	
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->H12_value = malloc(data->nsnps * sizeof(double*)); assert(data->H12_value != NULL);
	data->H2H1_value = malloc(data->nsnps * sizeof(double*)); assert(data->H2H1_value != NULL);

	////////////////////////
	// GET DATA FROM FILE //
	////////////////////////
	if (zipped == 1){
		inf = fopen(filename, "r");
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos= atoi(token);}	
				else if (itoken == 8) {thisH12 = atof(token);}
				else if (itoken == 9) {thisH2H1 = atof(token);}			
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->H12_value[isnp] = thisH12;
				data->H2H1_value[isnp] = thisH2H1;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))	
		fclose(inf);
	} // end unzipped
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos= atoi(token);}	
				else if (itoken == 1) {thisH12 = atof(token);}
				else if (itoken == 2) {thisH2H1 = atof(token);}			

				//else if (itoken == 8) {thisH12 = atof(token);}
				//else if (itoken == 9) {thisH2H1 = atof(token);}			
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->H12_value[isnp] = thisH12;
				data->H2H1_value[isnp] = thisH2H1;
				isnp++;
			} // end if-in-region
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))		
		gzclose(zinf);
	} // end zipped
	free(newLine);
} //end method
void free_H12_data(H12_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->H12_value);
	free(data->H2H1_value);
	data->nsnps = 0;
} //end method
void get_iSAFE_data(iSAFE_data* data, char filename[], int minPos, int maxPos) {
	const int line_size = 15000000; 
	FILE *inf=NULL;	
	gzFile zinf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;
	int thisPos;
	float thisISAFE;
	int zipped; // Boolean 0T 1F

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->iSAFE_value = NULL;

	if (strstr(filename, ".gz") == NULL){zipped = 1;} 	// unzipped file
	else {zipped = 0;} //zipped

	/////////////////
	// COUNT SNPS ///
	/////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		if (zinf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(zinf != NULL);
		gzgets(zinf, newLine, line_size); //strip header
		while (gzgets(zinf, newLine, line_size) != NULL) {
			assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}
				if (itoken == 1) {break;}	
			} //end for loop
			 if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			 if (thisPos > maxPos) {break;}
		} // end while loop
		gzclose(zinf);
	} // end if zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
		assert(inf != NULL);
		fgets(newLine, line_size, inf); //strip header
		while (fgets(newLine, line_size, inf) != NULL) {
				assert(strlen(newLine) < line_size);
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}
				if (itoken == 1) {break;}	
			} //end for loop
			 if (thisPos >= minPos && thisPos <= maxPos){data->nsnps++;}
			 if (thisPos > maxPos) {break;}
		} // end while loop
		fclose(inf);
	} // end if unzipped

	/////////////////////
	// ALLOCATE MEMORY //
	/////////////////////
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->iSAFE_value = malloc(data->nsnps * sizeof(double*)); assert(data->iSAFE_value != NULL);

	////////////////////////
	// GET DATA FROM FILE //
	////////////////////////
	if (zipped == 0){
		zinf = gzopen(filename, "rb");
		gzgets(zinf, newLine, line_size); //strip header
		isnp = 0;
		while (gzgets(zinf, newLine, line_size) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
				if (itoken == 0) {thisPos = atoi(token);}	
				else if (itoken == 1) {thisISAFE = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->iSAFE_value[isnp] = thisISAFE;
				isnp++;
			}
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))
		gzclose(zinf);
	} // end if zipped
	if (zipped == 1){
		inf = fopen(filename, "r");
		fgets(newLine, line_size, inf); //strip header
		isnp = 0;
		while (fgets(newLine, line_size, inf) != NULL) {
			for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
					if (itoken == 0) {thisPos = atoi(token);}	
				else if (itoken == 1) {thisISAFE = atof(token);}
			} // END for running=newLine
			if (thisPos >= minPos && thisPos <= maxPos){
				data->pos[isnp] = thisPos;
				data->iSAFE_value[isnp] = thisISAFE;
				isnp++;
			}
			if (thisPos > maxPos){break;}
		} //END while(fgets(newLine))
		fclose(inf);
	} // end if unzipped
	free(newLine);
} //end method
void free_iSAFE_data(iSAFE_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->iSAFE_value);
	data->nsnps = 0;
} //end method

/****************/
/***POP PAIR****/ 
/***************/
int get_num_anyData(int minPos, int maxPos, char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]){
	/* counts the number of unique SNPs that have a score value for at least one component score */
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;
	int nunique;
	int totNsnp, isnp, jsnp;
	int *allSnps;
	int nFiltered_iHS, nFiltered_delihh, nFiltered_nsl;
	int nFiltered_h12, nFiltered_isafe, nFiltered_freqs;
	int nFiltered_xp;

	////////////////////////////////////
	// LOAD IN DATA; FILTER TO REGION //
	////////////////////////////////////
	get_ihs_data(&ihs1, ihs_filename, minPos, maxPos);
	get_delihh_data(&delihh1, delihh_filename, minPos, maxPos);
	get_nsl_data(&nsl1, nsl_filename, minPos, maxPos);
	get_xpehh_data(&xp, xpehh_filename, minPos, maxPos);
	get_freqs_data(&freqs, freqs_filename, minPos, maxPos);
	get_H12_data(&H12, H12_filename, minPos, maxPos);
	get_iSAFE_data(&iSAFE, iSAFE_filename, minPos, maxPos);

	nFiltered_iHS = ihs1.nsnps;
	nFiltered_delihh = delihh1.nsnps;
	nFiltered_nsl = nsl1.nsnps;
	nFiltered_h12 = H12.nsnps;
	nFiltered_isafe = iSAFE.nsnps;
	nFiltered_freqs = freqs.nsnps;
	nFiltered_xp = xp.nsnps;
	totNsnp = nFiltered_iHS + nFiltered_delihh + nFiltered_nsl + nFiltered_h12 + nFiltered_isafe + nFiltered_freqs + nFiltered_xp;
	allSnps = calloc(totNsnp, sizeof(int));
	
	//////////////////////////////////////////////
	// LOAD REDUNDANT POSITIONS FROM ALL SCORES //
	//////////////////////////////////////////////
	jsnp = -1;
	for (isnp = 0; isnp < ihs1.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = ihs1.pos[isnp];
	}	// end ihs isnp
	for (isnp = 0; isnp < delihh1.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = delihh1.pos[isnp];
	}	// end delihh isnp
	for (isnp = 0; isnp < nsl1.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = nsl1.pos[isnp];
	}	// end nsl isnp
	for (isnp = 0; isnp < xp.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = xp.pos[isnp];
	}	// end xp isnp
	for (isnp = 0; isnp < freqs.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = freqs.pos[isnp];
	}	// end freqs isnp
	for (isnp = 0; isnp < H12.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = H12.pos[isnp];
	}	//end h12 isnp
	for (isnp = 0; isnp < iSAFE.nsnps; isnp++){
		jsnp++;
		allSnps[jsnp] = iSAFE.pos[isnp];
	}	// end isafe isnp

	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	nunique = count_unique_from_sorted(allSnps, totNsnp);

	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
	free(allSnps);
	return nunique;
} // end function
void get_popPair_anyData(int minPos, int maxPos, popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[], int preCountedSNPs){
	/*	loads in ANY data for a population pair  -- this implementation assumes all SNPs are present in freqs file*/
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;
	int nunique;
	int isnp, jsnp;
	int ksnp, lsnp, msnp, osnp, psnp, qsnp;

	/////////////////
	// INITIALIZE ///
	/////////////////
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
	data->nsl_normed = NULL;
	data->H12 = NULL;
	data->H2H1 = NULL;
	data->iSAFE = NULL;

	////////////////////////////////////////
	// LOAD COMPONENT SCORES BY POSITION ///
	////////////////////////////////////////
	get_ihs_data(&ihs1, ihs_filename, minPos, maxPos);
	get_delihh_data(&delihh1, delihh_filename, minPos, maxPos);
	get_nsl_data(&nsl1, nsl_filename, minPos, maxPos);
	get_xpehh_data(&xp, xpehh_filename, minPos, maxPos);
	get_freqs_data(&freqs, freqs_filename, minPos, maxPos);
	get_H12_data(&H12, H12_filename, minPos, maxPos);
	get_iSAFE_data( &iSAFE, iSAFE_filename, minPos, maxPos);
	
	nunique = freqs.nsnps;

	/////////////////////
	// RESERVE MEMORY ///
	/////////////////////
	data->nsnps = nunique;	//get_num_anyData(minPos, maxPos, ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename); 
	data->locus_id = malloc(data->nsnps * sizeof(char*));
	for (isnp = 0; isnp < data->nsnps; isnp++){
		data->locus_id[isnp] = malloc(256*sizeof(char));
		assert(data->locus_id[isnp] != NULL);
	}
	data->physpos = calloc(data->nsnps, sizeof(int));
	data->genpos = calloc(data->nsnps, sizeof(double)); 
	data->daf_selpop = calloc(data->nsnps, sizeof(double));
	data->delDAF = calloc(data->nsnps, sizeof(double));
	data->fst = calloc(data->nsnps, sizeof(double));
	data->xp_normed = calloc(data->nsnps, sizeof(double));
	data->ihs_normed = calloc(data->nsnps, sizeof(double));
	data->delihh_normed = calloc(data->nsnps, sizeof(double));
	data->nsl_normed = calloc(data->nsnps, sizeof(double)); 
	data->H12 = calloc(data->nsnps, sizeof(double));
	data->H2H1 = calloc(data->nsnps, sizeof(double));
	data->iSAFE = calloc(data->nsnps, sizeof(double)); 
	
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);
	assert(data->nsl_normed != NULL);
	assert(data->H12 != NULL);
	assert(data->H2H1 != NULL);	
	assert(data->iSAFE != NULL);

	for (isnp = 0; isnp < nunique; ++isnp){
		data->physpos[isnp] = freqs.pos[isnp];
	} //sone nan



	//DEBUG		---> HAVE THINGS GONE WRONG BY HERE? NON NAN VALS! 
	//fprintf(stderr, "\t%d\n", xp.nsnps);
	//for (isnp = 0; isnp < xp.nsnps; ++isnp){
	//	fprintf(stderr, "%f\n", xp.xpehh_normed[isnp]);
	//}

	//////////////////////////////////////////////////
	//filter non-redundant positions to data object // 
	//////////////////////////////////////////////////	

	jsnp = 0;	//ihs
	ksnp = 0;	//delihh
	lsnp = 0;	//nsl
	msnp = 0;	//xpehh	
	osnp = 0;	//freqs	
	psnp = 0; //H12
	qsnp = 0;	//iSAFE
	//if this database_snp is the one I'm looking at, record it and advance pointer.
	//advance pointer IF below all SNPs /AND/ this physpos has already been passed--!

	for (isnp = 0; isnp < data->nsnps; isnp++){
		//ihs (jsnp)
		while (data->physpos[isnp] > ihs1.pos[jsnp]) {jsnp++;}
		if(data->physpos[isnp] == ihs1.pos[jsnp]){data->ihs_normed[isnp] = ihs1.ihs_normed[jsnp]; if(jsnp < ihs1.nsnps - 1){jsnp++;}} 
		else{data->ihs_normed[isnp] = 0.0 / 0.0;} //NaN
		//if (data->physpos[isnp] > ihs1.pos[jsnp]) {jsnp++;}

		//delihh (ksnp)
		while (data->physpos[isnp] > delihh1.pos[ksnp]) {ksnp++;}
		if(data->physpos[isnp] == delihh1.pos[ksnp]){data->delihh_normed[isnp] = delihh1.delihh_normed[ksnp]; if(ksnp < delihh1.nsnps - 1){ksnp++;}}
		else{data->delihh_normed[isnp] = 0.0 / 0.0;} //NaN
		//if (data->physpos[isnp] > delihh1.pos[ksnp]) {ksnp++;}
	
		//nsl (lsnp)
		while (data->physpos[isnp] > nsl1.pos[lsnp]) {lsnp++;}
		if(data->physpos[isnp] == nsl1.pos[lsnp]){data->nsl_normed[isnp] = nsl1.nsl_normed[lsnp];if(lsnp < nsl1.nsnps -1){lsnp++;}}
		else{data->nsl_normed[isnp] = 0.0 / 0.0;} //NaN
		//if (data->physpos[isnp] > nsl1.pos[lsnp]) {lsnp++;}

		//xp-ehh (msnp)
		while (data->physpos[isnp] > xp.pos[msnp]) {msnp++;} //POSSIBLE THESE ARE NON-STRINCTLY MONOTONIC NOW?? no...
		if(data->physpos[isnp] == xp.pos[msnp]){data->xp_normed[isnp] = xp.xpehh_normed[msnp]; if(msnp < xp.nsnps - 1){msnp++;}} //;fprintf(stderr,"%f\n", xp.xpehh_normed[msnp]
		else{data->xp_normed[isnp] = 0.0 / 0.0;} //NaN
		//while (data->physpos[isnp] > xp.pos[msnp]) {msnp++;} //POSSIBLE THESE ARE NON-STRINCTLY MONOTONIC NOW?? no...
		//ELIF???
		//what if msnp ADVANCES (because we had a hit) but the next one SHOULD also be a hit, and we advance prematurely! (msnp++ twice on one round of isnp)


		//freqs (osnp)
		while (data->physpos[isnp] > freqs.pos[osnp]) {osnp++;}
		if(data->physpos[isnp] == freqs.pos[osnp]){data->fst[isnp] = freqs.fst[osnp]; data->delDAF[isnp] = freqs.deldaf[osnp]; data->genpos[isnp] = freqs.genpos[osnp]; data->daf_selpop[isnp] = freqs.popdaf[osnp]; if(osnp < freqs.nsnps - 1){osnp++;}}
		else{data->fst[isnp] = 0.0 / 0.0; data->delDAF[isnp] = 0.0 / 0.0; } //NaN
		//if (data->physpos[isnp] > freqs.pos[osnp]) {osnp++;}

		//h12 (psnp)
		while (data->physpos[isnp] > H12.pos[psnp]) {psnp++;}
		if(data->physpos[isnp] == H12.pos[psnp]){data->H12[isnp] = H12.H12_value[psnp]; data->H2H1[isnp] = H12.H2H1_value[psnp]; if(psnp < H12.nsnps - 1){psnp++;}}
		else{data->H12[isnp] = 0.0 / 0.0; data->H2H1[isnp] = 0.0 / 0.0;} //NaN
		//if (data->physpos[isnp] > H12.pos[psnp]) {psnp++;}

		//iSAFE (qsnp)
		while (data->physpos[isnp] > iSAFE.pos[qsnp]) {qsnp++;}
		if(data->physpos[isnp] == iSAFE.pos[qsnp]){ data->iSAFE[isnp] = iSAFE.iSAFE_value[qsnp]; if(qsnp < iSAFE.nsnps - 1){qsnp++;}}
		else{data->iSAFE[isnp] = 0.0 / 0.0;} //NaN
		//if (data->physpos[isnp] > iSAFE.pos[qsnp]) {qsnp++;}
	} // end isnp


	//DEBUG		---> THINGS HAVE GOTNE WRONG BY HERE!!
	//for (isnp = 0; isnp < nunique; ++isnp){
	//	fprintf(stderr, "%f\n", data->xp_normed[isnp]);
	//}



	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
} //end method
void free_popPair_data(popPair_data* data){
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
	free(data->nsl_normed);
	free(data->H12);
	free(data->H2H1);
	free(data->iSAFE);
	data->nsnps = 0;
} //end method

/***********************/
/***POP COMPARISONS*****/	
/***********************/
void get_popComp_anyData(int minPos, int maxPos, popComp_data_multiple* data, int nComparisons, int argc, char *argv[]){
	/*This function loads in all SNPs, even those for which there is incomplete data.
	argv is all files for pop-pairs (each of which points to further component score files)*/
	const int line_size = 15000000; 
	popPair_data data_sing;
	gzFile zinf=NULL;
	char *newLine;//, *token, *running;
	char infilename[512];
	int	isnp, jsnp, iComp, nunique; //thisPhysPos, itoken, totNsnp
	int *allUniqueSnps; //*allSnps, 
	//int *compCountSnps;
	int thisCompSnpCount;
	char ihs_filename[528], delihh_filename[528], nsl_filename[528], xpehh_filename[528], freqs_filename[528];
	char H12_filename[528], iSAFE_filename[528];
	int numExtraArgs;
	//int checked, flag, n;
	freqs_data freqs;

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
	data->nsl_normed = NULL;
	data->H12 = NULL;
	data->H2H1 = NULL;	
	data->iSAFE = NULL;

	numExtraArgs = argc - nComparisons;
	sprintf(infilename, "%s", argv[numExtraArgs]);
	zinf = gzopen(infilename, "rb");
	gzgets(zinf, ihs_filename, line_size); 
	//strtok(ihs_filename, "\n");
	gzgets(zinf, delihh_filename, line_size); 
	//strtok(delihh_filename, "\n");
	gzgets(zinf, nsl_filename, line_size); 
	//strtok(nsl_filename, "\n");
	gzgets(zinf, H12_filename, line_size); 
	//strtok(H12_filename, "\n");
	gzgets(zinf, iSAFE_filename, line_size); 
	//strtok(iSAFE_filename, "\n");
	gzgets(zinf, xpehh_filename, line_size); 
	//strtok(xpehh_filename, "\n");
	gzgets(zinf, freqs_filename, line_size); 
	strtok(freqs_filename, "\n");
	gzclose(zinf);
	get_freqs_data(&freqs, freqs_filename, minPos, maxPos);
	nunique = freqs.nsnps;
	allUniqueSnps = calloc(nunique, sizeof(int));
	for (isnp = 0; isnp < nunique; isnp++){
		allUniqueSnps[isnp] = freqs.pos[isnp];
	} // end isnp
	free_freqs_data(&freqs);
	//UP TO THIS POINT, we just need to get allUniqueSnps[] loaded. with len nunique

	///////////////////////
	/// ALLOCATE MEMORY ///
	//////////////////////
	data->nsnps = nunique;
	data->ncomp = nComparisons;
	data->physpos = malloc(nComparisons * sizeof(int*));
	data->genpos = malloc(nComparisons * sizeof(double*));
	data->daf_selpop = malloc(nComparisons * sizeof(double*));
	data->delDAF = malloc(nComparisons * sizeof(double*));
	data->fst = malloc(nComparisons * sizeof(double*));
	data->xp_normed = malloc(nComparisons * sizeof(double*));
	data->ihs_normed = malloc(nComparisons * sizeof(double*));
	data->delihh_normed = malloc(nComparisons * sizeof(double*));
	data->nsl_normed = malloc(nComparisons * sizeof(double*));
	data->H12 = malloc(nComparisons * sizeof(double*));
	data->H2H1 = malloc(nComparisons * sizeof(double*));	
	data->iSAFE = malloc(nComparisons * sizeof(double*));

	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);
	assert(data->nsl_normed != NULL);
	assert(data->H12 != NULL);
	assert(data->H2H1 != NULL);	
	assert(data->iSAFE != NULL);

	for (iComp = 0; iComp < nComparisons; iComp++){
		data->physpos[iComp] = calloc(nunique, sizeof(int));
		data->genpos[iComp] = calloc(nunique, sizeof(double)); // problem
		data->daf_selpop[iComp] = calloc(nunique, sizeof(double));
		data->delDAF[iComp] = calloc(nunique, sizeof(double));
		data->fst[iComp] = calloc(nunique, sizeof(double));
		data->xp_normed[iComp] = calloc(nunique, sizeof(double));		
		data->ihs_normed[iComp] = calloc(nunique, sizeof(double));
		data->delihh_normed[iComp] = calloc(nunique, sizeof(double));	
		data->nsl_normed[iComp] = calloc(nunique, sizeof(double));	
		data->H12[iComp] = calloc(nunique, sizeof(double));		
		data->H2H1[iComp] = calloc(nunique, sizeof(double));				
		data->iSAFE[iComp] = calloc(nunique, sizeof(double));		
			
		assert(data->physpos[iComp] != NULL);
		assert(data->genpos[iComp] != NULL);
		assert(data->daf_selpop[iComp] != NULL);
		assert(data->delDAF[iComp] != NULL);
		assert(data->fst[iComp] != NULL);
		assert(data->xp_normed[iComp] != NULL);
		assert(data->ihs_normed[iComp] != NULL);
		assert(data->delihh_normed[iComp] != NULL);
		assert(data->nsl_normed[iComp] != NULL);		
		assert(data->H12[iComp] != NULL);		
		assert(data->H2H1[iComp] != NULL);				
		assert(data->iSAFE[iComp] != NULL);		
	} // end for icomp

	/////////////////////////////////////////////
	// LOAD ALL COMPARISONS TO ONE DATA OBJECT //
	/////////////////////////////////////////////
	fprintf(stderr, "Loading all component scores to one data-object...\n");
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);

		zinf = gzopen(infilename, "rb");
		gzgets(zinf, ihs_filename, line_size); 
		strtok(ihs_filename, "\n");
		gzgets(zinf, delihh_filename, line_size); 
		strtok(delihh_filename, "\n");
		gzgets(zinf, nsl_filename, line_size); 
		strtok(nsl_filename, "\n");
		gzgets(zinf, H12_filename, line_size); 
		strtok(H12_filename, "\n");
		gzgets(zinf, iSAFE_filename, line_size); 
		strtok(iSAFE_filename, "\n");
		gzgets(zinf, xpehh_filename, line_size); 
		strtok(xpehh_filename, "\n");
		gzgets(zinf, freqs_filename, line_size); 
		strtok(freqs_filename, "\n");
		gzclose(zinf);
		thisCompSnpCount = -1;//compCountSnps[iComp];
		//fprintf(stderr, "start call to get_pair_any");
		get_popPair_anyData(minPos, maxPos, &data_sing, ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename, thisCompSnpCount); // PASS COUNTS
		fprintf(stderr, "\tloaded pop-pair object with %d snps\n", data_sing.nsnps);

		jsnp = 0; //isnp iterates (0, nunique) over allUnique Snps; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
			if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
				data->physpos[iComp][isnp] = data_sing.physpos[jsnp];	
				data->genpos[iComp][isnp] = data_sing.genpos[jsnp];	 
				data->daf_selpop[iComp][isnp] = data_sing.daf_selpop[jsnp];	 
				data->delDAF[iComp][isnp] = data_sing.delDAF[jsnp];	
				data->fst[iComp][isnp] = data_sing.fst[jsnp];	 
				data->xp_normed[iComp][isnp] = data_sing.xp_normed[jsnp]; //fprintf(stderr, "%f\n", data_sing.xp_normed[jsnp]);		//AT THIS POINT SOMETHING HAS GONE WRONG					 
				data->ihs_normed[iComp][isnp] = data_sing.ihs_normed[jsnp];	 
				data->delihh_normed[iComp][isnp] = data_sing.delihh_normed[jsnp];	
				data->nsl_normed[iComp][isnp] = data_sing.nsl_normed[jsnp];		 
				data->H12[iComp][isnp] = data_sing.H12[jsnp];		
				data->H2H1[iComp][isnp] = data_sing.H2H1[jsnp];						 
				data->iSAFE[iComp][isnp] = data_sing.iSAFE[jsnp];		 

				jsnp++; //assert(jsnp<=data_sing.nsnps);
				if (jsnp > data_sing.nsnps){break;}
			}
			else if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){jsnp++; if (jsnp > data_sing.nsnps){break;}}//assert(jsnp<=data_sing.nsnps);}
			//else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
		}// end for isnp loop

		//last 
		data->physpos[iComp][nunique-1] = data_sing.physpos[data_sing.nsnps-1];
		data->genpos[iComp][nunique-1] = data_sing.genpos[data_sing.nsnps-1];
		data->daf_selpop[iComp][nunique-1] = data_sing.daf_selpop[data_sing.nsnps-1];
		data->delDAF[iComp][nunique-1] = data_sing.delDAF[data_sing.nsnps-1];
		data->fst[iComp][nunique-1] = data_sing.fst[data_sing.nsnps-1];
		data->xp_normed[iComp][nunique-1] = data_sing.xp_normed[data_sing.nsnps-1];
		data->ihs_normed[iComp][nunique-1] = data_sing.ihs_normed[data_sing.nsnps-1];
		data->nsl_normed[iComp][nunique-1] = data_sing.nsl_normed[data_sing.nsnps-1];
		data->delihh_normed[iComp][nunique-1] = data_sing.delihh_normed[data_sing.nsnps-1];
		data->H12[iComp][nunique-1] = data_sing.H12[data_sing.nsnps-1];
		data->H2H1[iComp][nunique-1] = data_sing.H2H1[data_sing.nsnps-1];		
		data->iSAFE[iComp][nunique-1] = data_sing.iSAFE[data_sing.nsnps-1];

		free_popPair_data(&data_sing); 
	} // end for icomp
	free(allUniqueSnps); 	//problem??
	//free(allSnps);
	free(newLine);
	//free(compCountSnps);
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
		free(data->nsl_normed[iComp]);	
		free(data->H12[iComp]);
		free(data->H2H1[iComp]);	
		free(data->iSAFE[iComp]);
	}
	free(data->physpos);
	free(data->genpos);
	free(data->daf_selpop);
	free(data->delDAF);
	free(data->fst);
	free(data->xp_normed);
	free(data->ihs_normed);
	free(data->delihh_normed);
	free(data->nsl_normed);
	free(data->H12);
	free(data->H2H1);		
	free(data->iSAFE);		
	data->nsnps = 0;
	data->ncomp = 0;
} //end method
float compareXp(popComp_data_multiple* data, int isnp){//currently: takes max val
	double xp;
	int iComp;
	xp = data->xp_normed[0][isnp];
	for (iComp = 0; iComp < data->ncomp; iComp++){
		//fprintf(stderr, "%f\n", data->xp_normed[iComp][isnp]);
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
	}
	ave = deldaf / (double)theseComp;
	return ave;
} //end function
float get_outgroups_fst(popComp_data_multiple* data, int isnp, int iComp, int jComp){
	double daf_alt1, daf_alt2;
	double daf_mean;
	double fst, num, denom;
	double msp, msg;

	int nperPop = 172; // quick estimate (?) slash accurate for sims.
	//fprintf(stderr, "must validate");
	/*
            pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
            nic = ni - (double) ni * ni / (ni + nj);
            njc = nj - (double) nj * nj / (ni + nj);
            nc = nic + njc;
            msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
            msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
            num = msp - msg;
            denom = msp + (nc - 1) * msg;
	*/
			//thiscomp_deldaf = data->delDAF[iComp][isnp];
			//thiscomp_daf_sel = data->daf_selpop[iComp][isnp];
		//	thiscomp_daf_altpop = thiscomp_daf_sel - thiscomp_deldaf; //VALIDATE
		//	alt_daf += thiscomp_daf_altpop;

	daf_alt1 = data->delDAF[iComp][isnp] - data->daf_selpop[iComp][isnp];
	daf_alt2 = data->delDAF[jComp][isnp] - data->daf_selpop[jComp][isnp];

	daf_mean = (daf_alt1 + daf_alt2) / 2.;
	msp = nperPop * (daf_alt1 - daf_mean) * (daf_alt1 - daf_mean) + nperPop * (daf_alt2 - daf_mean) * (daf_alt2 - daf_mean);
    msg = (nperPop * daf_alt1 * (1. - daf_alt1) + nperPop * daf_alt2 * (1. - daf_alt2)) / (nperPop - 1 + nperPop - 1);
    num = msp - msg;
    denom = msp + ((nperPop) - 1) * msg; //really unsure about this. possibly don't use weir-hill estimator
	fst = (num / denom);
	return fst;
} // end function
float get_PBS(double in_t_1, double in_t_2, double out_t){
	/* population branch statistic; Yi et al	*/
	float PBS;
	PBS = ((in_t_1 + in_t_2 - out_t) / 2.);
	return PBS;
} // end function
float get_T(double fst){//helper method for PBS; transforms Fst cf Cavalli-Sforza 1969
	return -1 * log(1. - fst); 
}// end function
float compareFst_PBS(popComp_data_multiple* data, int isnp){ // population branch statistic (pop vs. two outgroups)
	// Population-Branch Statistic; an population-specific generalization of Fst for three populations
	// Yi et al., Science 2013
	//double fst;
	//double ave;
	double maxPbs = 0; //if nPop > 3; get PBS for each pair of outgroups and take the maximum
	int iComp, jComp;//, theseComp=0;
	double in_fst_1, in_fst_2, out_fst;
	double in_t_1, in_t_2, out_t;
	double thisPbs;
	//Maybe we don't even use the values calculated by earlier C program (calc_freqs) -- we need Fst between outgroups anyway.
	//Maybe instead we just get all DAF values for which we have data, and get PBS for selpop for each grouping with two outgroups?
	//(Then take maximum?)
	//every unique pairing of two files we can take from the args given makes a triad. 
	// nested loop.
	for (iComp = 0; iComp < data->ncomp; iComp++){
		for (jComp = iComp + 1; jComp < data->ncomp; jComp++){
			if ((iComp != jComp) && (data->fst[iComp][isnp] != 0) && (data->fst[jComp][isnp] != 0)){
				// calc PBS
				// get three dafs (fuck do I need n_pop????)
				in_fst_1 = data->fst[iComp][isnp];
				in_fst_2 = data->fst[jComp][isnp];
				out_fst = get_outgroups_fst(data, isnp, iComp, jComp);
				in_t_1 = get_T(in_fst_1);
				in_t_2 = get_T(in_fst_2);
				out_t = get_T(out_fst);
				thisPbs = get_PBS(in_t_1, in_t_2, out_t);
				if (thisPbs > maxPbs){maxPbs = thisPbs;}
			} // end if
		} // end jComp loop
	} // end icomp Loop
	return maxPbs;
} // end function
float comparedelDaf_outgroup_ave(popComp_data_multiple* data, int isnp){//daf_thispop - AVE(outgroup dafs) [==cms1.0]
	double deldaf;
	int iComp, theseComp=0;
	double 	thiscomp_deldaf, thiscomp_daf_sel, thiscomp_daf_altpop; //for each comp, retrieve altpop daf from datastructure
	double ave_daf, alt_daf;
	alt_daf = 0;
	for (iComp = 0; iComp < data->ncomp; iComp++){
		if (data->delDAF[iComp][isnp] !=0){
			thiscomp_deldaf = data->delDAF[iComp][isnp];
			thiscomp_daf_sel = data->daf_selpop[iComp][isnp];
			thiscomp_daf_altpop = thiscomp_daf_sel - thiscomp_deldaf; //VALIDATE
			alt_daf += thiscomp_daf_altpop;
			theseComp++;
		} // end if
	} // end for iComp
	ave_daf = alt_daf / (double)theseComp;
	deldaf = thiscomp_daf_sel - ave_daf;
	return deldaf;
} //end function

