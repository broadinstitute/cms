// functions for handling cms component(+composite) score datastructures
// last updated: 11.23.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
//#include <zlib.h>
#include "cms_data.h"

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
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
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
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
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
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
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
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
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

void get_nsl_data(nsl_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->freq1 = NULL;
	data->sl0 = NULL;
	data->sl1 = NULL;
	data->nsl_unnormed = NULL;
	data->nsl_normed = NULL;
	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->freq1 = malloc(data->nsnps * sizeof(double*)); assert(data->freq1 != NULL);
	data->sl0 = malloc(data->nsnps * sizeof(double*)); assert(data->sl0 != NULL);
	data->sl1 = malloc(data->nsnps * sizeof(double*)); assert(data->sl1 != NULL);
	data->nsl_unnormed = malloc(data->nsnps * sizeof(double*)); assert(data->nsl_unnormed != NULL);
	data->nsl_normed = malloc(data->nsnps * sizeof(double*)); assert(data->nsl_normed != NULL);

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
				data->sl0[isnp] = atof(token);
			}			
			else if (itoken == 4) {
				data->sl1[isnp] = atof(token);
			}
			else if (itoken == 5) {
				data->nsl_unnormed[isnp] = atof(token);
			}			
			else if (itoken == 6) {
				data->nsl_normed[isnp] = atof(token);
			}

		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method

void free_nsl_data(nsl_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->freq1);
	free(data->sl0);
	free(data->sl0);
	free(data->nsl_unnormed);
	free(data->nsl_normed);
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
	//fprintf(stderr, "\tloading likes tables from %s\n", filename);

	inf = fopen(filename, "r");
	assert(inf != NULL);
	fgets(newLine, line_size, inf); //strip
	while (fgets(newLine, line_size, inf) != NULL){
		data->nbins++;
	}
	data->nbins++; //one extra because 50 bins can be represented with 49 bounds

	data->start_bin = malloc(data->nbins * sizeof(double));
	data->end_bin = malloc(data->nbins * sizeof(double));
	data->probs = malloc(data->nbins * sizeof(double));
	fclose(inf);

	inf = fopen(filename, "r");
	assert(inf != NULL);
	ibin = 0;
	while (fgets(newLine, line_size, inf) != NULL){
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t ")) != NULL; itoken++) {
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
		//fprintf(stderr, "ibin is now %d\n", ibin);
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

/****************/
/***POP PAIR****/
/***************/

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

	fprintf(stderr, "\tloading from %s\n", filename);

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

	fprintf(stderr, "\tloading from %s\n", filename);

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






