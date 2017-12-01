// functions for handling cms component(+composite) score datastructures
// last updated: 10.19.2017 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "cms_data.h"

/*********************/
/***HELPER METHODS****/ 
/*********************/
int intcmp(const void *v1, const void *v2) {return (*(int *)v1 - *(int *)v2);}
int count_unique(int arr[], int len) {
    int counted[len], j, n, count, flag;
    counted[0] = arr[0]; 
    count = 1;/*one element is counted*/
        for(j=0; j <= len-1; ++j) {
        flag = 1;;
        /*the counted array will always have 'count' elements*/
        for(n=0; n < count; ++n) {
            if(arr[j] == counted[n]) {
                flag = 0;
            }
        }
        if(flag == 1) {
            ++count;
            counted[count-1] = arr[j];
        }
    }
    return count;
} //end method

/**********************/
/***COMPONENT SCORES***/
/**********************/
void get_freqs_data(freqs_data* data, char filename[]) {
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
	data->popdaf = NULL;
	data->genpos = NULL;

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
	assert(inf != NULL);
	fgets(newLine, line_size, inf); // strip header
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->fst = malloc(data->nsnps * sizeof(double*)); assert(data->fst != NULL);
	data->deldaf = malloc(data->nsnps * sizeof(double*)); assert(data->deldaf != NULL);
	data->genpos = malloc(data->nsnps * sizeof(double*)); assert(data->genpos != NULL);
	data->popdaf = malloc(data->nsnps * sizeof(double*)); assert(data->popdaf != NULL);

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
				data->genpos[isnp] = atof(token);
			}
			else if (itoken == 2) {
				data->popdaf[isnp] = atof(token);
			}	
			else if (itoken == 3) {
				data->deldaf[isnp] = atof(token);
			}			
			else if (itoken == 4) {
				data->fst[isnp] = atof(token);
			}			
			
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
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
	fgets(newLine, line_size, inf); //strip header
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
	int unnormed_index, normed_index; // selscan output may or may not include ihh decomposition; account for both

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

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
	assert(inf != NULL);
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// determine formatting based on line (selscan output may or may not include ihh decomposition)
	// https://github.com/szpiech/selscan
	unnormed_index = 0; //5 9
	normed_index = 0; //6 10
	int dummy_var = 0;
	for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++)
		{dummy_var ++;}

	if (itoken <= 8){unnormed_index = 5; normed_index = 6;}
	else {unnormed_index = 9; normed_index = 10;}


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
			else if (itoken == unnormed_index){
				data->ihs_unnormed[isnp] = atof(token);
			}			
			else if (itoken == normed_index){
				data->ihs_normed[isnp] = atof(token);
			}
			//else if (itoken == normed_index + 1) {
			//	data->lastcol[isnp] = atoi(token);
			//}
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
	free(data->sl1);
	free(data->nsl_unnormed);
	free(data->nsl_normed);
	data->nsnps = 0;
} //end method
void get_H12_data(H12_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->H12_value = NULL;
	data->H2H1_value = NULL;
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
	data->H12_value = malloc(data->nsnps * sizeof(double*)); assert(data->H12_value != NULL);
	data->H2H1_value = malloc(data->nsnps * sizeof(double*)); assert(data->H2H1_value != NULL);

	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 0) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 8) {
				data->H12_value[isnp] = atof(token);
			}
			else if (itoken == 9) {
				data->H2H1_value[isnp] = atof(token);
			}			
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method
void free_H12_data(H12_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->H12_value);
	free(data->H2H1_value);
	data->nsnps = 0;
} //end method
void get_iSAFE_data(iSAFE_data* data, char filename[]) {
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	int isnp, itoken;

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->pos = NULL; 
	data->iSAFE_value = NULL;

	inf = fopen(filename, "r");
	if (inf == NULL) {fprintf(stderr, "Missing file: %s\n", filename);}
	assert(inf != NULL);
	fgets(newLine, line_size, inf); //strip header
	while (fgets(newLine, line_size, inf) != NULL) {
			assert(strlen(newLine) < line_size);
			data->nsnps++;
		}
	fclose(inf);

	// Allocate memory; initialize
	data->pos = malloc(data->nsnps * sizeof(int*)); assert(data->pos != NULL);
	data->iSAFE_value = malloc(data->nsnps * sizeof(double*)); assert(data->iSAFE_value != NULL);

	/*******************
	GET DATA FROM FILE
	*******************/
	inf = fopen(filename, "r");
	fgets(newLine, line_size, inf); //strip header
	isnp = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, " \t")) != NULL; itoken++) {
			if (itoken == 0) {
				data->pos[isnp] = atoi(token);
			}	
			else if (itoken == 1) {
				data->iSAFE_value[isnp] = atof(token);
			}
		} // END for running=newLine
		isnp++;
	} //END while(fgets(newLine))
	
	fclose(inf);
	free(newLine);
} //end method
void free_iSAFE_data(iSAFE_data* data) {
	if (data == NULL) {return;}
	free(data->pos);
	free(data->iSAFE_value);
	data->nsnps = 0;
} //end method

/************************/
/***SCORE LIKELIHOODS***/
/************************/
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

	data->start_bin = calloc(data->nbins, sizeof(double));
	data->end_bin = calloc(data->nbins, sizeof(double));
	data->probs = calloc(data->nbins, sizeof(double));
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
void get_likes_data_multiple(likes_data_multiple* data, char filename[]){
	FILE *inf=NULL;
	const int line_size = 15000000; 
    char miss_probs_filename[256], hit_probs_hi_filename[256], hit_probs_mid_filename[256], hit_probs_low_filename[256];   
    likes_data onedist_data;
    int ibin;

	//////////////////
	/// INITIALIZE ///
	//////////////////
	data->nbins = 0;
	data->start_bin = NULL;
	data->end_bin = NULL;
	data->miss_probs = NULL;
	data->hit_probs = NULL;
	//data->hit_probs_mid = NULL;
	//data->hit_probs_low = NULL;
	data->hit_probs = malloc(3 * sizeof(double*)); //likesFreqs: hi, low, mid

	inf = fopen(filename, "r"); 
	assert(inf != NULL);
	// P(NEUT)
	fgets(miss_probs_filename, line_size, inf);
	strtok(miss_probs_filename, "\n");
	get_likes_data(&onedist_data, miss_probs_filename);
	data->nbins = onedist_data.nbins;
	data->start_bin = calloc(data->nbins, sizeof(double));
	data->end_bin = calloc(data->nbins, sizeof(double));
	data->miss_probs = calloc(data->nbins, sizeof(double));
	//data->hit_probs_hi = malloc(data->nbins * sizeof(double));	
	//data->hit_probs_mid = malloc(data->nbins * sizeof(double));	
	//data->hit_probs_low = malloc(data->nbins * sizeof(double));
	data->hit_probs[0] = calloc(data->nbins, sizeof(double));
	data->hit_probs[1] = calloc(data->nbins, sizeof(double));
	data->hit_probs[2] = calloc(data->nbins, sizeof(double));
				
	for (ibin = 0; ibin < data->nbins; ibin++){
		data->start_bin[ibin] = onedist_data.start_bin[ibin];
		data->end_bin[ibin] = onedist_data.end_bin[ibin];
		data->miss_probs[ibin] = onedist_data.probs[ibin];
	} // end ibin
	free_likes_data(&onedist_data);

	//P(SEL | DAF HI)
	fgets(hit_probs_hi_filename, line_size, inf);
	strtok(hit_probs_hi_filename, "\n");
	get_likes_data(&onedist_data, hit_probs_hi_filename);
	for (ibin = 0; ibin < data->nbins; ibin++){
		data->hit_probs[0][ibin] = onedist_data.probs[ibin];
	} // end ibin
	free_likes_data(&onedist_data);

	//P(SEL | DAF MID)
	fgets(hit_probs_mid_filename, line_size, inf);
	strtok(hit_probs_mid_filename, "\n");
	get_likes_data(&onedist_data, hit_probs_mid_filename);
	for (ibin = 0; ibin < data->nbins; ibin++){
		data->hit_probs[1][ibin] = onedist_data.probs[ibin];
	} // end ibin
	free_likes_data(&onedist_data);

	//P(SEL | DAF LOW)
	fgets(hit_probs_low_filename, line_size, inf);
	strtok(hit_probs_low_filename, "\n");
	get_likes_data(&onedist_data, hit_probs_low_filename);
	for (ibin = 0; ibin < data->nbins; ibin++){
		data->hit_probs[2][ibin] = onedist_data.probs[ibin];
	} // end ibin
	free_likes_data(&onedist_data);

	fclose(inf);
} // end function
void free_likes_data_multiple(likes_data_multiple* data){
    if (data == NULL) {return;}
	data->nbins = 0;
	free(data->start_bin);
	free(data->end_bin);
	free(data->miss_probs);
	//free(data->hit_probs_hi); //hit_probs
	//free(data->hit_probs_mid);
	//free(data->hit_probs_low);
	free(data->hit_probs[0]);
	free(data->hit_probs[1]);
	free(data->hit_probs[2]);
	free(data->hit_probs);
} // end function
float getHitProb(likes_data_multiple* data, int likesIndex, double value){
	int ibin;
	float hitProb = NAN;
	//fprintf(stderr, "gethitprob: %f\n", value);
	for (ibin = 0; ibin < data->nbins; ibin++){
		//fprintf(stderr, "ibin %d \t", ibin);
		if (value >= data->start_bin[ibin] && value <= data->end_bin[ibin]){return data->hit_probs[likesIndex][ibin];}
	}
	if (value < data->start_bin[0]){return data->hit_probs[likesIndex][0];}
	if (value > data->end_bin[data->nbins - 1]){return data->hit_probs[likesIndex][data->nbins - 1];}
	return hitProb;
} //end function
float getMissProb(likes_data_multiple* data, double value){
	int ibin;
	float missProb = NAN;
	//fprintf(stderr, "getmissprob: %f\n", value);
	for (ibin = 0; ibin < data->nbins; ibin++){
		if (value >= data->start_bin[ibin] && value <= data->end_bin[ibin]){return data->miss_probs[ibin];}
	}
	if (value < data->start_bin[0]){return data->miss_probs[0];}
	if (value > data->end_bin[data->nbins - 1]){return data->miss_probs[data->nbins - 1];}
	return missProb;
} //end function
float getMaxBf(likes_data_multiple* data, int likesIndex){
	int ibin;
	float thisBf;
	float maxBf = 0.;
	for (ibin = 0; ibin < data->nbins; ibin++){
		if(data->hit_probs[likesIndex][ibin] > 1e-10 && data->miss_probs[ibin] > 1e-10){
			thisBf = data->hit_probs[likesIndex][ibin] / data->miss_probs[ibin];
			if (thisBf > maxBf){maxBf = thisBf;}
		}
	}//end ibin
	//fprintf(stderr, "found max bf: %f\n", maxBf);
	return maxBf;
}//end function
float getMinBf(likes_data_multiple* data, int likesIndex){
	int ibin;
	float thisBf;
	float minBf = 1.;
	for (ibin = 0; ibin < data->nbins; ibin++){
		if(data->hit_probs[likesIndex][ibin] > 1e-10 && data->miss_probs[ibin] > 1e-10){
			thisBf = data->hit_probs[likesIndex][ibin] / data->miss_probs[ibin];
			if (thisBf < minBf && thisBf != 0){minBf = thisBf;}
		}
	}//end ibin
	//fprintf(stderr, "found min bf: %f\n", minBf);	
	return minBf;
}//end function
float getMaxProb(likes_data_multiple* data, int likesIndex, double prior){
	int ibin;
	float thisProb;
	float maxProb = 0.;
	for (ibin = 0; ibin < data->nbins; ibin++){
		if(data->hit_probs[likesIndex][ibin] > 1e-10 && data->miss_probs[ibin] > 1e-10){
			thisProb = (prior * data->hit_probs[likesIndex][ibin]) / ((prior * data->hit_probs[likesIndex][ibin]) + ((1. - prior) * data->miss_probs[ibin]));
			if (thisProb > maxProb){maxProb = thisProb;}
		}
	}//end ibin
	//fprintf(stderr, "found max prob: %f\n", maxProb);
	return maxProb;
}//end function
float getMinProb(likes_data_multiple* data, int likesIndex, double prior){
	int ibin;
	float thisProb;
	float minProb = 1.;
	for (ibin = 0; ibin < data->nbins; ibin++){
		if(data->hit_probs[likesIndex][ibin] > 1e-10 && data->miss_probs[ibin] > 1e-10){
			thisProb = (prior * data->hit_probs[likesIndex][ibin]) / ((prior * data->hit_probs[likesIndex][ibin]) + ((1. - prior) * data->miss_probs[ibin]));
			if (thisProb < minProb){minProb = thisProb;}
		}
	}//end ibin
	//fprintf(stderr, "found min prob: %f\n", minProb);	
	return minProb;
}//end function

/****************/
/***POP PAIR****/ //ADD ISAFE HERE
/***************/
int get_num_completeData(char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]){
	/* 	counts the number of (unique) SNPs that have ALL component scores present */
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;
	int nsnps;
	int ihs1_index, delihh1_index, nsl1_index, xp_index, freqs_index, H12_index, iSAFE_index;
	int ihs1pos, delihh1pos, nsl1pos, xppos, freqspos, h12pos, isafepos, minimum;

	get_ihs_data(&ihs1, ihs_filename);
	get_delihh_data(&delihh1, delihh_filename);
	get_nsl_data(&nsl1, nsl_filename);
	get_xpehh_data(&xp, xpehh_filename);
	get_freqs_data(&freqs, freqs_filename);
	get_H12_data(&H12, H12_filename);
	get_iSAFE_data(&iSAFE, iSAFE_filename);

	ihs1_index=0;
	nsl1_index=0;
	delihh1_index=0;
	xp_index=0;
	freqs_index=0;
	H12_index=0;
	iSAFE_index=0;
	nsnps = 0;
	while (ihs1_index < ihs1.nsnps && delihh1_index < delihh1.nsnps && xp_index < xp.nsnps && freqs_index < freqs.nsnps && H12_index < H12.nsnps && iSAFE_index < iSAFE.nsnps)
	{
		ihs1pos = ihs1.pos[ihs1_index];
		nsl1pos = nsl1.pos[nsl1_index];
		h12pos = H12.pos[H12_index];
		isafepos = iSAFE.pos[iSAFE_index];
		delihh1pos = delihh1.pos[delihh1_index];
		xppos = xp.pos[xp_index];
		freqspos = freqs.pos[freqs_index];
		if (ihs1pos == delihh1pos && delihh1pos == xppos && xppos == freqspos && nsl1pos == freqspos && h12pos == freqspos && isafepos == h12pos){
			nsnps++;
		}
		//If not at the same point, find out which position is lowest and advance its pointer.
		minimum = 2147483647;
		if (ihs1pos < minimum){minimum = ihs1pos;}
		if (nsl1pos < minimum){minimum = nsl1pos;}
		if (delihh1pos < minimum){minimum = delihh1pos;}
		if (xppos < minimum){minimum = xppos;}		
		if (freqspos < minimum){minimum = freqspos;}
		if (h12pos < minimum){minimum = h12pos;}
		if (isafepos < minimum){minimum = isafepos;}

		if (ihs1pos == minimum){ihs1_index++;}
		if (nsl1pos == minimum){nsl1_index++;}
		if (delihh1pos == minimum){delihh1_index++;}
		if (xppos == minimum){xp_index++;}		
		if (freqspos == minimum){freqs_index++;}
		if (h12pos == minimum){H12_index++;}	
		if (isafepos == minimum){iSAFE_index++;}
	} // end while loop
	//fprintf(stderr, "\tfound complete data for %d SNPs\n", data->nsnps);
	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
	return nsnps;
} // end function
int get_num_anyData(char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]){
	/* counts the number of unique SNPs that have a score value for at least one component score */
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;
	int nunique;
	int totNsnp, isnp;//, jsnp;
	int *allSnps, *allUniqueSnps;
	int cumulativeSNPCount = 0;

	//////////////////////
	// LOAD IN ALL DATA //
	//////////////////////
	get_ihs_data(&ihs1, ihs_filename);
	get_delihh_data(&delihh1, delihh_filename);
	get_nsl_data(&nsl1, nsl_filename);
	get_xpehh_data(&xp, xpehh_filename);
	get_freqs_data(&freqs, freqs_filename);
	get_H12_data(&H12, H12_filename);
	get_iSAFE_data(&iSAFE, iSAFE_filename);
	/*
	fprintf(stderr, "\tloading component score data from: %s\n", ihs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n", ihs1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", delihh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  delihh1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", nsl_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  nsl1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", xpehh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  xp.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", freqs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  freqs.nsnps);
	*/
	totNsnp = ihs1.nsnps + delihh1.nsnps + nsl1.nsnps + xp.nsnps + freqs.nsnps + H12.nsnps + iSAFE.nsnps;
	allSnps = malloc(totNsnp * sizeof(int));

	for (isnp = 0; isnp < ihs1.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = ihs1.pos[isnp];
	}	// end ihs isnp
	cumulativeSNPCount += ihs1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+delihh1.nsnps; isnp++){
		//fprintf(stderr, "\tisnp: %d\n", isnp);
		allSnps[isnp] = delihh1.pos[isnp-cumulativeSNPCount];
	}	//end delihh isno
	cumulativeSNPCount += delihh1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+nsl1.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = nsl1.pos[isnp-cumulativeSNPCount];
	}	//end nsl isnp
	cumulativeSNPCount += nsl1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+xp.nsnps; isnp++){
		//fprintf(stderr, "\tisnp: %d\n", isnp);
		allSnps[isnp] = xp.pos[isnp-cumulativeSNPCount];
	}	//end xp isnp
	cumulativeSNPCount += xp.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+freqs.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = freqs.pos[isnp-cumulativeSNPCount];
	}	//end freq isnp
	cumulativeSNPCount += freqs.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount + H12.nsnps; isnp++){
		allSnps[isnp] = H12.pos[isnp-cumulativeSNPCount];
	} // end H12 isnp
	cumulativeSNPCount += H12.nsnps;
	for (isnp = cumulativeSNPCount; isnp < totNsnp; isnp++){
		allSnps[isnp] = H12.pos[isnp-cumulativeSNPCount];
	} // end H12 isnp

	//fprintf(stderr, "get_num_anyData returns totnsnp: %d\n", totNsnp);

	/*
	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	nunique = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{nunique++;}
	} // end for isnp
	if (allSnps[totNsnp-1] != allSnps[totNsnp-2]){nunique++;} //edge case
	*/
	nunique = count_unique(allSnps, totNsnp);
	allUniqueSnps= malloc(nunique * sizeof(int)); //is this step necessary?
	/*
	jsnp = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
	} // end for isnp
	if (allSnps[totNsnp-1] != allSnps[totNsnp-2]){allUniqueSnps[jsnp] = allSnps[totNsnp-1];} //edge case
	//BUT IF THIS IS NOT THE CASE WHAT HAPPENS TO LAST SNP
	*/
	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
	free(allSnps);
	free(allUniqueSnps);
	return nunique;
} // end function
void get_popPair_completeData(popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]){
	/* loads in COMPLETE data for a population pair */
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;

	char *locus;
	int isnp;
	int ihs1_index, delihh1_index, nsl1_index, xp_index, freqs_index, H12_index, iSAFE_index;
	int ihs1pos, delihh1pos, nsl1pos, xppos, freqspos, h12pos, isafepos, minimum;
	double thisDeldaf, thisXp;

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

	locus = malloc(256 * sizeof(char));


	////////////////////////////////
	// LOAD EACH COMPONENT SCORE ///
	////////////////////////////////
	get_ihs_data(&ihs1, ihs_filename);
	get_delihh_data(&delihh1, delihh_filename);
	get_nsl_data(&nsl1, nsl_filename);
	get_xpehh_data(&xp, xpehh_filename);
	get_freqs_data(&freqs, freqs_filename);
	get_H12_data(&H12, H12_filename);
	get_iSAFE_data(&iSAFE, iSAFE_filename);

	/*
	fprintf(stderr, "\tloading component score data from: %s\n", ihs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n", ihs1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", delihh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  delihh1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", nsl_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  nsl1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", xpehh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  xp.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", freqs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  freqs.nsnps);
	*/

	///////////////////
	// COLLATE SNPS ///
	///////////////////
	ihs1_index=0;
	nsl1_index=0;
	delihh1_index=0;
	xp_index=0;
	freqs_index=0;
	H12_index=0;	
	iSAFE_index=0;

	while (ihs1_index < ihs1.nsnps && delihh1_index < delihh1.nsnps && xp_index < xp.nsnps && freqs_index < freqs.nsnps && H12_index < freqs.nsnps && iSAFE_index < iSAFE.nsnps)
	{
		ihs1pos = ihs1.pos[ihs1_index];
		nsl1pos = nsl1.pos[nsl1_index];
		delihh1pos = delihh1.pos[delihh1_index];
		xppos = xp.pos[xp_index];
		freqspos = freqs.pos[freqs_index];
		h12pos = H12.pos[H12_index];
		isafepos = iSAFE.pos[iSAFE_index];
		if (ihs1pos == delihh1pos && delihh1pos == xppos && xppos == freqspos && nsl1pos == freqspos && freqspos == h12pos && h12pos == isafepos){
			data->nsnps++;
		}
		//If not at the same point, find out which position is lowest and advance its pointer.
		minimum = 2147483647;
		if (ihs1pos < minimum){minimum = ihs1pos;}
		if (nsl1pos < minimum){minimum = nsl1pos;}
		if (delihh1pos < minimum){minimum = delihh1pos;}
		if (xppos < minimum){minimum = xppos;}		
		if (freqspos < minimum){minimum = freqspos;}
		if (h12pos < minimum){minimum = h12pos;}
		if (isafepos < minimum){minimum = isafepos;}

		if (ihs1pos == minimum){ihs1_index++;}
		if (nsl1pos == minimum){nsl1_index++;}
		if (delihh1pos == minimum){delihh1_index++;}
		if (xppos == minimum){xp_index++;}		
		if (freqspos == minimum){freqs_index++;}
		if (h12pos == minimum){H12_index++;}
		if (isafepos == minimum){iSAFE_index++;}		
	} // end while loop
	//fprintf(stderr, "\tfound complete data for %d SNPs\n", data->nsnps);

	/////////////////////
	// RESERVE MEMORY ///
	/////////////////////
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
	data->nsl_normed = malloc(data->nsnps * sizeof(double));
	data->H12 = malloc(data->nsnps * sizeof(double));
	data->H2H1 = malloc(data->nsnps * sizeof(double));
	data->iSAFE = malloc(data->nsnps * sizeof(double));
	
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
	//fprintf(stderr, "reserved\n"); //debug

	/////////////////////////
	// LOAD COLLATED DATA ///
	/////////////////////////
	ihs1_index=0;
	nsl1_index=0;
	delihh1_index=0;
	xp_index=0;
	freqs_index=0;
	H12_index=0;
	iSAFE_index=0;

	isnp = -1;
	while (ihs1_index < ihs1.nsnps && delihh1_index < delihh1.nsnps && xp_index < xp.nsnps && freqs_index < freqs.nsnps && H12_index < H12.nsnps && iSAFE_index < iSAFE.nsnps)
	{
		ihs1pos = ihs1.pos[ihs1_index];
		nsl1pos = nsl1.pos[nsl1_index];
		delihh1pos = delihh1.pos[delihh1_index];
		xppos = xp.pos[xp_index];
		freqspos = freqs.pos[freqs_index];
		h12pos = H12.pos[H12_index];
		isafepos = iSAFE.pos[iSAFE_index];

		if (ihs1pos == delihh1pos && delihh1pos == xppos && xppos == freqspos && nsl1pos == freqspos && freqspos == h12pos && h12pos == isafepos){
			thisXp=xp.xpehh_normed[xp_index];
			thisDeldaf=freqs.deldaf[freqs_index];
			//fprintf(stderr, "found match\n"); //debug
			isnp +=1;
			sprintf(locus, "%d", ihs1pos); //leaving this for now; issues with copying string/pointer-to-char
			strcpy(data->locus_id[isnp], locus);	
			//fprintf(stderr, "isnp is%d\n", isnp); //debug
			data->physpos[isnp] = ihs1pos;
			data->genpos[isnp] = xp.genpos[xp_index];
			data->daf_selpop[isnp] = (1. - ihs1.freq1[ihs1_index]);
			data->delDAF[isnp] = thisDeldaf;
			data->fst[isnp] = freqs.fst[freqs_index];
			data->xp_normed[isnp] = thisXp;			
			data->ihs_normed[isnp] = ihs1.ihs_normed[ihs1_index];
			data->nsl_normed[isnp] = nsl1.nsl_normed[nsl1_index];
			data->delihh_normed[isnp] = delihh1.delihh_normed[delihh1_index];
			data->H12[isnp] = H12.H12_value[H12_index];
			data->H2H1[isnp] = H12.H2H1_value[H12_index];
			data->iSAFE[isnp] = iSAFE.iSAFE_value[iSAFE_index];
			
			//fprintf(stderr, "finished data entry\n");
		}
		//If not at the same point, find out which position is lowest and advance its pointer.
		minimum = 2147483647;
		if (ihs1pos < minimum){minimum = ihs1pos;}
		if (nsl1pos < minimum){minimum = nsl1pos;}
		if (delihh1pos < minimum){minimum = delihh1pos;}
		if (xppos < minimum){minimum = xppos;}		
		if (freqspos < minimum){minimum = freqspos;}
		if (h12pos < minimum){minimum = h12pos;}
		if (isafepos < minimum){minimum = isafepos;}

		//fprintf(stderr, "advancing\n"); //debug
		if (ihs1pos == minimum){ihs1_index++;}
		if (nsl1pos == minimum){nsl1_index++;}
		if (delihh1pos == minimum){delihh1_index++;}
		if (xppos == minimum){xp_index++;}		
		if (freqspos == minimum){freqs_index++;}
		if (h12pos == minimum){H12_index++;}	
		if (isafepos == minimum){iSAFE_index++;}	

	} //end while
	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
	free(locus);
	//fprintf(stderr, "loaded all data to object\n");
} //end method
void get_popPair_anyData(popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]){
	/*	loads in ANY data for a population pair */
	ihs_data ihs1;
	delihh_data delihh1;
	nsl_data nsl1;
	xpehh_data xp;
	freqs_data freqs;
	H12_data H12;
	iSAFE_data iSAFE;
	//char *locus;
	int isnp, jsnp;
	int ksnp, lsnp, msnp, osnp, psnp, qsnp;
	int totNsnp; //ie, redundant
	int *allSnps;
	int cumulativeSNPCount = 0;

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

	//locus = malloc(256 * sizeof(char));

	/////////////////////
	// RESERVE MEMORY ///
	/////////////////////
	data->nsnps = get_num_anyData(ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename); //fprintf(stderr, "loaded %d snps for pop-pair_any\n", data->nsnps);
	data->locus_id = malloc(data->nsnps * sizeof(char*));
	for (isnp = 0; isnp < data->nsnps; isnp++){
		data->locus_id[isnp] = malloc(256*sizeof(char));
		assert(data->locus_id[isnp] != NULL);
	}
	data->physpos = calloc(data->nsnps, sizeof(int)); //******** issue of matching zeros when this gets handled by ->multiple object...
	data->genpos = calloc(data->nsnps, sizeof(double)); //call ? quick fix
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

	////////////////////////////////////////
	// LOAD COMPONENT SCORES BY POSITION ///
	////////////////////////////////////////
	get_ihs_data(&ihs1, ihs_filename);
	get_delihh_data(&delihh1, delihh_filename);
	get_nsl_data(&nsl1, nsl_filename);
	get_xpehh_data(&xp, xpehh_filename);
	get_freqs_data(&freqs, freqs_filename);
	get_H12_data(&H12, H12_filename);
	get_iSAFE_data(&iSAFE, iSAFE_filename);
	/*fprintf(stderr, "\tloading component score data from: %s\n", ihs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n", ihs1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", delihh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  delihh1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", nsl_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  nsl1.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", xpehh_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  xp.nsnps);
	fprintf(stderr, "\tloading component score data from: %s\n", freqs_filename);
	fprintf(stderr, "\t\t found values for %d SNPs\n",  freqs.nsnps);*/

	totNsnp = ihs1.nsnps + delihh1.nsnps + nsl1.nsnps + xp.nsnps + freqs.nsnps + H12.nsnps + iSAFE.nsnps;
	allSnps = malloc(totNsnp * sizeof(int));
	for (isnp = 0; isnp < ihs1.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = ihs1.pos[isnp];
	}	// end ihs isnp
	cumulativeSNPCount += ihs1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+delihh1.nsnps; isnp++){
		//fprintf(stderr, "\tisnp: %d\n", isnp);
		allSnps[isnp] = delihh1.pos[isnp-cumulativeSNPCount];
	}	//end delihh isno
	cumulativeSNPCount += delihh1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+nsl1.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = nsl1.pos[isnp-cumulativeSNPCount];
	}	//end nsl isnp
	cumulativeSNPCount += nsl1.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+xp.nsnps; isnp++){
		//fprintf(stderr, "\tisnp: %d\n", isnp);
		allSnps[isnp] = xp.pos[isnp-cumulativeSNPCount];
	}	//end xp isnp
	cumulativeSNPCount += xp.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount+freqs.nsnps; isnp++){
		//fprintf(stderr, "isnp: %d\n", isnp);
		allSnps[isnp] = freqs.pos[isnp-cumulativeSNPCount];
	}	//end freq isnp
	cumulativeSNPCount += freqs.nsnps;
	for (isnp = cumulativeSNPCount; isnp < cumulativeSNPCount + H12.nsnps; isnp++){
		allSnps[isnp] = H12.pos[isnp-cumulativeSNPCount];
	} // end H12 isnp
	cumulativeSNPCount += H12.nsnps;
	for (isnp = cumulativeSNPCount; isnp < totNsnp; isnp++){
		allSnps[isnp] = H12.pos[isnp-cumulativeSNPCount];
	} // end H12 isnp


	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	//fprintf(stderr,"\nsorted\n");

	//////////////////////////////
	/// load all redundant snps //
	//////////////////////////////
	jsnp = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{data->physpos[jsnp] = allSnps[isnp]; jsnp++;}
	} // end isnp
	data->physpos[data->nsnps-1]  = allSnps[totNsnp-1]; // edge snp
	
	//fprintf(stderr, "loaded unique positions: %d\t%d\t%d\t%d\n", data->physpos[0], data->physpos[1], data->physpos[data->nsnps-2], data->physpos[data->nsnps-1]);
	//fprintf(stderr, "loaded unique scores: %f\t%f\t%f\t%f\t%f\t%f\n", ihs1.ihs_normed[0], delihh1.delihh_normed[0], nsl1.nsl_normed[0], xp.xpehh_normed[0], freqs.fst[0], freqs.deldaf[0]);

	jsnp = 0;	//ihs
	ksnp = 0;	//delihh
	lsnp = 0;	//nsl
	msnp = 0;	//xpehh	
	osnp = 0;	//freqs	
	psnp = 0; //H12
	qsnp = 0;	//iSAFE
	for (isnp = 0; isnp < data->nsnps; isnp++){

		//ihs (jsnp)
		//if this database_snp is the one I'm looking at, record it and advance pointer
		if(data->physpos[isnp] == ihs1.pos[jsnp]){data->ihs_normed[isnp] = ihs1.ihs_normed[jsnp]; if(jsnp < ihs1.nsnps - 1){jsnp++;}} 
		//if this database_snp < the snp I'm looking at, then record null.
		//if(data->physpos[isnp] < ihs1.pos[jsnp])
		else{data->ihs_normed[isnp] = 0.0 / 0.0;} //NaN

		//delihh (ksnp)
		if(data->physpos[isnp] == delihh1.pos[ksnp]){data->delihh_normed[isnp] = delihh1.delihh_normed[ksnp]; if(ksnp < delihh1.nsnps - 1){ksnp++;}}
		//if(data->physpos[isnp] < delihh1.pos[ksnp])
		else{data->delihh_normed[isnp] = 0.0 / 0.0;} //NaN
	
		//nsl (lsnp)
		if(data->physpos[isnp] == nsl1.pos[lsnp]){data->nsl_normed[isnp] = nsl1.nsl_normed[lsnp];if(lsnp < nsl1.nsnps -1){lsnp++;}}
		//if(data->physpos[isnp] < nsl1.pos[lsnp])
		else{data->nsl_normed[isnp] = 0.0 / 0.0;} //NaN

		//xp-ehh (msl)
		if(data->physpos[isnp] == xp.pos[msnp]){data->xp_normed[isnp] = xp.xpehh_normed[msnp];if(msnp < xp.nsnps - 1){msnp++;}}
		//if(data->physpos[isnp] < xp.pos[msnp])
		else{data->xp_normed[isnp] = 0.0 / 0.0;} //NaN

		//freqs (osnp)
		if(data->physpos[isnp] == freqs.pos[osnp]){data->fst[isnp] = freqs.fst[osnp]; data->delDAF[isnp] = freqs.deldaf[osnp]; data->genpos[isnp] = freqs.genpos[osnp]; data->daf_selpop[isnp] = freqs.popdaf[osnp]; if(osnp < freqs.nsnps - 1){osnp++;}}
		//if(data->physpos[isnp] < freqs.pos[osnp])
		else{data->fst[isnp] = 0.0 / 0.0; data->delDAF[isnp] = 0.0 / 0.0; } //NaN

		//h12 (psnp)
		if(data->physpos[isnp] == H12.pos[psnp]){data->H12[isnp] = H12.H12_value[psnp]; data->H2H1[isnp] = H12.H2H1_value[psnp]; if(psnp < H12.nsnps - 1){psnp++;}}
		//if(data->physpos[isnp] < freqs.pos[osnp])
		else{data->H12[isnp] = 0.0 / 0.0; data->H2H1[isnp] = 0.0 / 0.0;} //NaN

		//iSAFE (qsnp)
		if(data->physpos[isnp] == iSAFE.pos[qsnp]){data->iSAFE[isnp] = iSAFE.iSAFE_value[qsnp]; if(qsnp < iSAFE.nsnps - 1){qsnp++;}}
		//if(data->physpos[isnp] < freqs.pos[osnp])
		else{data->iSAFE[isnp] = 0.0 / 0.0;} //NaN

		//fprintf(stderr, "isnp: %d\tjsnp: %d\tksnp: %d\tlsnp: %d\tmsnp:%d\tosnp:%d\n", isnp, jsnp, ksnp, lsnp, msnp, osnp); //DEBUG
	} // end isnp
	
	//DEBUG
	/*for (isnp = 0; isnp < data->nsnps; isnp++){
		fprintf(stderr, "%d\t%f\n", data->physpos[isnp], data->ihs_normed[isnp]);
	}*/

	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_freqs_data(&freqs);
	free_H12_data(&H12);
	free_iSAFE_data(&iSAFE);
	//free(locus);
	free(allSnps);
	

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
void get_popComp_completeData(popComp_data_multiple* data, int nComparisons, int argc, char *argv[]){
	/*
	06.12.17: 	preserve this function - only loads into data those SNPs for which there exist ALL component scores.
				updating cms_local to cms_local_amend theoretically should obviate this.
	argv is all files for pop-pairs (each of which points to further component score files)
	*/
	const int line_size = 15000000; 
	popPair_data data_sing;
	FILE *inf=NULL;
	char *newLine;//, *token, *running;
	char infilename[512];
	int	isnp, jsnp, iComp, totNsnp, nunique; //thisPhysPos, itoken
	int *allSnps, *allUniqueSnps;
	char ihs_filename[528], delihh_filename[528], nsl_filename[528], xpehh_filename[528], freqs_filename[528];
	char H12_filename[528], iSAFE_filename[528]; 
	int numExtraArgs;

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

	////////////////////
	/// COLLATE LOCI ///
	////////////////////
	totNsnp = 0; //pass over first argument (run params)
	numExtraArgs = argc - nComparisons;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);
		fprintf(stderr, "loading pop-pair from file: ");
		fprintf(stderr, infilename);
		fprintf(stderr, "\n");
		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");
		fgets(H12_filename, line_size, inf);
		strtok(H12_filename, "\n");
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);
		totNsnp += get_num_completeData(ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename);
	} // end iComp
	fprintf(stderr, "found a total of %d complete-info snps present in any of %d comparisons.\n", totNsnp,  nComparisons);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		//get each infilename
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);
		//fprintf(stderr, "loading pop-pair from file: ");
		//fprintf(stderr, infilename);
		//fprintf(stderr, "\n");
		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");		
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");
		fgets(H12_filename, line_size, inf);
		strtok(H12_filename, "\n");
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);
		get_popPair_completeData(&data_sing, ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename); 
			
		for (jsnp = 0; jsnp < data_sing.nsnps; jsnp++){
			allSnps[isnp] = data_sing.physpos[jsnp];
			isnp +=1;
		}
		free_popPair_data(&data_sing); 
	}//end icomp
	//fprintf(stderr, "loaded all positions!\n");
	//for (ivar = 0; ivar < totNsnp; ivar++){
	//	fprintf(stderr, " %d ", allSnps[ivar]);
	//}

	//////////////////////////
	/// COLLATE LOCI: SORT ///
	/////////////////////////
	//fprintf(stderr, "now sorting:\n");
	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	//fprintf(stderr, "Sorted SNPs from pos %d to %d\n", allSnps[0], allSnps[totNsnp-1]);
	nunique = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{nunique++;}
	} // end for isnp
	//fprintf(stderr, "Found %d SNPs with values for at least one pop comparison...\n", nunique);

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

	//fprintf(stderr, "Allocating memory...\n");
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
		data->genpos[iComp] = calloc(nunique, sizeof(double));
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

	//fprintf(stderr, "Loading all component scores...\n");
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);
		//fprintf(stderr, "loading pop-pair from file: ");
		//fprintf(stderr, infilename);
		//fprintf(stderr, "\n");

		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");		
		fgets(H12_filename, line_size, inf);
		strtok(H12_filename, "\n");		
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");		
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);

		get_popPair_completeData(&data_sing, ihs_filename, delihh_filename, nsl_filename, H12_filename, iSAFE_filename, xpehh_filename, freqs_filename);
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
				data->nsl_normed[iComp][isnp] = data_sing.nsl_normed[jsnp];		 
				data->H12[iComp][isnp] = data_sing.H12[jsnp];		 
				data->H2H1[iComp][isnp] = data_sing.H2H1[jsnp];		 
				data->iSAFE[iComp][isnp] = data_sing.iSAFE[jsnp];		 

				jsnp++; //assert(jsnp<=data_sing.nsnps);
				if (jsnp >= data_sing.nsnps){break;}
			}
			else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; if (jsnp >= data_sing.nsnps){break;}}//assert(jsnp<=data_sing.nsnps);}
			//else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
		}// end for isnp loop
		free_popPair_data(&data_sing); 
	} // end for icomp
	free(allUniqueSnps);
	free(allSnps);
	free(newLine);
	//fprintf(stderr, "loaded multiple pop-pair comparisons to data object.\n");
} //end method
void get_popComp_anyData(popComp_data_multiple* data, int nComparisons, int argc, char *argv[]){
	/*This function loads in all SNPs, even those for which there is incomplete data.
	argv is all files for pop-pairs (each of which points to further component score files)*/
	const int line_size = 15000000; 
	popPair_data data_sing;
	FILE *inf=NULL;
	char *newLine;//, *token, *running;
	char infilename[512];
	int	isnp, jsnp, iComp, totNsnp, nunique; //thisPhysPos, itoken
	int *allSnps, *allUniqueSnps;
	char ihs_filename[528], delihh_filename[528], nsl_filename[528], xpehh_filename[528], freqs_filename[528];
	char h12_filename[528], iSAFE_filename[528];
	int numExtraArgs;

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

	////////////////////
	/// COLLATE LOCI ///
	////////////////////
	totNsnp = 0; //pass over first argument (run params)
	numExtraArgs = argc - nComparisons;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);
		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");
		fgets(h12_filename, line_size, inf);
		strtok(h12_filename, "\n");		
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");		
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);
		totNsnp += get_num_anyData(ihs_filename, delihh_filename, nsl_filename, h12_filename, iSAFE_filename, xpehh_filename, freqs_filename);
	} // end iComp
	//fprintf(stderr, "totNsnp: %d\n", totNsnp);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, "%s", argv[iComp + numExtraArgs]);
		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");
		fgets(h12_filename, line_size, inf);
		strtok(h12_filename, "\n");
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);
		get_popPair_anyData(&data_sing, ihs_filename, delihh_filename, nsl_filename, h12_filename, iSAFE_filename, xpehh_filename, freqs_filename); 
		//fprintf(stderr, "checkpoint: completed get_popPair_anyData\n");
		//fprintf(stderr, "%d\t%d\t%d\t%d\n", data_sing.physpos[0], data_sing.physpos[1], data_sing.physpos[data_sing.nsnps-2],data_sing.physpos[data_sing.nsnps-1]);


		for (jsnp = 0; jsnp < data_sing.nsnps; jsnp++){
			allSnps[isnp] = data_sing.physpos[jsnp];
			isnp +=1;
		}
		free_popPair_data(&data_sing); 
	}//end icomp

	//////////////////////////
	/// COLLATE LOCI: SORT ///	Problem here?
	/////////////////////////
	qsort(allSnps, totNsnp, sizeof(int), intcmp);

	nunique = count_unique(allSnps, totNsnp);
	//nunique = 0;
	//for (isnp = 0; isnp <= totNsnp-2; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
	//	if (allSnps[isnp] == allSnps[isnp+1]){continue;}
	//	else{nunique++;}
	//} //right?// end for isnp


	allUniqueSnps	= malloc(nunique * sizeof(int));
	jsnp = 0;
	for (isnp = 0; isnp <= totNsnp-2; isnp++){
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
	} // end for isnp
	allUniqueSnps[totNsnp-1] = allSnps[totNsnp-1]; 

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

		inf = fopen(infilename, "r");
		fgets(ihs_filename, line_size, inf);
		strtok(ihs_filename, "\n");
		fgets(delihh_filename, line_size, inf);
		strtok(delihh_filename, "\n");
		fgets(nsl_filename, line_size, inf);
		strtok(nsl_filename, "\n");
		fgets(h12_filename, line_size, inf);
		strtok(h12_filename, "\n");
		fgets(iSAFE_filename, line_size, inf);
		strtok(iSAFE_filename, "\n");		
		fgets(xpehh_filename, line_size, inf);
		strtok(xpehh_filename, "\n");
		fgets(freqs_filename, line_size, inf);
		strtok(freqs_filename, "\n");
		fclose(inf);

		get_popPair_anyData(&data_sing, ihs_filename, delihh_filename, nsl_filename, h12_filename, iSAFE_filename, xpehh_filename, freqs_filename);
		fprintf(stderr, "loaded pop-pair object with %d snps\n", data_sing.nsnps);
		//fprintf(stderr, "DEBUG: %d\n ", data_sing.physpos[data_sing.nsnps-1]);

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
	free(allSnps);
	free(newLine);
	//fprintf(stderr, "loaded multiple pop-pair comparisons to data object.\n");
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
	//fprintf(stderr, "%f\t", fst); //FOR DEBUG
	return -1 * log(1. - fst); //validate?
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
	//fprintf(stderr, "\n%d\n", data->physpos[0][isnp]);
	alt_daf = 0;
	for (iComp = 0; iComp < data->ncomp; iComp++){
		if (data->delDAF[iComp][isnp] !=0){
			thiscomp_deldaf = data->delDAF[iComp][isnp];
			thiscomp_daf_sel = data->daf_selpop[iComp][isnp];
			thiscomp_daf_altpop = thiscomp_daf_sel - thiscomp_deldaf; //VALIDATE
			alt_daf += thiscomp_daf_altpop;
			theseComp++;
			//fprintf(stderr, "thiscomp_deldaf %f thiscomp_daf_sel %f\n", thiscomp_deldaf, thiscomp_daf_sel);
			//fprintf(stderr, "alt: %f\n", thiscomp_daf_altpop);
		} // end if
	} // end for iComp
	ave_daf = alt_daf / (double)theseComp;
	deldaf = thiscomp_daf_sel - ave_daf;
	//fprintf(stderr, "made %d comparisons and found average outgroup DAF: %f\n", theseComp, ave_daf); //FOR DEBUG
	//fprintf(stderr, "selpop daf: %f ; delDAF: %f\n", thiscomp_daf_sel, deldaf); // FOR DEBUG
	return deldaf;
} //end function
