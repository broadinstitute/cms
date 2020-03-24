// data structure to facilitate parsing of pop ihh files
// last updated: 10.21.16   vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
//#include <zlib.h>
#include "pop_ihh_data.h"

void get_pop_ihh_data(ihh_data* data, char filename[]) {
  const int line_size = 15000000;
  FILE *inf=NULL;
  char *newLine, *token, *running;
  int isnp, itoken;
    
  newLine = malloc((line_size+1) * sizeof(char));
  assert(newLine != NULL); 

  data->nsnps = 0;
  data->pos = NULL; 
  data->genpos = NULL;
  data->freq = NULL;
  data->ihh = NULL;

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
  data->freq = malloc(data->nsnps * sizeof(double*)); assert(data->freq != NULL);
  data->ihh = malloc(data->nsnps * sizeof(double*)); assert(data->ihh != NULL);

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
      else if (itoken == 1) {
        data->genpos[isnp] = atof(token);
      }
      else if (itoken == 2) {
        data->freq[isnp] = atof(token);
      }      
      else if (itoken == 3) {
        data->ihh[isnp] = atof(token);
      }
    } // END for running=newLine
    isnp++;
  } //END while(fgets(newLine))
  
  fclose(inf);
  free(newLine);
}

void free_pop_ihh_data(ihh_data* data) {
  if (data == NULL) {return;}
  free(data->pos);
  free(data->genpos);
  free(data->freq);
  free(data->ihh);
  data->nsnps = 0;
}