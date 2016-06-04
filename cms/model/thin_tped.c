// {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
// thins a tped down to nchroms given as input. samples without replacement 11.09.15 
// ORIG: gcc -o thin_tped -O0 -ggdb3 -lm -Wall thin_tped.c coal_data_1kg_autosome_tped.c
//## last updated: 06.04.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "coal_data_tped_vers.h"

int intcmp(const void *v1, const void *v2);
int intcmp(const void *v1, const void *v2){return (*(int *)v1 - *(int *)v2);}

/***********/
/***MAIN****/
/***********/

int main(int argc, char **argv){
  coal_data data;
  FILE *outf=NULL;
  char outfilename[528], popstr[3], chromstr[3], nchromstotakestr[3], genotypestring[2];
  char towritestring[5000];
  int chromosome, nChromsToTake, isamp, ihapindiv, r, ichrom, isnp;
  int *toTakeIndices;
  int duplicate; //0 FALSE 1 TRUE
  int thisindex;

  if (argc != 4) {
    fprintf(stderr, "Usage: thin_tped <pop> <chrom> <nChromsToTake>\n");
    exit(0);
  }

  strcpy(popstr, argv[1]);
  chromosome = atoi(argv[2]);
  sprintf(chromstr, "%d", chromosome);
  nChromsToTake = atoi(argv[3]);
  sprintf(nchromstotakestr, "%d", nChromsToTake);

  get_coal_data_1kg_autosome_tped(&data, popstr, "strictMask", chromosome, "062215_2pm_nomultialleles");

  ////////////////////////////////
  /// CHOOSE A SUBSET OF INDIVS //
  ////////////////////////////////

  fprintf(stderr, "sampling ");
  fprintf(stderr, nchromstotakestr);
  fprintf(stderr, " haploid individuals...\n");
  toTakeIndices = calloc(nChromsToTake, sizeof(int)); 
  assert (toTakeIndices != NULL);
  for (ichrom = 0; ichrom < nChromsToTake; ichrom++){//initialize all elements to negative values 
    toTakeIndices[ichrom] = -1;
  } 

  srand(time(NULL));
  for (ihapindiv = 0; ihapindiv < nChromsToTake; ihapindiv++){
    duplicate = 1; // start off assuming we have a duplicate
    while(duplicate == 1){
      //pick a random index
      r = rand() % data.nsample;
      //fprintf (stderr, "r: %d\n", r);
      //check if it's been taken
      duplicate = 0; //false
      for (ichrom = 0; ichrom < nChromsToTake; ichrom++){
        if (toTakeIndices[ichrom] == r){duplicate = 1; break;}
      } // end for ichrom loop
    } // end while duplicate loop
    toTakeIndices[ihapindiv] = r;
  } // end for ihapindiv loop

  ///////////////////////////////////
  /// WRITE SAMPLED CHROMS TO FILE //
  ///////////////////////////////////

  qsort(toTakeIndices, nChromsToTake, sizeof(int), intcmp);

  strcpy(outfilename, "/idi/sabeti-scratch/jvitti/cms_venv/tpeds/thinned/v5a_chr");
  strcat(outfilename, chromstr);
  strcat(outfilename, "_strictMask_");
  strcat(outfilename, popstr);
  strcat(outfilename, "_062215_2pm_nomultialleles");
  strcat(outfilename, "_thinned_");
  strcat(outfilename, nchromstotakestr);
  strcat(outfilename, ".tped");
  fprintf(stderr, "writing to: ");
  fprintf(stderr, outfilename);
  fprintf(stderr, "\n");
  outf = fopen(outfilename, "w");
  assert(outf != NULL);

  //for each snp in the tped, generate a new list of genotypes for our filtered file
  for (isnp = 0; isnp < data.nsnp; isnp++){

    //iterate over all indivs and only retain those in the sorted array, toTakeIndices
    sprintf(towritestring, "%d %d %.12f %d ", chromosome, data.pos[isnp], data.genloc[isnp], data.pos[isnp]);
    thisindex = 0;
    for (isamp = 0; isamp < data.nsample; isamp++){
      //retaining this individual
      if (isamp == toTakeIndices[thisindex]){
        sprintf(genotypestring, "%d", data.genotype[isamp][isnp]);
        strcat(towritestring, genotypestring);

        if (thisindex != nChromsToTake-1){strcat(towritestring, " ");} //is not last:
        //advance thisindex
        thisindex +=1;
      } // end match
    } // end loop over all indivs
    fprintf(outf, towritestring);
    fprintf(outf, "\n");

  } // end snp loop

  fclose(outf);
  free_coal_data(&data);
  return 0;
} // end main