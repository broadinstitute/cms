// for a given pop pair, parse genome-wide iHH to calculate and record XP-EHH files.
// last updated: 06.12.17 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "pop_ihh_data.h"

int main(int argc, char **argv) {
	ihh_data data1, data2;
	FILE *outf=NULL;
	char n1[3], n2[3];
	char infilename1[264], infilename2[264], outfilename[264];
	char *newLine, *running, *token;
	int line_size = 5000;
	int itoken, inf1_snp, inf2_snp, pos1, pos2;
	double gdpos1, gdpos2, ihh1, ihh2, p1, p2, xpehh;
	int index_1, index_2, ichrom;

	if (argc != 4) {
		fprintf(stderr, "Usage: ./write_xpehh_fromihh <infilename1> <infilename2> <outfilename>\n");
		exit(0);
	}

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 
	strcpy(infilename1,argv[1]);
	strcpy(infilename2,argv[2]);
	strcpy(outfilename,argv[3]);
	fprintf(stderr, "writing to: %s\n", outfilename);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL);

    fprintf(stderr, "getting data from: %s and %s\n", infilename1, infilename2);
	
	////////////////
	// LOAD SCORES//
	////////////////

	get_pop_ihh_data(&data1, infilename1);
	get_pop_ihh_data(&data2, infilename2);

	inf1_snp = data1.nsnps;
	inf2_snp = data2.nsnps;
	index_1 = 0;
	index_2 = 0;

	/////////////////
	// COLLATE SNPS//
	/////////////////

	while (index_1 < inf1_snp && index_2 < inf2_snp){
        pos1 = data1.pos[index_1];
        pos2 = data2.pos[index_2];
        ihh1 = data1.ihh[index_1];
        ihh2 = data2.ihh[index_2];
        p1 = data1.freq[index_1];
        p2 = data2.freq[index_2];
        gdpos1 = data1.genpos[index_1];
        gdpos2 = data1.genpos[index_2];
 		//calculate xpehh and write to file
		if (pos1 == pos2){				
			if (ihh1 == 0){ihh1 = 1e-08;}
			if (ihh2 == 0){ihh2 = 1e-08;}
			xpehh = log(ihh1/ihh2); 
			fprintf(outf, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", pos1, pos1, gdpos1, p1, ihh1, p2, ihh2, xpehh);
		     index_1++;
             index_2++;
        }
		//advance indices
		if (pos1 < pos2){
            index_1++;
            //fprintf(stderr, "%d\t%d\n", pos1, pos2);
            }
        if (pos2 < pos1){
            index_2++;
            //fprintf(stderr, "%d\t%d\n", pos1, pos2);
        }			
	} // end while loop

	free_pop_ihh_data(&data1);
	free_pop_ihh_data(&data2);
	

	return 0;
} // end main