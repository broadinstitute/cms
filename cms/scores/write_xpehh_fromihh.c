// {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
// for a given pop pair, parse genome-wide iHH to calculate and record XP-EHH files.
// lifted from cms_venv : // gcc -o write_xpehh_fromihh -O0 -ggdb3 -lm -Wall write_xpehh_fromihh.c pop_ihh_data.c
//#JV: must fix file locations globally
// last updated: 06.04.16 	vitti@broadinstitute.org

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
	char chromstr[3], pop1[3], pop2[3], n1[3], n2[3];
	char infilename1[264], infilename2[264], outfilename[264];
	char *newLine, *running, *token;
	int line_size = 5000;
	int itoken, inf1_snp, inf2_snp, pos1, pos2;
	double gdpos1, gdpos2, ihh1, ihh2, p1, p2, xpehh;
	int index_1, index_2, ichrom;

	if (argc != 3) {
		fprintf(stderr, "Usage: ./write_xpehh_fromihh <pop1> <pop2>\n");
		exit(0);
	}

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	strcpy(pop1, argv[1]);
	strcpy(pop2, argv[2]);

	strcpy(outfilename, "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/from_ihh/merged/v5a_autosomes_strictMask_062215_2pm_nomultialleles_merged_121915_");
	strcat(outfilename, pop1);
	strcat(outfilename, "_");
	strcat(outfilename, pop2);
	strcat(outfilename, ".xpehh");
	fprintf(stderr, "writing to: ");
	fprintf(stderr, outfilename);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL);

	////////////////////////
	// LOOP OVER AUTOSOMES//
	////////////////////////

	for (ichrom = 1; ichrom < 23; ichrom++){
		sprintf(chromstr, "%d", ichrom);
		if (ichrom <= 20){ // only needed to redo calculation for chr1-20 (>100cM)
		strcpy(infilename1, "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/merged/v5a_chr");
		strcat(infilename1, chromstr);
		strcat(infilename1, "_strictMask_");
		strcpy(infilename2, infilename1);
		strcat(infilename1, pop1);
		strcat(infilename2, pop2);
		strcat(infilename1, "_062215_2pm_nomultialleles_merged_121915.xpehh.out.ihh");
		strcat(infilename2, "_062215_2pm_nomultialleles_merged_121915.xpehh.out.ihh");
		}
		else{ 
		strcpy(infilename1, "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/xp_run/v5a_chr");
		strcat(infilename1, chromstr);
		strcat(infilename1, "_strictMask_");
		strcpy(infilename2, infilename1);
		strcat(infilename1, pop1);
		strcat(infilename2, pop2);
		strcat(infilename1, "_062215_2pm_nomultialleles_xpehh.out.ihh");
		strcat(infilename2, "_062215_2pm_nomultialleles_xpehh.out.ihh");
		}
        fprintf(stderr, "getting data from: ");
        fprintf(stderr, "\n");
        fprintf(stderr, infilename1);
        fprintf(stderr, "\n");
        fprintf(stderr, infilename2);
        fprintf(stderr, "\n");
		
		////////////////
		// LOAD SCORES//
		////////////////

		get_pop_ihh_data(&data1, infilename1);
		get_pop_ihh_data(&data2, infilename2);

		inf1_snp = data1.nsnps;
		inf2_snp = data2.nsnps;
        fprintf(stderr, "%d\t%d\n", inf1_snp, inf2_snp);
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
			//fprintf(stderr, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", pos1, pos1, gdpos1, ihh1, p1, ihh2, p2, xpehh);
			//calculate xpehh and write to file
			if (pos1 == pos2){				
				if (ihh1 == 0){ihh1 = 1e-08;}
				if (ihh2 == 0){ihh2 = 1e-08;}
				xpehh = log(ihh1/ihh2) * -1;
				fprintf(outf, chromstr);
				fprintf(outf, "_");
				fprintf(outf, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", pos1, pos1, gdpos1, ihh1, p1, ihh2, p2, xpehh);
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
	}

	return 0;
} // end main