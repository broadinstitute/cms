// {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
// parse out iHH from XP-EHH files. Configured for autosomes
// original gcc -o rewrite_xpehh_files_merged -O0 -ggdb3 -lm -Wall rewrite_xpehh_files_merged.c 12.27.15
// JV: FILELOCS
// last updated 06.04.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char **argv) {
    FILE *outf1=NULL;
    FILE *outf2=NULL;    
    FILE *inf=NULL;
    char chromstr[3], pop1[3], pop2[3];
    char infilename[264];
    char outfilename1[264];
    char outfilename2[264];
    char *newLine, *running, *token;
    int line_size = 5000;
    int itoken;
    double gdpos, p1, p2, ihh1, ihh2, xpehh;
    int pos, name;

    if (argc != 4) {
        fprintf(stderr, "Usage: ./rewrite_xpehh_files_merged <chrom> <pop1> <pop2>\n");
        exit(0);
    }

    strcpy(chromstr, argv[1]);
    strcpy(pop1, argv[2]);
    strcpy(pop2, argv[3]);

    strcpy(infilename, "/idi/sabeti-scratch/jvitti/cms_venv/xpehh/merged/v5a_chr");
    strcat(infilename, chromstr);
    strcat(infilename, "_strictMask_");
    strcpy(outfilename1, infilename);
    strcpy(outfilename2, infilename);
    strcat(infilename, pop1);
    strcat(infilename, "_");
    strcat(infilename, pop2);
    strcat(infilename, "_062215_2pm_nomultialleles_merged_121915.xpehh.out");
    strcat(outfilename1, pop1);
    strcat(outfilename1, "_062215_2pm_nomultialleles_merged_121915.xpehh.out.ihh");
    strcat(outfilename2, pop2);
    strcat(outfilename2, "_062215_2pm_nomultialleles_merged_121915.xpehh.out.ihh");

    inf = fopen(infilename, "r");
    fprintf(stderr, "reading in from: ");
    fprintf(stderr, infilename);
    fprintf(stderr, "\n");
    assert(inf != NULL);

    outf1 = fopen(outfilename1, "w");
    fprintf(stderr, "\twriting to: ");
    fprintf(stderr, outfilename1);
    fprintf(stderr, "\n");
    assert(outf1 != NULL);

    outf2 = fopen(outfilename2, "w");
    fprintf(stderr, "\twriting to: ");
    fprintf(stderr, outfilename2);
    fprintf(stderr, "\n");
    assert(outf2 != NULL);

	newLine = malloc((line_size+1) * sizeof(char));
  	assert(newLine != NULL);

    /*********
    PARSE FILE
    **********/

    fgets(newLine, line_size, inf); //strip header
    while (fgets(newLine, line_size, inf) != NULL){
    	for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
    		if (itoken == 0) {pos = atoi(token);}
    		else if (itoken == 1) {name = atoi(token);}
            else if (itoken == 2) {gdpos = atof(token);}            
			else if (itoken == 3) {p1 = atof(token);}
    		else if (itoken == 4) {ihh1 = atof(token);}
    		else if (itoken == 5) {p2 = atof(token);}
            else if (itoken == 6) {ihh2 = atof(token);}
            else if (itoken == 7) {xpehh = atof(token);}
    	} //end for(running=newLine)

        fprintf(outf1, "%d\t%f\t%f\t%f\n", pos, gdpos, p1, ihh1);
        fprintf(outf2, "%d\t%f\t%f\t%f\n", pos, gdpos, p2, ihh2);

	} //end while(fgets(newLine)) --- i.e., on to the next SNP
    fclose(inf);
    fclose(outf1);
    fclose(outf2);
    return 0;
} // end main