// {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
// last updated: 06.04.16   vitti@broadinstitute.org
//## JV-> orig: gcc -o bootstrap_freq_popstats_1kg_regions_thinned -O0 -ggdb3 -lm -Wall bootstrap_freq_popstats_1kg_regions_thinned.c coal_data_1kg_autosome_tped_thinned.c

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h"

int main(int argc, char **argv){
	coal_data data;
	FILE *outf=NULL, *inf=NULL;
	char filename[264], bedfilename[264], popstr[3], regionfilename[264];
	int chromosome, isnp, isamp, ns, nderiv, seqlen, startpos, endpos, itoken, nanc, npoly;
	double pi, pi_sum=0.;
	char *newLine, *token, *running;
	const int line_size = 15000000; 
    int nlimiting;
    
	if (argc != 5) {
		fprintf(stderr, "Usage: bootstrap_freq_popstats_1kg_regions_thinned <pop> <chrom> <preprocessed bedfile from regions folder> <nlimiting>\n");
		exit(0);
	}

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	//initialize sums to zero
	pi_sum = 0;
    seqlen = 0;
    npoly = 0;
	

	/************
	 OUTFILE PREP
	 ************/

	strcpy(popstr, argv[1]);
    chromosome = atoi(argv[2]);
	strcpy(bedfilename, argv[3]);
    nlimiting = atoi(argv[4]);

    strcpy(regionfilename, "regions/");
    strcat(regionfilename, bedfilename);
    strcat(regionfilename, ".txt");
	inf = fopen(regionfilename, "r");
	assert(inf != NULL);
	fprintf(stderr, "Getting regions from file: ");
	fprintf(stderr, regionfilename);
	fprintf(stderr, "\n");

	if (inf == NULL) {fprintf(stderr, "Missing regions file: "); fprintf(stderr, regionfilename);}

    get_coal_data_1kg_autosome_tped_thinned(&data, popstr, "strictMask", chromosome, "062215_2pm_nomultialleles", nlimiting);

    strcpy(filename,"bootstrap/thinned/");
    strcat(filename, popstr);
    strcat(filename, "_");
     strcat(filename, argv[4]);
    strcat(filename, "_");    
    strcat(filename, bedfilename);
    strcat(filename, "_freqs_by_snp.txt");
    outf = fopen(filename, "w");
    assert(outf != NULL);
    fprintf(stderr, "Writing to file: ");
    fprintf(stderr, filename);
    fprintf(stderr, "\n");

	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
			if (itoken == 0) {chromosome = atoi(token);}
			if (itoken == 1) {startpos = atoi(token);}
			if (itoken == 2) {endpos = atoi(token);}
		}
		
    	fprintf(stderr, "chrom: %d\t", chromosome);
		fprintf(stderr, "start: \t%d\t", startpos);
		fprintf(stderr, "end: \t%d\n", endpos);

		//augment seqlen
		seqlen += (endpos-startpos);

		
    	//loop over SNPs
    	for (isnp = 0; isnp < data.nsnp; isnp++) {
    		//bounds for this region
    		if (data.pos[isnp] < startpos){continue;}
    		if (data.pos[isnp] > endpos){break;}
    
    		nderiv = 0;
    		ns = 0;
    		for (isamp = 0; isamp < data.nsample; isamp++) {
    			if (data.genotype[isamp][isnp] != 2) {
    				ns++;
    				if (data.genotype[isamp][isnp] == 0) {nderiv++;}
    				}
    			}
    		assert(ns == data.nsample);
    		nanc = ns - nderiv;

    		if (nanc != 0 && nderiv != 0) {
                npoly++;
                pi = 2. * nderiv * nanc / ns / (ns-1);
                fprintf(outf, "%d\t%d\t%.9f\t%d\t%d\n", chromosome, data.pos[isnp], pi, nderiv, nanc);
    		}// end if polymorphic
    
        } // end loop over snps
        fprintf(outf, "%d\n\n", endpos - startpos);
            //must print region len
     //mark off blocks with line break
	} // end regions

free_coal_data(&data);
fclose(outf);
return 0;
}
