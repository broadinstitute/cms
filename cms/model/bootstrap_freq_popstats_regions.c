// last updated: 07.27.16   vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h" 

int main(int argc, char **argv){
    const int line_size = 15000000; 
    int chromosome, isnp, isamp, ns, nderiv, seqlen, startpos, endpos, itoken, nanc, npoly;
    char recomfilename[264], tpedfilename[264], regionfilename[264], outfilename[264];
    char *newLine, *token, *running;
    double pi;//, pi_sum=0.;
    FILE *outf=NULL, *inf=NULL;
	coal_data data;
    
	if (argc != 5) {
		fprintf(stderr, "Usage: bootstrap_freq_popstats_regions <tpedfilename> <recomfilename> <regionfilename> <outfilename>\n");
		exit(0);
	}

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	//initialize sums to zero
	//pi_sum = 0;
    seqlen = 0;
    npoly = 0;
	
    strcpy(tpedfilename, argv[1]);
    strcpy(recomfilename, argv[2]);
    strcpy(regionfilename, argv[3]);
    strcpy(outfilename, argv[4]);

    get_coal_data_tped_vers(&data, tpedfilename, recomfilename);

    inf = fopen(regionfilename, "r");
	assert(inf != NULL);
	if (inf == NULL) {fprintf(stderr, "Missing regions file: %s\n", regionfilename);}

    outf = fopen(outfilename, "w");
    assert(outf != NULL);
    fprintf(stderr, "Writing to file: %s\n", outfilename);

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
        fprintf(outf, "%d\n\n", endpos - startpos); //must print region len; mark off blocks with line break
	} // end regions

    free_coal_data(&data);
    fclose(outf);
    return 0;
} //end main