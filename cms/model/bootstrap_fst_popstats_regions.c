// last updated: 07.14.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h"

int main(int argc, char **argv) {
	const int line_size = 15000000; 
	int nai[2], naj[2], na_both[2], ni, nj, startpos, endpos, seqlen, chromosome, itoken, isamp, isnp, jsnp;
	int *nall0[2]={NULL}, *nall1[2]={NULL}; // count 'zero' alleles; count 'one' alleles
	double p[2], pmean, nic, njc, nc, msp, msg, num, denom;
	char recomfilename[264], tped1filename[264], tped2filename[264], regionfilename[264], outfilename[264];
	char *newLine, *token, *running;
	FILE *outf=NULL, *inf=NULL;
	coal_data data_i, data_j;

	if (argc != 6) {
		fprintf(stderr, "Usage: bootstrap_fst_popstats_regions <tped1filename> <tped2filename> <recomfilename> <regionfilename> <outfilename>\n\n");
		exit(0);
	}
	
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	strcpy(tped1filename, argv[1]);
	strcpy(tped2filename, argv[2]);
	strcpy(recomfilename, argv[3]);
	strcpy(regionfilename, argv[4]);
	strcpy(outfilename, argv[5]);

	//count alleles for pop_i
	get_coal_data_tped_vers(&data_i, tped1filename, recomfilename);
	nall0[0] = calloc(data_i.nsnp, sizeof(int));
	nall1[0] = calloc(data_i.nsnp, sizeof(int));
	for (isnp = 0; isnp < data_i.nsnp; isnp++) {
		for (isamp = 0; isamp < data_i.nsample; isamp++) {
			if (data_i.genotype[isamp][isnp] == 0) {nall0[0][isnp]++;}
			else if (data_i.genotype[isamp][isnp] == 1) {nall1[0][isnp]++;}
		}
	}
	fprintf(stderr, "%d snps, %d samples.\n", data_i.nsnp, data_i.nsample);	

	//count alleles for pop_j
	get_coal_data_tped_vers(&data_j, tped2filename, recomfilename);
	nall0[1] = calloc(data_j.nsnp, sizeof(int));
	nall1[1] = calloc(data_j.nsnp, sizeof(int));
	for (isnp = 0; isnp < data_j.nsnp; isnp++) {
		for (isamp = 0; isamp < data_j.nsample; isamp++) {
			if (data_j.genotype[isamp][isnp] == 0) {nall0[1][isnp]++;}
			else if (data_j.genotype[isamp][isnp] == 1) {nall1[1][isnp]++;}
		}
	}
	fprintf(stderr, "%d snps, %d samples.\n", data_j.nsnp, data_j.nsample);
	
	inf = fopen(regionfilename, "r");
	assert(inf != NULL);
	if (inf == NULL) {fprintf(stderr, "Missing regions file: %s\n", regionfilename);}
	seqlen = 0;

	outf = fopen(outfilename, "w");
	assert(outf != NULL);

	isnp = 0;
	jsnp = 0;

	//////////////////////////
	// ITERATE OVER REGIONS //
	//////////////////////////

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

		//if not in region, advance 
		while (data_i.pos[isnp] < startpos || data_j.pos[jsnp] < startpos){
			if (data_i.pos[isnp] < data_j.pos[jsnp]){isnp++;}
			if (data_i.pos[isnp] > data_j.pos[jsnp]){jsnp++;}
			if (data_i.pos[isnp] == data_j.pos[jsnp]){isnp++; jsnp++;}
		}

		//advance snp indices, stay within region
		while ((data_i.pos[isnp] >= startpos) && (data_i.pos[isnp] <= endpos) && (data_j.pos[jsnp] >= startpos) && (data_j.pos[jsnp] <= endpos)){
			//found match
			if (data_i.pos[isnp] == data_j.pos[jsnp]){
				if (nall0[0] == NULL) {continue;}
				nai[0] = nall0[0][isnp]; //number of 0 alleles in pop_i
				nai[1] = nall1[0][isnp]; //number of 1 alleles in pop_i
				ni = nai[0] + nai[1]; //number of alleles, pop_i
				p[0] = (double) nai[0] / ni; //frequency of allele 0 in pop_i
				
				if (nall0[1] == NULL) {continue;}
				naj[0] = nall0[1][jsnp]; //number of 0 alleles in pop_j
				naj[1] = nall1[1][jsnp]; //number of 1 alleles in pop_j
				nj = naj[0] + naj[1]; //number of alleles, pop_j
				p[1] = (double) naj[0] / nj; //frequency of allele 0 in pop_j

				na_both[0] = nai[0] + naj[0]; //combined number of 0 alleles
				na_both[1] = nai[1] + naj[1]; //combined number of 1 alleles

				//fixed in both pop
				if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) {isnp++; jsnp++; continue;}
				//not fixed
				if (ni > 0 && nj > 0) {
					// Weir-Hill estimator
					pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
					nic = ni - (double) ni * ni / (ni + nj);
					njc = nj - (double) nj * nj / (ni + nj);
					nc = nic + njc;
					msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
					msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
					num = msp - msg;
					denom = msp + (nc - 1) * msg;
					if (denom != 0) {fprintf(outf, "%d\t%d\t%f\n", chromosome, data_i.pos[isnp], (num/denom));} 
				} //end Weir Hill estimator
				isnp++;
				jsnp++;
			} // end found match
			//if no match, advance the lower index
			if (data_i.pos[isnp] < data_j.pos[jsnp]){isnp++;}
			if (data_i.pos[isnp] > data_j.pos[jsnp]){jsnp++;}
		} // end while in region
		fprintf(outf, "\n"); //separate regions with newline
	} // end fgets newline
	fprintf(stderr, "Wrote to file: %s\n", outfilename);
	free_coal_data(&data_i);
	free_coal_data(&data_j);
	return 0;
} // end main