// {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
// last updated: 06.04.16 	vitti@broadinstitute.org
//#JV- orig: gcc -o bootstrap_fst_1kg_regions_thinned -O0 -ggdb3 -lm -Wall bootstrap_fst_1kg_regions_thinned.c coal_data_1kg_autosome_tped_thinned.c

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h"
#define POPSIZE 2

int main(int argc, char **argv) {
	coal_data data;
	char filename[264];
	FILE *outf=NULL, *inf=NULL;
	char pop0[3], pop1[3];
	char mask[25];
	char bedfilename[264], regionfilename[528];
	int chrom;
	int *nall0[POPSIZE]={NULL}, *nall1[POPSIZE]={NULL};
	int isamp, isnp, ichrom, nsnp;
	int nai[2], naj[2], na_both[2];
	int samp0, samp1;
	char runtime[264];
	int ni, nj;
	double p[2];
	double pmean, nic, njc, nc, msp, msg, num, denom;
	double fst_sum; //per chrom
	int nfst; // per chrom
	char *newLine, *token, *running, *token2;
	const int line_size = 15000000; 
	int startpos, endpos, seqlen, chromosome, itoken;
    int nlimiting;	
	if (argc != 6) {
		fprintf(stderr, "Usage: bootstrap_fst_1kg_regions_thinned <pop0> <pop1> <chrom> <preprocessed bedfile from regions folder> <nlimiting>\n");
		exit(0);
	}
	
	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	strcpy(pop0, argv[1]);
	strcpy(pop1, argv[2]);
	chromosome = atoi(argv[3]);
	strcpy(bedfilename, argv[4]);
    nlimiting = atoi(argv[5]);

    //count alleles for pop0
	get_coal_data_1kg_autosome_tped_thinned(&data, pop0, "strictMask", chromosome, "062215_2pm_nomultialleles", nlimiting);
	nall0[0] = calloc(data.nsnp, sizeof(int));
	nall1[0] = calloc(data.nsnp, sizeof(int));
	for (isnp = 0; isnp < data.nsnp; isnp++) {
		for (isamp = 0; isamp < data.nsample; isamp++) {
			if (data.genotype[isamp][isnp] == 0) {nall0[0][isnp]++;}
			else if (data.genotype[isamp][isnp] == 1) {nall1[0][isnp]++;}
		}
	}
	fprintf(stderr, "%d snps, %d samples.\n", data.nsnp, data.nsample);
	free_coal_data(&data);	
	//count alleles for pop1
	get_coal_data_1kg_autosome_tped_thinned(&data, pop1, "strictMask", chromosome, "062215_2pm_nomultialleles", nlimiting);	
	nall0[1] = calloc(data.nsnp, sizeof(int));
	nall1[1] = calloc(data.nsnp, sizeof(int));
	for (isnp = 0; isnp < data.nsnp; isnp++) {
		for (isamp = 0; isamp < data.nsample; isamp++) {
			if (data.genotype[isamp][isnp] == 0) {nall0[1][isnp]++;}
			else if (data.genotype[isamp][isnp] == 1) {nall1[1][isnp]++;}
		}
	}
	fprintf(stderr, "%d snps, %d samples.\n", data.nsnp, data.nsample);
	nsnp = data.nsnp;

	fst_sum = 0;
	nfst = 0;
	
	strcpy(regionfilename, "regions/");
	strcat(regionfilename, bedfilename);
	strcat(regionfilename, ".txt");
	inf = fopen(regionfilename, "r");
	assert(inf != NULL);
    if (inf == NULL) {fprintf(stderr, "Missing regions file: "); fprintf(stderr, regionfilename);}

	fprintf(stderr, "Getting information from file: ");
	fprintf(stderr, regionfilename);
	fprintf(stderr, "\n");

	strcpy(filename,"bootstrap/thinned/");
	strcat(filename, pop0);
	strcat(filename, "_");
	strcat(filename, pop1);
	strcat(filename, "_");
	strcat(filename, argv[5]);
	strcat(filename, "_");
	strcat(filename, bedfilename);
	strcat(filename, "_fst_by_snp.txt");
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

		
		//DEBUGGING
		for (isnp = 0; isnp < nsnp; isnp++) {
			//bounds for this region
			if (data.pos[isnp] < startpos){continue;}
			if (data.pos[isnp] > endpos){break;}

			//fprintf(stderr, "snp: %d\n", isnp);
			if (nall0[0] == NULL) {continue;}
			nai[0] = nall0[0][isnp]; //number of 0 alleles in pop 0
			//fprintf(stderr, "number of 0 alleles in pop 0: %d\n", nai[0]);
			nai[1] = nall1[0][isnp]; //number of 1 alleles in pop 0
			//fprintf(stderr, "number of 1 alleles in pop 0: %d\n", nai[1]);
			ni = nai[0] + nai[1]; //number of alleles, pop 0
			//fprintf(stderr, "number of alleles in pop 0: %d\n", ni);
			p[0] = (double) nai[0] / ni; //frequency of allele 0 in pop 0
			//fprintf(stderr, "frequency of 0 in pop 0: %f\n", p[0]);
			
			if (nall0[1] == NULL) {continue;}
			naj[0] = nall0[1][isnp]; //number of 0 alleles in pop 1
			//fprintf(stderr, "number of 0 alleles in pop 1: %d\n", naj[0]);
			naj[1] = nall1[1][isnp]; //number of 1 alleles in pop 1
			//fprintf(stderr, "number of 1 alleles in pop 1: %d\n", naj[1]);
			na_both[0] = nai[0] + naj[0]; //combined number of 0 alleles
			//fprintf(stderr, "number of 0 alleles total: %d\n", na_both[0]);
			na_both[1] = nai[1] + naj[1]; //combined number of 1 alleles
			//fprintf(stderr, "number of 1 alleles total: %d\n", na_both[1]);
			nj = naj[0] + naj[1]; //number of alleles, pop 1
			//fprintf(stderr, "number of alleles in pop 1: %d\n", nj);
			p[1] = (double) naj[0] / nj; //frequency of allele 0 in pop 1
		   // fprintf(stderr, "frequency of 0 in pop 1: %f\n", p[1]);
			
			if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) {continue;}
			
			if (ni > 0 && nj > 0) {
				//fprintf(stderr, "WHE\t");
				// Weir-Hill estimator
				pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
				nic = ni - (double) ni * ni / (ni + nj);
				njc = nj - (double) nj * nj / (ni + nj);
				nc = nic + njc;
				msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
				msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
				num = msp - msg;
				denom = msp + (nc - 1) * msg;
				
				if (denom != 0) {
					fprintf(outf, "%d\t%d\t%f\n", chromosome, data.pos[isnp], (num/denom));
					//fst_sum += (num/denom);
					//nfst++;
				} // end if denom != 0
			} //end Weir Hill estimator
		} //end snp loop

		fprintf(outf, "\n"); //separate regions with newline
	} // end region
	free_coal_data(&data);
	return 0;
} // end main