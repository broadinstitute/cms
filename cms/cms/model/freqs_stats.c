//  alternate/replacement to calc_fst_deldaf; pass information along in a way that anticipate the needs of cms_data structures
//  06.20.2017  	vitti@broadinstitute.org   //for now, preserve both.

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include "coal_data_tped_vers.h"

int main(int argc, char **argv) {
	coal_data data;
	char filename[264];
	FILE *outf=NULL;
	char inTped1[264], inTped2[264], inRecomfile[264];
	int *nall0[2]={NULL}, *nall1[2]={NULL};
	int isamp, isnp, nsnp;
	int nai[2], naj[2], na_both[2];
	int ni, nj;
	double p[2];
	double pmean, nic, njc, nc, msp, msg, num, denom;
	double fst_sum=0; //per chrom(/file)
	int nfst=0; // per chrom(/file)
	double fst, delDAF; //per site
	double T; // transformation on Fst cf Cavalli-Sforza 1969
	
	if (argc != 5) {
		fprintf(stderr, "Usage: freqs_stats <inTped1> <inTped2> <recomfilename> <writefilename>\n");
		exit(0);
	}
	strcpy(inTped1, argv[1]);
	strcpy(inTped2, argv[2]);
	strcpy(inRecomfile, argv[3]);
	strcpy(filename, argv[4]);
	outf = fopen(filename, "w");
	assert(outf != NULL);
	fprintf(outf, "physPos\tgenPos\tpopDAF\tdelDAF\tFst\tT\n"); //header 

	////////////////////////////
	// COUNT ALLELES FOR pop0 //
	////////////////////////////
	get_coal_data_tped_vers_gz(&data, inTped1, inRecomfile);   
	nall0[0] = calloc(data.nsnp, sizeof(int));
	nall1[0] = calloc(data.nsnp, sizeof(int));
	for (isnp = 0; isnp < data.nsnp; isnp++) {
		for (isamp = 0; isamp < data.nsample; isamp++) {
			if (data.genotype[isamp][isnp] == 0) {nall0[0][isnp]++;}
			else if (data.genotype[isamp][isnp] == 1) {nall1[0][isnp]++;}
		} //end for isamp
	} //end for isnp
	//fprintf(stderr, "%d snps, %d samples.\n", data.nsnp, data.nsample);
	free_coal_data(&data);
	
	////////////////////////////
	// COUNT ALLELES FOR pop1 //
	////////////////////////////
	get_coal_data_tped_vers_gz(&data, inTped2, inRecomfile);
	nall0[1] = calloc(data.nsnp, sizeof(int));
	nall1[1] = calloc(data.nsnp, sizeof(int));
	for (isnp = 0; isnp < data.nsnp; isnp++) {
		for (isamp = 0; isamp < data.nsample; isamp++) {
			if (data.genotype[isamp][isnp] == 0) {nall0[1][isnp]++;}
			else if (data.genotype[isamp][isnp] == 1) {nall1[1][isnp]++;}
		}   //end for isamp
	} //end for isnp
	//fprintf(stderr, "%d snps, %d samples.\n", data.nsnp, data.nsample);	

	///////////////////////////////
	// COMPARE FREQS, SNP BY SNP //
	///////////////////////////////
	nsnp = data.nsnp;
	for (isnp = 0; isnp < nsnp; isnp++) {
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
		//fprintf(stderr, "frequency of 0 in pop 1: %f\n", p[1]);
		
		//if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) {continue;} 
        //PREVIOUS: if the same allele is fixed in both populations, we skip it. 
        //Here, we include it, to facilitate loading data downstream for composite calcs.

        if (ni > 0 && nj > 0) {
			/////////////////////////
			// Weir-Hill estimator //
			////////////////////////
			//fprintf(stderr, "WHE\t");
			pmean = (ni * p[0] + nj * p[1]) / (ni + nj);
			nic = ni - (double) ni * ni / (ni + nj);
			njc = nj - (double) nj * nj / (ni + nj);
			nc = nic + njc;
			msp = ni * (p[0] - pmean) * (p[0] - pmean) + nj * (p[1] - pmean) * (p[1] - pmean);
			msg = (ni * p[0] * (1. - p[0]) + nj * p[1] * (1. - p[1])) / (ni - 1 + nj - 1);
			num = msp - msg;
			denom = msp + (nc - 1) * msg;
			
			delDAF = p[0] - p[1];//allele 0 is encoded as derived; we treat the first pop as the putative selected pop
			if (denom != 0) {
				fst = (num/denom);
				fst_sum += fst;
				nfst++;
			} // end if denom != 0
			else{fst = NAN;}

			T= -1 * log(1. - fst);

			fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\n", data.pos[isnp], data.genloc[isnp], p[0], delDAF, fst, T); // "physPos\tgenPos\tpopDAF\tdelDAF\tFst\tT\n
		} //end Weir Hill estimator
	} //end snp loop
	free_coal_data(&data);
	//fprintf(stderr, "fst_sum: %f\n", fst_sum);
	//fprintf(stderr, "nfst: %d\n", nfst);
	//fprintf(stderr, "chrom ave Fst: %.8f\n", (fst_sum/nfst));
	fprintf(stderr, "wrote to %s\n", filename);
	return 0;
} //end main
