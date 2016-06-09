// {{POP GEN DATA --> COMPONENT SCORES} --> COMPOSITE SCORES}
// lifted: gcc -o TESTER -O0 -ggdb3 -lm -Wall TESTER.c coal_data_1kg_autosome_tped.c 01.14.16
// 06.04.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h"

int main(int argc, char **argv) {
    coal_data data;
    char filename[264];
    FILE *outf=NULL;
    char pop0[3], pop1[3];
    char mask[25];
    char chromstr[3];
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
    //double fst_complete_sum; //across all sites
    //int nfst_all; //across all
    int nfst; // per chrom
    double fst, delDAF; //per site
    int nSkip; 
    
    if (argc != 6) {
        fprintf(stderr, "Usage: calc_Fst <pop0> <pop1> <chrom> <runtime> <mask>\n");
        exit(0);
    }
    
    strcpy(pop0, argv[1]);
    strcpy(pop1, argv[2]);
    strcpy(chromstr, argv[3]);
    strcpy(runtime, argv[4]);
    strcpy(mask, argv[5]);

    strcat(filename, "TEST_TEST.txt");
    nSkip = 0;

    chrom = atoi(chromstr);
    strcpy(filename, "/idi/sabeti-scratch/jvitti/cms_venv/fst/chr");
    strcat(filename, chromstr);
    strcat(filename, "/");
    strcat(filename, pop0);
    strcat(filename, "_");
    strcat(filename, pop1);
    strcat(filename, "_chr");
    strcat(filename, chromstr);
    strcat(filename, "_");
    strcat(filename, runtime);
    strcat(filename, "_");
    strcat(filename, mask);
    strcat(filename, "_fst_daf.txt"); 
    outf = fopen(filename, "w");
    assert(outf != NULL);
    
    //nfst_all = 0;
    //fst_complete_sum = 0;

    fst_sum = 0;
    nfst = 0;
    //count alleles for pop0
    get_coal_data_1kg_autosome_tped(&data, pop0, mask, chrom, runtime);
    //if (data.nsample == 0) {continue;}
    
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
    get_coal_data_1kg_autosome_tped(&data, pop1, mask, chrom, runtime);
    //if (data.nsample == 0) {continue;}
    
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
    
    fprintf(outf, "pos\tfst\tdelDAF\n");

    //DEBUGGING
    for (isnp = 0; isnp < 8000; isnp++) {
        fprintf(stderr, "snp: %d\n", isnp);
        if (nall0[0] == NULL) {continue;}
        nai[0] = nall0[0][isnp]; //number of 0 alleles in pop 0
        fprintf(stderr, "number of 0 alleles in pop 0: %d\n", nai[0]);
        nai[1] = nall1[0][isnp]; //number of 1 alleles in pop 0
        fprintf(stderr, "number of 1 alleles in pop 0: %d\n", nai[1]);
        ni = nai[0] + nai[1]; //number of alleles, pop 0
        fprintf(stderr, "number of alleles in pop 0: %d\n", ni);
        p[0] = (double) nai[0] / ni; //frequency of allele 0 in pop 0
        fprintf(stderr, "frequency of 0 in pop 0: %f\n", p[0]);
        
        if (nall0[1] == NULL) {continue;}
        naj[0] = nall0[1][isnp]; //number of 0 alleles in pop 1
        fprintf(stderr, "number of 0 alleles in pop 1: %d\n", naj[0]);
        naj[1] = nall1[1][isnp]; //number of 1 alleles in pop 1
        fprintf(stderr, "number of 1 alleles in pop 1: %d\n", naj[1]);
        na_both[0] = nai[0] + naj[0]; //combined number of 0 alleles
        fprintf(stderr, "number of 0 alleles total: %d\n", na_both[0]);
        na_both[1] = nai[1] + naj[1]; //combined number of 1 alleles
        fprintf(stderr, "number of 1 alleles total: %d\n", na_both[1]);
        nj = naj[0] + naj[1]; //number of alleles, pop 1
        fprintf(stderr, "number of alleles in pop 1: %d\n", nj);
        p[1] = (double) naj[0] / nj; //frequency of allele 0 in pop 1
        fprintf(stderr, "frequency of 0 in pop 1: %f\n", p[1]);
        
        if ((nai[0] == 0 && naj[0] == 0) || (nai[1] == 0 && naj[1] == 0)) {nSkip++; continue;}
        
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
            
            delDAF = p[0] - p[1];//allele 0 is encoded as derived; we treat the first pop as the putative selected pop

            if (denom != 0) {
                fst = (num/denom);
                fst_sum += fst;

                nfst++;
            } // end if denom != 0
            else{fst = NAN;}
            fprintf(stderr, "%d\t%f\t%f\n", data.pos[isnp],fst,delDAF);
            fprintf(outf, "%d\t%f\t%f\n", data.pos[isnp],fst,delDAF);
        } //end Weir Hill estimator

    } //end snp loop

    free_coal_data(&data);

    fprintf(stderr, "skipped on: %d\n", nSkip);
    //fprintf(stderr, "fst_sum: %f\n", fst_sum);
    //fprintf(stderr, "nfst: %d\n", nfst);
    fprintf(outf, "chrom ave Fst: %.8f\n", (fst_sum/nfst));
    return 0;
}


