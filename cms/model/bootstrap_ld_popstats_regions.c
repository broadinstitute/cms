// {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
// COMPILE: gcc -o bootstrap_ld_popstats_regions -O0 -ggdb3 -lm -Wall bootstrap_ld_popstats_regions.c coal_data_tped_vers.c
// last updated: 06.07.16 	vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "coal_data_tped_vers.h" 

int main(int argc, char **argv){
	const double maxgendist = .017;
	const int maxdist = 70000; //only calculate LD for 70kb
	const int min_minor = 3; //don't calculate LD at singletons or doubletons
	//const double min_freq = .2; //only look at P(D'=1) above this freq // NB must filter downstream; 11.2015 bootstrap program does r2 and d' at once
	const int line_size = 15000000; 
	int isnp, jsnp, isamp, ns, nderiv, seqlen, startpos, endpos, itoken, nanc, iregion;
    char recomfilename[264], tpedfilename[264], regionfilename[264], outfilename[264];
	int hap[2][2], dist, ai, aj, *nminor, *nmajor;
	double dprime, ddenom, r2denom, r2, genDist;
	double p_a, p_b, q_a, q_b; //allele frequencies
	double p_00, p_11, p_01, p_10; //haplotype frequencies
	char *newLine, *token, *running;
	FILE *outf=NULL, *inf=NULL;
	double d, freqi, freqj;
	int *starts, *ends;
	coal_data data;
	int nregions = 0;
	int neut; //BOOLEAN: 0 = False, 1 = True

	if (argc != 5) {
		fprintf(stderr, "Usage: bootstrap_ld_popstats_regions <tpedfilename> <recomfilename> <regionfilename> <outfilename>\n");
		exit(0);
	}

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 
	
    strcpy(tpedfilename, argv[1]);
    strcpy(recomfilename, argv[2]);
    strcpy(regionfilename, argv[3]);
    strcpy(outfilename, argv[4]);

    get_coal_data_tped_vers(&data, tpedfilename, recomfilename);

    //////////////////////
    //// LOAD REGIONS ////
    //////////////////////

	inf = fopen(regionfilename, "r");
	assert(inf != NULL);
	if (inf == NULL) {fprintf(stderr, "Missing regions file: %s\n", regionfilename);}
	//first time through file: get # regions
	while (fgets(newLine, line_size, inf) != NULL) {
		nregions++;
		}
	fclose(inf);
	starts = malloc(nregions * sizeof(int));
	ends = malloc(nregions * sizeof(int));
	//second time through: load all region coordinates into memory
	inf = fopen(regionfilename, "r");
	iregion = 0;
	while (fgets(newLine, line_size, inf) != NULL) {
		for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
			if (itoken == 0) {chromosome = atoi(token);}
			if (itoken == 1) {startpos = atoi(token);}
			if (itoken == 2) {endpos = atoi(token);}
		}
			starts[iregion] = startpos;
			ends[iregion] = endpos;
			iregion++;
	}
	fclose(inf);
	fprintf(stderr, "Loaded info for %d regions\n", nregions);

    /////////////////////////
    //// LOAD GENOTYPES ////
    ////////////////////////

    get_coal_data_tped_vers(&data, tpedfilename, recomfilename);
    nminor = calloc(data.nsnp, sizeof(int));
    nmajor = calloc(data.nsnp, sizeof(int));
    //loop over SNPs
    for (isnp = 0; isnp < data.nsnp; isnp++) {
        if (data.nallele[isnp] != 2) {continue;}
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
        nminor[isnp] = (nderiv > nanc) ? nanc : nderiv;
        nmajor[isnp] = (nderiv > nanc) ? nderiv : nanc;
    } // end loop over snps

    //////////////////////
    //// OUTFILE PREP ///
    /////////////////////

    outf = fopen(outfilename, "w");
    assert(outf != NULL);
    fprintf(stderr, "Writing to file: %s\n", outfilename);

	/**************************************
	 DATA ANALYSIS - LINKAGE DISEQUILIBRIUM
	 **************************************/

	// OUTER LOOP OVER SNPS (isnp)
	for (isnp = 0; isnp < data.nsnp; isnp++) {
		//check if in putative neutral region
		neut = 0;
		for (iregion = 0; iregion < nregions; iregion++){
			if (data.pos[isnp] > starts[iregion] && data.pos[isnp] < ends[iregion]){neut = 1;break;}
		}
		if (neut == 0){continue;}
    	if (nminor[isnp] < min_minor){continue;}

    	// INNER LOOP OVER SNPS (jsnp)
		for (jsnp = isnp+1; jsnp < data.nsnp; jsnp++) {
			//check if in putative neutral region
       		neut = 0;
			for (iregion = 0; iregion < nregions; iregion++){
				if (data.pos[jsnp] > starts[iregion] && data.pos[jsnp] < ends[iregion]){neut = 1;break;}
			}
			if (neut == 0){continue;}
			if (nminor[jsnp] < min_minor) {continue;}
			dist = data.pos[jsnp] - data.pos[isnp];
			genDist = getGenDist(&data, data.pos[isnp], data.pos[jsnp]);
			if (dist > maxdist && genDist > maxgendist) {break;}

			//loop over samples at these two SNPs and count haplotypes
			hap[0][0] = hap[0][1] = hap[1][0] = hap[1][1] = 0;
			for (isamp = 0; isamp < data.nsample; isamp++) {
				ai = data.genotype[isamp][isnp];
				aj = data.genotype[isamp][jsnp];
				hap[ai][aj]++;
			}
			
			p_00 = ((double)hap[0][0]/data.nsample);
			p_01 = ((double)hap[0][1]/data.nsample);
			p_10 = ((double)hap[1][0]/data.nsample);
			p_11 = ((double)hap[1][1]/data.nsample);
			p_a = ((double)(hap[0][0] + hap[0][1])/data.nsample);//locus a, indexed by isnp
			q_a = ((double)(hap[1][0] + hap[1][1])/data.nsample);//locus a, indexed by isnp
			p_b = ((double)(hap[0][0] + hap[1][0])/data.nsample);//locus b, indexed by jsnp
			q_b = ((double)(hap[0][1] + hap[1][1])/data.nsample);//locus b, indexed by isnp
			
			//only continue if both sites polymorphic in this pop
			if (((p_a != 0) && (q_a !=0)) && ((p_b !=0) && (q_b !=0))){
				d = (p_00 * p_11) - (p_01 * p_10);
				r2denom = (p_a * q_a * p_b * q_b);
				r2 = d * d / r2denom;
				
				//d pos
				if (d > 0.){ddenom = ((q_a * p_b) > (p_a * q_b)) ? (p_a * q_b) : (q_a * p_b);}
				//d neg
				else{ddenom = ((-p_a * p_b) > (-q_a*q_b)) ? (-p_a * p_b) : (-q_a*q_b);}
				dprime = d / ddenom;
			
				freqi = (double)nminor[isnp]/data.nsample;
				freqj = (double)nminor[jsnp]/data.nsample;

				fprintf(outf, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n", chromosome, data.pos[isnp], data.pos[jsnp], genDist, freqi, freqj, r2, dprime);
				
			} //end if polymorphic

		} // jsnp loop
	}  // end isnp loop
	free_coal_data(&data);
	fclose(outf);
	return 0;
} // end main