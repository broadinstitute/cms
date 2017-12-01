// 	stripped-down CMS functionality: don't calculate composite, just collate scores by SNP
//	last updated 10.19.2017 	vitti@broadinstitute.org

//gcc -c collate_scores.c
//gcc -O0 -ggdb3 -lm -lz -Wall -o collate_scores_fromzipped collate_scores.o cms_data_zipped.c
//gcc -O0 -ggdb3 -lm -Wall -o collate_scores collate_scores.o cms_data.c

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"

/**********/
/***MAIN***/
/**********/
int main(int argc, char **argv) {
	const int line_size = 15000000; 
	popComp_data_multiple score_data;
	FILE *inf=NULL, *outf=NULL, *outf2=NULL;
	char *token, *running;
	char cms_param_filename[528], paramline[528], outfilename[256];
	double thisihs, thisihh, thisnsl; // per-pop
	double thisfst, thisxpehh, thisdelDaf, thisdaf;
	double thisH12, thisH2H1, thisiSAFE;
	int isnp, iComp, itoken, thisPos, likesFreqIndex, nComparisons, maxPos, minPos;
	int istart, iend;
	//int ibin;  //for debug

	if (argc <= 2) {
		fprintf(stderr, "Usage: ./collate_scores <savefilename> <input_pair_file1> ...\n");
		exit(0);
	}
	nComparisons = argc - 2;
	
	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "\nPreparing to load component scores...\n");
	get_popComp_anyData(&score_data, nComparisons, argc, argv);  //BUILD IN H12 SUPPORT 
	fprintf(stderr, "\tloaded data object with %d snps and %d population comparisons.\n", score_data.nsnps, score_data.ncomp);
	
	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	strcpy(outfilename, argv[1]);
	fprintf(stderr, "Preparing to write to: %s\n", outfilename);
	fprintf(stderr, "debug: score_data.nsnps = %d\n", score_data.nsnps);
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(outf, "physPos\tpopDAF\tnormed_iHS\tnormed_deliHH\tnormed_nsl\tH12\tH2H1\tiSAFE\tnormed_xp-ehh\tfst\tdelDAF\n");
	for (isnp = 0; isnp < score_data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////
		iComp = 0; 
		for (iComp = 0; iComp < score_data.ncomp; iComp++){
			if (score_data.physpos[iComp][isnp] != 0){break;}
		} //advance to the first comparison for which we have any data
		if (iComp >= score_data.ncomp){iComp = 0;} //catch SNPs at position 0
		thisihs = score_data.ihs_normed[iComp][isnp];
		thisihh = score_data.delihh_normed[iComp][isnp];
		thisnsl = score_data.nsl_normed[iComp][isnp];
		thisH12 = score_data.H12[iComp][isnp];
		thisH2H1 = score_data.H2H1[iComp][isnp];
		thisiSAFE = score_data.iSAFE[iComp][isnp];
		thisxpehh = compareXp(&score_data, isnp);
		thisfst = compareFst_PBS(&score_data, isnp);	
		thisdelDaf = comparedelDaf_outgroup_ave(&score_data, isnp);	
		thisdaf = score_data.daf_selpop[iComp][isnp];

		fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", score_data.physpos[iComp][isnp], thisdaf, thisihs, thisihh, thisnsl, thisH12, thisH2H1, thisiSAFE, thisxpehh, thisfst, thisdelDaf);
		
	} // end isnp
	fclose(outf);
	fprintf(stderr, "\nWrote to %s\n", outfilename);
	free_popComp_data_multiple(&score_data);
	return 0;
} // end main