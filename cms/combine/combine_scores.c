// 12.28.16: new top-level program for compositing 		vitti@broadinstitute.org

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
	popComp_data_multiple data;
	FILE *outf=NULL;
	char outfilename[256]; 

	if (argc != 2) {
		fprintf(stderr, "Usage: ./combine_scores <cms_run_paramfile>\n");
		exit(0);
	}

	//<cms_run_paramfile> has:
	// selpop and altpops and dem model and other options
	//	some number of (pop-pairs), with pointers to score files
	//maybe each popComp data passes a file that is just the list of filenames 

	//////////////////
	// LOAD SCORES ///
	//////////////////
	fprintf(stderr, "Preparing to load component scores...\n");
	get_popComp_data_multiple(&data, argc, argv); 

	/////////////////////////////
	// LOAD SCORE LIKELIHOODS ///
	/////////////////////////////
	fprintf(stderr, "Preparing to load score likelihoods...\n");
	//get_likes_data_multiple(&data, argc, argv); 

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: %s\n", outfilename);
	for (isnp = 0; isnp < data.nsnps; isnp++){
		//////////////////////////////////
		//HANDLE POPULATION COMPARISONS //
		//////////////////////////////////


		//compLikeRatio = delihh_bf * nsl_bf  * fst_bf * deldaf_bf * xpehh_bf; //* ihs_bf
		//fprintf(outf, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", data.physpos[iComp][isnp], data.genpos[iComp][isnp], thisihs, thisihh, thisnsl, thisxpehh, thisfst, thisdelDaf, compLikeRatio);
	} // end isnp
	fclose(outf);
	free_popComp_data_multiple(&data);
	return 0;

} // end main