// PHASING OUT 12.29.16 (->combine_scores.c)

// for a given population pair, pulls and collates all component score statistics. 
// last updated: 11.23.16   vitti@broadinstitute.org

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"

int main(int argc, char **argv) {
	delihh_data delihh1;
	ihs_data ihs1;
	nsl_data nsl1;
	xpehh_data xp;
	fst_deldaf_data fst_deldaf;
	FILE *outf=NULL;
	int xp_rev, deldaf_rev; // {}_rev are Booleans 0T 1F to track putative selpop. 
	char ihs1filename[528], delihh1filename[528];
	char nsl1filename[528];
	char xpfilename[528], fst_deldaffilename[528], outfilename[528];
	int ihs1_index, delihh1_index, xp_index, fst_deldaf_index, nsl1_index;
	int ihs1pos, delihh1pos, xppos, fst_deldafpos, nsl1pos;
	float thisXp, thisDeldaf; //in case we need to invert to maintain consistent putative selpop
	double minimum;

	if (argc != 9) {
		fprintf(stderr, "Usage: ./combine_cms_scores_poppairs <ihs1filename> <nsl1filename> <delihh1filename> <xpfilename> <xp reversed? 0T 1F> <fst_deldaffilename> <deldaf reversed? 0T 1F> <outfilename>\n");
		exit(0);
	}

	sprintf(ihs1filename, "%s", argv[1]);
	sprintf(nsl1filename, "%s", argv[2]);
	sprintf(delihh1filename, "%s", argv[3]);	
	sprintf(xpfilename, "%s", argv[4]);
	xp_rev = atoi(argv[5]);
	sprintf(fst_deldaffilename, "%s", argv[6]);
	deldaf_rev = atoi(argv[7]);

	fprintf(stderr, "\n"); 
	fprintf(stderr, "loading data from: %s\n", ihs1filename);
	get_ihs_data(&ihs1, ihs1filename);
	fprintf(stderr, "\t nsnps: %d\n", ihs1.nsnps);

	fprintf(stderr, "loading data from: %s\n", nsl1filename);
	get_nsl_data(&nsl1, nsl1filename);
	fprintf(stderr, "\t nsnps: %d\n", nsl1.nsnps);

	fprintf(stderr, "loading data from: %s\n", delihh1filename);
	get_delihh_data(&delihh1, delihh1filename);
	fprintf(stderr, "\t nsnps: %d\n", delihh1.nsnps);

	fprintf(stderr, "loading data from: %s\n", xpfilename);
	get_xpehh_data(&xp, xpfilename);
	fprintf(stderr, "\t nsnps: %d\n", xp.nsnps);

	fprintf(stderr, "loading data from: %s\n", fst_deldaffilename);
	get_fst_deldaf_data(&fst_deldaf, fst_deldaffilename);
	fprintf(stderr, "\t nsnps: %d\n", fst_deldaf.nsnps);

	sprintf(outfilename, "%s", argv[8]);

	////////////////////////
	// ITERATE OVER SNPS ///
	////////////////////////
	ihs1_index=0;
	nsl1_index=0;
	delihh1_index=0;
	xp_index=0;
	fst_deldaf_index=0;

	outf = fopen(outfilename, "w");
	assert(outf != NULL);
	fprintf(stderr, "writing to: %s\n", outfilename);

	fprintf(outf, "locus\tphyspos\tgenpos\tDAF_selpop\tdelDAF\tfst\txp_normed\tihs1_normed\tnsl1_normed\tdelihh1_normed\n");//header
	while (ihs1_index < ihs1.nsnps && delihh1_index < delihh1.nsnps && xp_index < xp.nsnps && fst_deldaf_index < fst_deldaf.nsnps)
	{
		ihs1pos = ihs1.pos[ihs1_index];
		nsl1pos = nsl1.pos[nsl1_index];
		delihh1pos = delihh1.pos[delihh1_index];
		xppos = xp.pos[xp_index];
		fst_deldafpos = fst_deldaf.pos[fst_deldaf_index];

		if (ihs1pos == delihh1pos && delihh1pos == xppos && xppos == fst_deldafpos && nsl1pos == fst_deldafpos){
			thisXp=xp.xpehh_normed[xp_index];
			thisDeldaf=fst_deldaf.deldaf[fst_deldaf_index];
			if(xp_rev == 0){thisXp*=-1;}
			if(deldaf_rev == 0){thisDeldaf*=-1;}

			//fprintf(outf, "chr%s_", argv[1]);
			fprintf(outf, "%d\t%d\t", ihs1pos, ihs1pos);
			fprintf(outf, "%f\t", xp.genpos[xp_index]);
			fprintf(outf, "%f\t", 1 - ihs1.freq1[ihs1_index]);
			fprintf(outf, "%f\t", thisDeldaf);
			fprintf(outf, "%f\t", fst_deldaf.fst[fst_deldaf_index]);
			fprintf(outf, "%f\t", thisXp);			
			fprintf(outf, "%f\t", ihs1.ihs_normed[ihs1_index]);
			fprintf(outf, "%f\t", nsl1.nsl_normed[nsl1_index]);
			fprintf(outf, "%f\n", delihh1.delihh_normed[delihh1_index]);
		}

		//If not at the same point, find out which position is lowest and advance its pointer.
		minimum = 1000000000000;
		if (ihs1pos < minimum){minimum = ihs1pos;}
		if (nsl1pos < minimum){minimum = nsl1pos;}
		if (delihh1pos < minimum){minimum = delihh1pos;}
		if (xppos < minimum){minimum = xppos;}		
		if (fst_deldafpos < minimum){minimum = fst_deldafpos;}

		if (ihs1pos == minimum){ihs1_index++;}
		if (nsl1pos == minimum){nsl1_index++;}
		if (delihh1pos == minimum){delihh1_index++;}
		if (xppos == minimum){xp_index++;}		
		if (fst_deldafpos == minimum){fst_deldaf_index++;}
			
	} // end while loop

	free_ihs_data(&ihs1);
	free_nsl_data(&nsl1);
	free_delihh_data(&delihh1);
	free_xpehh_data(&xp);
	free_fst_deldaf_data(&fst_deldaf);

	return 0;
} // end main