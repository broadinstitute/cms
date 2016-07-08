// NEEDS TO BE VALIDATED/cleaned up. last updated 07.08.16 	vitti@broadinstitute.org

/********************/
/***MULTIPLE POPS****/
/********************/

void get_popComp_data_multiple(popComp_data_multiple* data, int argc, char *argv[]){
  ///collates snps for an arbitrary number of populations compared to a putative selpop,
  ///and returns component scores
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilebase[256], infilesuffix[256], infilename[512];
	int  isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique; //i,maxNsnp
	popComp_data data_sing;
	int *allSnps, *allUniqueSnps;

	//////////////////
	/// INITIALIZE ///
	//////////////////

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

  //////////////////////////
  /// COLLATE LOCI: LOAD ///
  /////////////////////////

	nComparisons = argc; //MUST CHECK THIS IS BEING PASSED CORRECTLY
	totNsnp = 0;
	for (iComp = 1; iComp <= nComparisons; iComp++){
		//sprintf(infilename, infilebase);
		//strcat(infilename, argv[iComp + 6]);
		//strcat(infilename, infilesuffix);
		sprintf(infilename, argv[iComp])
		inf = fopen(infilename, "r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header
		while (fgets(newLine, line_size, inf) != NULL){totNsnp++;}
		fclose(inf);
	} // end iComp
	//fprintf(stderr, "Found %d loci \n", totNsnp);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 1; iComp <= nComparisons; iComp++){
		//sprintf(infilename, infilebase);
		//strcat(infilename, argv[iComp + 6]);
		//strcat(infilename, infilesuffix);
		sprintf(infilename, argv[iComp])
		inf = fopen(infilename,"r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header

		while (fgets(newLine, line_size, inf) != NULL){
			for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 1) {
				 thisPhysPos = atoi(token);
				 allSnps[isnp] =thisPhysPos;
				 break;
				 } // end if
		  } // end for 
      isnp +=1;   
    } // end while
    fclose(inf);
  }//end icomp

  //////////////////////////
  /// COLLATE LOCI: SORT ///
  /////////////////////////

	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	fprintf(stderr, "Sorted SNPs from pos %d to %d\n", allSnps[0], allSnps[totNsnp-1]);
	nunique = 0;
	for (isnp = 0; isnp <= totNsnp-1; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		else{nunique++;}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison...\n", nunique);

  allUniqueSnps  = malloc(nunique * sizeof(int));
  jsnp = 0;
  for (isnp = 0; isnp <= totNsnp-1; isnp++){
    //fprintf(stderr, "%d\t", allSnps[isnp]);
    if (allSnps[isnp] == allSnps[isnp+1]){continue;}
    else{allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
  } // end for isnp

  ///////////////////////
  /// ALLOCATE MEMORY ///
  //////////////////////

	fprintf(stderr, "Allocating memory...\n");
	data->nsnps = nunique;
  data->ncomp =nComparisons;
	data->physpos = malloc(nComparisons * sizeof(int*));
	data->genpos = malloc(nComparisons * sizeof(double*));
	data->daf_selpop = malloc(nComparisons * sizeof(double*));
	data->delDAF = malloc(nComparisons * sizeof(double*));
	data->fst = malloc(nComparisons * sizeof(double*));
	data->xp_normed = malloc(nComparisons * sizeof(double*));
	data->ihs_normed = malloc(nComparisons * sizeof(double*));
	data->delihh_normed = malloc(nComparisons * sizeof(double*));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);
 
	for (iComp = 0; iComp < nComparisons; iComp++){
		data->physpos[iComp] = calloc(nunique, sizeof(int));
		data->genpos[iComp] = calloc(nunique, sizeof(double));
		data->daf_selpop[iComp] = calloc(nunique, sizeof(double));
		data->delDAF[iComp] = calloc(nunique, sizeof(double));
		data->fst[iComp] = calloc(nunique, sizeof(double));
		data->xp_normed[iComp] = calloc(nunique, sizeof(double));		
		data->ihs_normed[iComp] = calloc(nunique, sizeof(double));
		data->delihh_normed[iComp] = calloc(nunique, sizeof(double));		
		assert(data->physpos[iComp] != NULL);
		assert(data->genpos[iComp] != NULL);
		assert(data->daf_selpop[iComp] != NULL);
		assert(data->delDAF[iComp] != NULL);
		assert(data->fst[iComp] != NULL);
		assert(data->xp_normed[iComp] != NULL);
		assert(data->ihs_normed[iComp] != NULL);
		assert(data->delihh_normed[iComp] != NULL);
	} // end for icomp

	/////////////////////////////////////////////
	// LOAD ALL COMPARISONS TO ONE DATA OBJECT //
	/////////////////////////////////////////////

	fprintf(stderr, "Loading all component scores...\n");
	for (iComp = 0; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 6]);
		strcat(infilename, infilesuffix);
		get_popComp_data(&data_sing, infilename);

            //isnp iterates (0, nunique) over allUnique Snps 
    jsnp = 0; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
      //fprintf(stderr, "%d\t%d\t%d\t%d\n", isnp, jsnp, allUniqueSnps[isnp], data_sing.physpos[jsnp]);

      if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
        data->physpos[iComp][isnp] = data_sing.physpos[jsnp];  
        data->genpos[iComp][isnp] = data_sing.genpos[jsnp];   
        data->daf_selpop[iComp][isnp] = data_sing.daf_selpop[jsnp];   
        data->delDAF[iComp][isnp] = data_sing.delDAF[jsnp];  
        data->fst[iComp][isnp] = data_sing.fst[jsnp];   
        data->xp_normed[iComp][isnp] = data_sing.xp_normed[jsnp];               
        data->ihs_normed[iComp][isnp] = data_sing.ihs_normed[jsnp];   
        data->delihh_normed[iComp][isnp] = data_sing.delihh_normed[jsnp];       
        jsnp++; assert(jsnp<=data_sing.nsnps);
      }
      else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; assert(jsnp<=data_sing.nsnps);}
      //else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
    }// end for isnp loop

		free_popComp_data(&data_sing); 
	} // end for icomp
} //end method

void free_popComp_data_multiple(popComp_data_multiple* data){
	int iComp;
	if (data == NULL) {return;}
	for (iComp = 0; iComp < data->ncomp; iComp++){
		free(data->physpos[iComp]);
		free(data->genpos[iComp]);
		free(data->daf_selpop[iComp]);
		free(data->delDAF[iComp]);
		free(data->fst[iComp]);
		free(data->xp_normed[iComp]);
		free(data->ihs_normed[iComp]);
		free(data->delihh_normed[iComp]);
	}
	free(data->physpos);
	free(data->genpos);
	free(data->daf_selpop);
	free(data->delDAF);
	free(data->fst);
	free(data->xp_normed);
	free(data->ihs_normed);
	free(data->delihh_normed);
	data->nsnps = 0;
	data->ncomp = 0;
} //end method

void get_popComp_data_multiple_region(popComp_data_multiple* data, int argc, char *argv[]){
  ///collates snps for an arbitrary number of populations compared to a putative selpop,
  ///and returns component scores
	const int line_size = 15000000; 
	FILE *inf=NULL;
	char *newLine, *token, *running;
	char infilebase[256], infilesuffix[256], infilename[512];
	int isnp, jsnp, itoken, iComp, nComparisons, totNsnp, thisPhysPos, nunique; //maxNsnp, i, 
	int startBp, endBp;
	popComp_data data_sing;
	int *allSnps, *allUniqueSnps;
	int iLine;

	//////////////////
	/// INITIALIZE ///
	//////////////////

	newLine = malloc((line_size+1) * sizeof(char));
	assert(newLine != NULL); 

	data->nsnps = 0;
	data->physpos = NULL;
	data->genpos = NULL;
	data->daf_selpop = NULL;
	data->delDAF = NULL;
	data->fst = NULL;
	data->xp_normed = NULL;
	data->ihs_normed = NULL;
	data->delihh_normed = NULL;

	startBp = atoi(argv[1]);
	endBp = atoi(argv[2]);

  //////////////////////////
  /// COLLATE LOCI: LOAD ///
  /////////////////////////

	nComparisons = argc - 5;
	totNsnp = 0;
	for (iComp = 1; iComp <= nComparisons; iComp++){
		//sprintf(infilename, infilebase);
		//strcat(infilename, argv[iComp + 6]);
		//strcat(infilename, infilesuffix);
		sprintf(infilename, argv[iComp])
		//fprintf(stderr, infilename);
		inf = fopen(infilename, "r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header
		while (fgets(newLine, line_size, inf) != NULL){totNsnp++;}
		fclose(inf);
	} // end iComp
	//fprintf(stderr, "Found %d loci \n", totNsnp);

	//then get array for all of them
	allSnps = malloc(totNsnp * sizeof(int));
	isnp = 0;
	for (iComp = 1; iComp < nComparisons; iComp++){
		sprintf(infilename, infilebase);
		strcat(infilename, argv[iComp + 5]);
		strcat(infilename, infilesuffix);

		inf = fopen(infilename,"r");
		assert(inf != NULL);
		fgets(newLine, line_size, inf); // header

		while (fgets(newLine, line_size, inf) != NULL){
			for (running = newLine, itoken = 0; (token = strsep(&running, "\t")) != NULL; itoken++) {
				if (itoken == 1) {
				 thisPhysPos = atoi(token);
				 allSnps[isnp] =thisPhysPos;
				 break;
				 } // end if
		  } // end for 
      isnp +=1;   
    } // end while
    fclose(inf);
  }//end icomp

  //////////////////////////
  /// COLLATE LOCI: SORT ///
  /////////////////////////

	qsort(allSnps, totNsnp, sizeof(int), intcmp);
	fprintf(stderr, "Sorted SNPs from pos %d to %d\n", allSnps[0], allSnps[totNsnp-1]);
	nunique = 0;
	for (isnp = 0; isnp < totNsnp-1; isnp++){
		//fprintf(stderr, "%d\t", allSnps[isnp]);
		if (allSnps[isnp] == allSnps[isnp+1]){continue;}
		//else{nunique++;}
		else{
			if (allSnps[isnp] >= startBp && allSnps[isnp] <= endBp){nunique++; iLine = isnp;}
		}
	} // end for isnp
	fprintf(stderr, "Found %d SNPs with values for at least one pop comparison within region...\n", nunique);

  allUniqueSnps  = malloc(nunique * sizeof(int));
  jsnp = 0;
  for (isnp = 0; isnp < totNsnp-1; isnp++){
    //fprintf(stderr, "%d\t", allSnps[isnp]);
    if (allSnps[isnp] == allSnps[isnp+1]){continue;}
    else{
    	if (allSnps[isnp] >= startBp && allSnps[isnp] <= endBp)
    		{ //fprintf(stderr, " %d ", allSnps[isnp]);
    		allUniqueSnps[jsnp] = allSnps[isnp]; jsnp++;}
   }
  } // end for isnp

  ///////////////////////
  /// ALLOCATE MEMORY ///
  //////////////////////

	fprintf(stderr, "Allocating memory...\n");
	data->nsnps = nunique;
  data->ncomp =nComparisons;
	data->physpos = malloc(nComparisons * sizeof(int*));
	data->genpos = malloc(nComparisons * sizeof(double*));
	data->daf_selpop = malloc(nComparisons * sizeof(double*));
	data->delDAF = malloc(nComparisons * sizeof(double*));
	data->fst = malloc(nComparisons * sizeof(double*));
	data->xp_normed = malloc(nComparisons * sizeof(double*));
	data->ihs_normed = malloc(nComparisons * sizeof(double*));
	data->delihh_normed = malloc(nComparisons * sizeof(double*));
	assert(data->physpos != NULL);
	assert(data->genpos != NULL);
	assert(data->daf_selpop != NULL);
	assert(data->delDAF != NULL);
	assert(data->fst != NULL);
	assert(data->xp_normed != NULL);
	assert(data->ihs_normed != NULL);
	assert(data->delihh_normed != NULL);
 
	for (iComp = 0; iComp < nComparisons; iComp++){
		data->physpos[iComp] = calloc(nunique, sizeof(int));
		data->genpos[iComp] = calloc(nunique, sizeof(double));
		data->daf_selpop[iComp] = calloc(nunique, sizeof(double));
		data->delDAF[iComp] = calloc(nunique, sizeof(double));
		data->fst[iComp] = calloc(nunique, sizeof(double));
		data->xp_normed[iComp] = calloc(nunique, sizeof(double));		
		data->ihs_normed[iComp] = calloc(nunique, sizeof(double));
		data->delihh_normed[iComp] = calloc(nunique, sizeof(double));		
		assert(data->physpos[iComp] != NULL);
		assert(data->genpos[iComp] != NULL);
		assert(data->daf_selpop[iComp] != NULL);
		assert(data->delDAF[iComp] != NULL);
		assert(data->fst[iComp] != NULL);
		assert(data->xp_normed[iComp] != NULL);
		assert(data->ihs_normed[iComp] != NULL);
		assert(data->delihh_normed[iComp] != NULL);
	} // end for icomp

	/////////////////////////////////////////////
	// LOAD ALL COMPARISONS TO ONE DATA OBJECT //
	/////////////////////////////////////////////

	fprintf(stderr, "Loading all component scores...\n");
	for (iComp = 1; iComp <= nComparisons; iComp++){
		//sprintf(infilename, infilebase);
		//strcat(infilename, argv[iComp + 6]);
		//strcat(infilename, infilesuffix);
		sprintf(infilename, argv[iComp])

		//fprintf(stderr, "CHECK THIS ");
		//fprintf(stderr, infilename);
		//fprintf(stderr, "\n");
		get_popComp_data_region(&data_sing, infilename, startBp, endBp);

		fprintf(stderr, "Region starts with SNP at pos: %d\n", data_sing.physpos[0]);

            //isnp iterates (0, nunique) over allUnique Snps 
    jsnp = 0; // jsnp runs (0, data_sing.nsnp) over data_sing.physpos, smaller range.
		for (isnp = 0; isnp < nunique; isnp++){
      //fprintf(stderr, "%d\t%d\t%d\t%d\n", isnp, jsnp, allUniqueSnps[isnp], data_sing.physpos[jsnp]);

      if (allUniqueSnps[isnp] == data_sing.physpos[jsnp]){ // the snp matches; load all data
        data->physpos[iComp-1][isnp] = data_sing.physpos[jsnp];  
        data->genpos[iComp-1][isnp] = data_sing.genpos[jsnp];   
        data->daf_selpop[iComp-1][isnp] = data_sing.daf_selpop[jsnp];   
        data->delDAF[iComp-1][isnp] = data_sing.delDAF[jsnp];  
        data->fst[iComp-1][isnp] = data_sing.fst[jsnp];   
        data->xp_normed[iComp-1][isnp] = data_sing.xp_normed[jsnp];               
        data->ihs_normed[iComp-1][isnp] = data_sing.ihs_normed[jsnp];   
        data->delihh_normed[iComp-1][isnp] = data_sing.delihh_normed[jsnp];       
        jsnp++; assert(jsnp<=data_sing.nsnps);
      }
      else if (allUniqueSnps[isnp] > data_sing.physpos[jsnp]){jsnp++; assert(jsnp<=data_sing.nsnps);}
      //else if (allUniqueSnps[isnp] < data_sing.physpos[jsnp]){pass;}
    }// end for isnp loop

		free_popComp_data(&data_sing); 
	} // end for icomp
} //end method

void free_popComp_data_multiple_region(popComp_data_multiple* data){
	int iComp;
	if (data == NULL) {return;}
	for (iComp = 0; iComp < data->ncomp; iComp++){
		free(data->physpos[iComp]);
		free(data->genpos[iComp]);
		free(data->daf_selpop[iComp]);
		free(data->delDAF[iComp]);
		free(data->fst[iComp]);
		free(data->xp_normed[iComp]);
		free(data->ihs_normed[iComp]);
		free(data->delihh_normed[iComp]);
	}
	free(data->physpos);
	free(data->genpos);
	free(data->daf_selpop);
	free(data->delDAF);
	free(data->fst);
	free(data->xp_normed);
	free(data->ihs_normed);
	free(data->delihh_normed);
	data->nsnps = 0;
	data->ncomp = 0;
} //end method