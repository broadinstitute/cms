//last updated 11.15.16 	vitti@broadinstitute.org

/**********************/
/***POP COMPARISONS****/
/**********************/

//VARIABLE ARGS FUNCTION
typedef struct popComp_data_multiple{
    int nsnps;
    int ncomp; // how many non-sel pops
    int **physpos;
    double **genpos;
    double **daf_selpop;
    double **delDAF;
    double **fst;
    double **xp_normed;
    double **ihs_normed;
    double **delihh_normed;
} popComp_data_multiple;
void get_popComp_data_multiple(popComp_data_multiple* data, int argc, char *argv[]);
void get_popComp_data_multiple_region(popComp_data_multiple* data, int argc, char *argv[]);
void free_popComp_data_multiple(popComp_data_multiple* data);

float getMinBf(likes_data* data_hit, likes_data* data_miss);
float getProb(likes_data* data, double value);
float compareXp(popComp_data_multiple* data, int isnp);
float compareFst(popComp_data_multiple* data, int isnp);
float comparedelDaf(popComp_data_multiple* data, int isnp);
