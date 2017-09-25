// datastructures and function declarations for handling cms component(+composite) score datastructures
// last updated: 09.17.2017 	vitti@broadinstitute.org

int intcmp(const void *v1, const void *v2);
int count_unique(int arr[], int len);

/**********************/
/***COMPONENT SCORES***/
/**********************/
typedef struct xpehh_data {
    int nsnps;
    int *pos; 
    double *genpos;
    double *freq1;
    double *ihh1;
    double *freq2;
    double *ihh2;
    double *xpehh_unnormed;
    double *xpehh_normed;
    int *lastcol;
} xpehh_data; 
void get_xpehh_data(xpehh_data* data, char filename[]);
void free_xpehh_data(xpehh_data* data);

typedef struct delihh_data {
    int nsnps;
    int *pos; 
    double *freq1;
    double *ihs_unnormed;
    double *delihh_unnormed;
    double *delihh_normed;
    int *lastcol;    
} delihh_data;
void get_delihh_data(delihh_data* data, char filename[]);
void free_delihh_data(delihh_data* data);

typedef struct ihs_data {
    int nsnps;
    int *pos; 
    int *lastcol; 
    double *freq1;
    double *ihh0;
    double *ihh1;
    double *ihs_unnormed;
    double *ihs_normed; 
} ihs_data;
void get_ihs_data(ihs_data* data, char filename[]);
void free_ihs_data(ihs_data* data);

typedef struct nsl_data {
    int nsnps;
    int *pos; 
    double *freq1;
    double *sl0;
    double *sl1;
    double *nsl_unnormed;
    double *nsl_normed; 
} nsl_data;
void get_nsl_data(nsl_data* data, char filename[]);
void free_nsl_data(nsl_data* data);

typedef struct freqs_data {
    int nsnps;
    int *pos; 
    double *popdaf;
    double *fst;
    double *deldaf;
    double *genpos;
} freqs_data;
void get_freqs_data(freqs_data* data, char filename[]);
void free_freqs_data(freqs_data* data);

/*************************/
/***SCORE LIKELIHOODS***/
/*************************/
typedef struct likes_data{
    int nbins;
    double *start_bin;
    double *end_bin;
    double *probs;
} likes_data;
void get_likes_data(likes_data* data, char filename[]);
void free_likes_data(likes_data* data);
typedef struct likes_data_multiple{
    int nbins;
    double *start_bin;
    double *end_bin;
    double *miss_probs;      
    double **hit_probs; 
} likes_data_multiple;
void get_likes_data_multiple(likes_data_multiple* data, char filename[]);
void free_likes_data_multiple(likes_data_multiple* data);
float getHitProb(likes_data_multiple* data, int likesIndex, double value);
float getMissProb(likes_data_multiple* data, double value);
float getMaxBf(likes_data_multiple* data, int likesIndex);
float getMinBf(likes_data_multiple* data,  int likesIndex);
float getMaxProb(likes_data_multiple* data, int likesIndex, double prior);
float getMinProb(likes_data_multiple* data,  int likesIndex, double prior);

/*************************/
/***TWO-POP COMPARISON***/
/*************************/
int get_num_completeData(char ihs_filename[], char delihh_filename[], char nsl_filename[], char xpehh_filename[], char freqs_filename[]); 
int get_num_anyData(char ihs_filename[], char delihh_filename[], char nsl_filename[], char xpehh_filename[], char freqs_filename[]); 
typedef struct popPair_data{
    int nsnps;
    char **locus_id;
    int *physpos;
    double *genpos;
    double *daf_selpop;
    double *delDAF;
    double *fst;
    double *xp_normed;
    double *ihs_normed;
    double *delihh_normed;
    double *nsl_normed;
} popPair_data;
void get_popPair_completeData(popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[], char xpehh_filename[], char freqs_filename[]); 
void get_popPair_anyData(popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[], char xpehh_filename[], char freqs_filename[]); 
void free_popPair_data(popPair_data* data);
void get_popPair_data_region(popPair_data* data, char infilename[], int startBp, int endBp);

/**************************/
/***MULTIPOP COMPARISON****/
/**************************/
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
    double **nsl_normed;
} popComp_data_multiple;
void get_popComp_anyData(popComp_data_multiple* data, int nComparisons, int argc, char *argv[]);
void get_popComp_completeData(popComp_data_multiple* data, int nComparisons, int argc, char *argv[]);
void free_popComp_data_multiple(popComp_data_multiple* data);
float compareXp(popComp_data_multiple* data, int isnp);
float compareFst(popComp_data_multiple* data, int isnp);
float comparedelDaf(popComp_data_multiple* data, int isnp);
float compareFst_PBS(popComp_data_multiple* data, int isnp);
float comparedelDaf_outgroup_ave(popComp_data_multiple* data, int isnp);
float get_T(double fst);
float get_PBS(double in_t_1, double in_t_2, double out_t);
float get_outgroups_fst(popComp_data_multiple* data, int isnp, int iComp, int jComp);