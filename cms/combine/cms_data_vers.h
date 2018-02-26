
int intcmp(const void *v1, const void *v2);
int count_unique_from_sorted(int arr[], int len);

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
void get_xpehh_data(xpehh_data* data, char filename[], int minPos, int maxPos);
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
void get_delihh_data(delihh_data* data, char filename[], int minPos, int maxPos);
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
void get_ihs_data(ihs_data* data, char filename[], int minPos, int maxPos);
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
void get_nsl_data(nsl_data* data, char filename[], int minPos, int maxPos);
void free_nsl_data(nsl_data* data);

typedef struct freqs_data {
    int nsnps;
    int *pos; 
    double *popdaf;
    double *fst;
    double *deldaf;
    double *genpos;
} freqs_data;
void get_freqs_data(freqs_data* data, char filename[], int minPos, int maxPos);
void free_freqs_data(freqs_data* data);

typedef struct H12_data {
    int nsnps;
    int *pos; 
    double *H12_value;
    double *H2H1_value;
} H12_data;
void get_H12_data(H12_data* data, char filename[], int minPos, int maxPos);
void free_H12_data(H12_data* data);

typedef struct iSAFE_data {
    int nsnps;
    int *pos; 
    double *iSAFE_value;
} iSAFE_data;
void get_iSAFE_data(iSAFE_data* data, char filename[], int minPos, int maxPos);
void free_iSAFE_data(iSAFE_data* data);

/*************************/
/***TWO-POP COMPARISON***/
/*************************/
int get_num_anyData(int minPos, int maxPos, char ihs_filename[], char delihh_filename[], char nsl_filename[], char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[]); 
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
    double *H12;
    double *H2H1;
    double *iSAFE;
} popPair_data;
void get_popPair_anyData(int minPos, int maxPos, popPair_data* data, char ihs_filename[], char delihh_filename[], char nsl_filename[],char H12_filename[], char iSAFE_filename[], char xpehh_filename[], char freqs_filename[], int preCountedSNPs); 
void free_popPair_data(popPair_data* data);

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
    double **H12;
    double **H2H1;
    double **iSAFE;    
} popComp_data_multiple;
void get_popComp_anyData(int minPos, int maxPos, popComp_data_multiple* data, int nComparisons, int argc, char *argv[]);
void free_popComp_data_multiple(popComp_data_multiple* data);
float compareXp(popComp_data_multiple* data, int isnp);
float compareFst(popComp_data_multiple* data, int isnp);
float comparedelDaf(popComp_data_multiple* data, int isnp);
float compareFst_PBS(popComp_data_multiple* data, int isnp);
float comparedelDaf_outgroup_ave(popComp_data_multiple* data, int isnp);
float get_T(double fst);
float get_PBS(double in_t_1, double in_t_2, double out_t);
float get_outgroups_fst(popComp_data_multiple* data, int isnp, int iComp, int jComp);