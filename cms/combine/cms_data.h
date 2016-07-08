// datastructures and function declarations for handling cms component(+composite) score datastructures
// last updated: 07.08.16   vitti@broadinstitute.org

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
    int *lastcol; //not sure what information this field contains, selscan documentation is sparse. 0/1
    double *freq1;
    double *ihh0;
    double *ihh1;
    double *ihs_unnormed;
    double *ihs_normed; 
} ihs_data;
void get_ihs_data(ihs_data* data, char filename[]);
void free_ihs_data(ihs_data* data);

typedef struct fst_deldaf_data {
    int nsnps;
    int *pos; 
    double *fst;
    double *deldaf;
} fst_deldaf_data;
void get_fst_deldaf_data(fst_deldaf_data* data, char filename[]);
void free_fst_deldaf_data(fst_deldaf_data* data);

/*************************/
/***SCORE DISTRIBUTIONS***/
/*************************/

typedef struct likes_data{
    int nbins;
    double *start_bin;
    double *end_bin;
    double *probs;
} likes_data;
void get_likes_data(likes_data* data, char filename[]);
void free_likes_data(likes_data* data);

typedef struct popComp_data{
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
} popComp_data;

void get_popComp_data(popComp_data* data, char filename[]); //_multiple?
//void get_popComp_data_region(popComp_data* data, char filename[], int startBp, int endBp);
void free_popComp_data(popComp_data* data);

void get_popComp_data_region(popComp_data* data, char infilename[], int startBp, int endBp);

