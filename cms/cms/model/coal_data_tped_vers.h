// {{POP GEN DATA --> DEM MODEL}: CALC POP SUMMARY STATS}
// {{POP GEN DATA --> COMPONENT SCORES}}
// last updated: 06.20.17   vitti@broadinstitute.org

/************************/
/**DEFINE DATASTRUCTURE**/
/************************/
typedef struct coal_data {
    int nsample;
    char **samp_id;
    int nsnp;
    char **snp_id;
    int *pos;
    int *chrom;
    int **genotype;  // haploid
    char **snp_base[4];
    int *anc_base;
    int *nallele;
    double *genPos; //for determining gen dist
    double *genloc;
    int *physPos; //parallel w above
    int nRecom; //# of lines in recombination file, as nsample or nsnp
} coal_data;

/************************/
/***DECLARE FUNCTIONS****/
/************************/
/*helper methods*/
double getGenDist(coal_data* data, int pos_i, int pos_j);
int getIndexOfItem(int *values, int numVals, int itemToFind);
int getLowerIndexOfItem(int *values, int numVals, int itemToFind);
int getUpperIndexOfItem(int *values, int numVals, int itemToFind);
/*manipulating datastrucutre*/
void get_coal_data_tped_vers(coal_data* data, char tpedfilename[], char recomfilename[]);
void get_coal_data_tped_vers_gz(coal_data* data, char tpedfilename[], char recomfilename[]);
void free_coal_data(coal_data* data);