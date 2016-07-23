// data structure to facilitate parsing of pop ihh files
// /idi/sabeti-scratch/jvitti/cms_venv
// 10.12.15

typedef struct ihh_data {
    int nsnps;
    int *pos; 
    double *genpos;
    double *freq;
    double *ihh; 
} ihh_data;

void get_pop_ihh_data(ihh_data* data, char filename[]);
void free_pop_ihh_data(ihh_data* data);
