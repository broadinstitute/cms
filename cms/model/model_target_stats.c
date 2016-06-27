// experimental; trying to wrap current bootstrap framework, patterning after selscan-main.cpp 
// last updated: 06.21.16
// possibly preferable to wrap in python
#include <string>
#include "param_t.h"

using namespace std;

const string VERSION = "test";
const string PREAMBLE = "\ncms_modeller v" + VERSION + " -- a program for exploratory fitting of demographic models to population genetic data as part of CMS 2.0 package.\n\
Source code and binaries can be found at https://github.com/broadinstitute/cms>.\n\
\n\
To calculate bootstrap values of per-pop target statistics:\n\
\n\
./cms_modeller --target --tped <pop inputfile> --map <mapfile> --out <outfile>\n\
\n\
To calculate XP-EHH:\n\
\n\
./cms_modeller --xpehh --hap <pop1 haps> --ref <pop2 haps> --map <mapfile> --out <outfile>\n";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = 1;
const string HELP_THREAD = "The number of threads to spawn during the calculation. Partitions loci across threads.";

const string ARG_FILENAME_POP1_TPED = "--tped";
const string DEFAULT_FILENAME_POP1_TPED = "__hapfile1";
const string HELP_FILENAME_POP1_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1";


