## Makefile for CMS 2.0
## last updated 07.14.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : combine/combine_scores_poppair, calc_fst_deldaf
#combine_scores_poppair combine_scores_multiplepops bootstrap_freq_popstats_regions bootstrap_fst_popstats_regions bootstrap_ld_popstats_regions

calc_fst_deldaf : model/coal_data_tped_vers.c model/coal_data_tped_vers.h calc_fst_deldaf.o 
	$(CC) $(CCFLAG) -o calc_fst_deldaf calc_fst_deldaf.o model/coal_data_tped_vers.c

calc_fst_deldaf.o : calc_fst_deldaf.c
	$(CC) $(CCFLAG) -c calc_fst_deldaf.c 

combine/combine_scores_poppair : combine_scores_poppair.o combine/cms_data.c combine/cms_data.h
	$(CC) $(CCFLAG) -o combine/combine_scores_poppair combine_scores_poppair.o combine/cms_data.c

combine_scores_poppair.o : combine/combine_scores_poppair.c
	$(CC) $(CCFLAG) -c combine/combine_scores_poppair.c

#combine_scores_multiplepops : 

#combine_scores_multiplepops.o : 

#bootstrap_freq_popstats_regions :

#bootstrap_freq_popstats_regions.o : 

#bootstrap_fst_popstats_regions : 

#bootstrap_fst_popstats_regions.o : 

#bootstrap_ld_popstats_regions : 

#bootstrap_ld_popstats_regions.o:


clean :
	rm *.o



