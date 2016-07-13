## leaf Makefile for cms2.0/model
## last updated 07.13.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : calc_fst_deldaf

calc_fst_deldaf.o : calc_fst_deldaf.c
	$(CC) $(CCFLAG) -c calc_fst_deldaf.c

calc_fst_deldaf : calc_fst_deldaf.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o calc_fst_deldaf calc_fst_deldaf.o coal_data_tped_vers.c
	
clean :
	rm *.o



