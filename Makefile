## top-level Makefile for cms2.0
## last updated 07.13.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : cms/model/calc_fst_deldaf, cms/composite/combine_scores_poppair

model :
	cd cms/model && $(MAKE)

composite :
	cd cms/composite && $(MAKE)
	
clean :
	rm *.o



