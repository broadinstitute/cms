## leaf Makefile for cms2.0/combine
## last updated 07.13.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : combine_scores_poppair, combine_scores_multiplepops 

combine_scores_poppair : combine_scores_poppair.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_poppair combine_scores_poppair.o cms_data.c

combine_scores_poppair.o : combine_scores_poppair.c
	$(CC) $(CCFLAG) -c combine_scores_poppair.c

combine_scores_multiplepops : combine_scores_multiplepops.o cms_data.c cms_data.h  pop_comparison.c pop_comparison.h
	$(CC) $(CCFLAG) -o combine_cms_scores_multiplepops combine_scores_multiplepops.o cms_data.c pop_comparison.c

combine_scores_multiplepops.o : 
	$(CC) $(CCFLAG) -c combine_scores_multiplepops.c
	
clean :
	rm *.o



