## leaf Makefile for cms2.0/combine
## last updated 11.04.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 


##################
## DEFINE RULES ##
##################

all : combine_scores_poppair combine_scores_multiplepops combine_scores_multiplepops_region write_xpehh_from_ihh combine_scores_multiplepops_decompose #install

combine_scores_poppair.o : combine_scores_poppair.c
	$(CC) -c combine_scores_poppair.c

combine_scores_poppair : combine_scores_poppair.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_poppair combine_scores_poppair.o cms_data.c

combine_scores_multiplepops_decompose.o : 
	$(CC) -c combine_scores_multiplepops_decompose.c

combine_scores_multiplepops_decompose : combine_scores_multiplepops_decompose.o cms_data.c cms_data.h  pop_comparison.c pop_comparison.h
	$(CC) $(CCFLAG) -o combine_scores_multiplepops_decompose combine_scores_multiplepops_decompose.o cms_data.c pop_comparison.c


combine_scores_multiplepops.o : 
	$(CC) -c combine_scores_multiplepops.c

combine_scores_multiplepops : combine_scores_multiplepops.o cms_data.c cms_data.h  pop_comparison.c pop_comparison.h
	$(CC) $(CCFLAG) -o combine_scores_multiplepops combine_scores_multiplepops.o cms_data.c pop_comparison.c


combine_scores_multiplepops_region.o : 
	$(CC) -c combine_scores_multiplepops_region.c

combine_scores_multiplepops_region : combine_scores_multiplepops_region.o cms_data.c cms_data.h  pop_comparison.c pop_comparison.h
	$(CC) $(CCFLAG) -o combine_scores_multiplepops_region combine_scores_multiplepops_region.o cms_data.c pop_comparison.c


write_xpehh_from_ihh.o : 
	$(CC) -c write_xpehh_from_ihh.c

write_xpehh_from_ihh : write_xpehh_from_ihh.o cms_data.c pop_ihh_data.h  pop_ihh_data.c 
	$(CC) $(CCFLAG) -o write_xpehh_from_ihh write_xpehh_from_ihh.o pop_ihh_data.c 

#$(phony install): all
#	install -m 0755 combine_scores_poppair /usr/local/bin/


clean :
	rm *.o && rm combine_scores_poppair && rm combine_scores_multiplepops && rm combine_scores_multiplepops_region && rm write_xpehh_from_ihh && rm combine_scores_multiplepops_decompose &&  rm -R *.dSYM/
install:
	cp combine_scores_poppair /usr/local/ && cp combine_scores_multiplepops /usr/local && cp combine_scores_multiplepops_region /usr/local && cp write_xpehh_from_ihh /usr/local

