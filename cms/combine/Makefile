## leaf Makefile for cms2.0/combine
## last updated 06.14.17 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : combine_scores_local combine_scores_local_amend combine_scores_gw write_xpehh_from_ihh  #install

combine_scores_local.o : combine_scores_local.c
	$(CC) -c combine_scores_local.c

combine_scores_local : combine_scores_local.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_local combine_scores_local.o cms_data.c

combine_scores_local_amend.o : combine_scores_local_amend.c
	$(CC) -c combine_scores_local_amend.c

combine_scores_local_amend : combine_scores_local_amend.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_local_amend combine_scores_local_amend.o cms_data.c

combine_scores_gw.o : combine_scores_gw.c
	$(CC) -c combine_scores_gw.c

combine_scores_gw : combine_scores_gw.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_gw combine_scores_gw.o cms_data.c

write_xpehh_from_ihh.o : 
	$(CC) -c write_xpehh_from_ihh.c

write_xpehh_from_ihh : write_xpehh_from_ihh.o cms_data.c pop_ihh_data.h  pop_ihh_data.c 
	$(CC) $(CCFLAG) -o write_xpehh_from_ihh write_xpehh_from_ihh.o pop_ihh_data.c 


clean :
	rm *.o && rm combine_scores_local && rm combine_scores_local_amend && rm combine_scores_gw && rm write_xpehh_from_ihh #&& rm -R *.dSYM/
#install:
#	cp combine_scores /usr/local/ && cp write_xpehh_from_ihh /usr/local

