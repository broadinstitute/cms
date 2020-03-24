## leaf Makefile for cms2.0/combine
## last updated 12.01.17 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall -lz 

##################
## DEFINE RULES ##
##################

all : combine_scores_local_zipped  combine_scores_gw_zipped write_xpehh_from_ihh  #install

# combine_scores_local combine_scores_gw

#collate_scores.o : collate_scores.c
#	$(CC) -c collate_scores.c

#collate_scores : collate_scores.o cms_data.c cms_data.h
#	$(CC) $(CCFLAG) -o collate_scores collate_scores.o cms_data.c

#collate_scores_zipped : collate_scores.o cms_data_zipped.c cms_data.h
#	$(CC) $(CCFLAG) -o collate_scores_zipped collate_scores.o cms_data_zipped.c


combine_scores_local.o : combine_scores_local.c
	$(CC) -c combine_scores_local.c

combine_scores_local : combine_scores_local.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_local combine_scores_local.o cms_data.c

combine_scores_gw.o : combine_scores_gw.c
	$(CC) -c combine_scores_gw.c

combine_scores_gw : combine_scores_gw.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_gw combine_scores_gw.o cms_data.c

combine_scores_local_zipped : combine_scores_local.o cms_data_zipped.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_local_zipped combine_scores_local.o cms_data_zipped.c

combine_scores_gw_zipped : combine_scores_gw.o cms_data_zipped.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores_gw_zipped combine_scores_gw.o cms_data_zipped.c

write_xpehh_from_ihh.o : 
	$(CC) -c write_xpehh_from_ihh.c

write_xpehh_from_ihh : write_xpehh_from_ihh.o cms_data.c pop_ihh_data.h  pop_ihh_data.c 
	$(CC) $(CCFLAG) -o write_xpehh_from_ihh write_xpehh_from_ihh.o pop_ihh_data.c 


clean :
	rm *.o && rm combine_scores_local_zipped && rm combine_scores_gw_zipped
	# && rm write_xpehh_from_ihh && rm combine_scores_local_zipped && rm combine_scores_gw_zipped
	#&& rm -R *.dSYM/
#install:
#	cp combine_scores /usr/local/ && cp write_xpehh_from_ihh /usr/local

