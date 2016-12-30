## leaf Makefile for cms2.0/combine
## last updated 12.30.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 


##################
## DEFINE RULES ##
##################

all : combine_scores  write_xpehh_from_ihh  #install

combine_scores.o : combine_scores.c
	$(CC) -c combine_scores.c

combine_scores : combine_scores.o cms_data.c cms_data.h
	$(CC) $(CCFLAG) -o combine_scores combine_scores.o cms_data.c

write_xpehh_from_ihh.o : 
	$(CC) -c write_xpehh_from_ihh.c

write_xpehh_from_ihh : write_xpehh_from_ihh.o cms_data.c pop_ihh_data.h  pop_ihh_data.c 
	$(CC) $(CCFLAG) -o write_xpehh_from_ihh write_xpehh_from_ihh.o pop_ihh_data.c 



clean :
	rm *.o && rm combine_scores && rm write_xpehh_from_ihh &&  rm -R *.dSYM/
install:
	cp combine_scores /usr/local/ && cp write_xpehh_from_ihh /usr/local

