## leaf Makefile for cms2.0/model
## last updated 06.13.17 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : freqs_stats calc_fst_deldaf  bootstrap_freq_popstats_regions bootstrap_ld_popstats_regions bootstrap_fst_popstats_regions #calc_fst_deldaf_fromzipped


freqs_stats.o : freqs_stats.c
	$(CC) $(CCFLAG) -c freqs_stats.c

freqs_stats : freqs_stats.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o freqs_stats freqs_stats.o coal_data_tped_vers.c




calc_fst_deldaf.o : calc_fst_deldaf.c
	$(CC) $(CCFLAG) -c calc_fst_deldaf.c

calc_fst_deldaf : calc_fst_deldaf.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o calc_fst_deldaf calc_fst_deldaf.o coal_data_tped_vers.c


#calc_fst_deldaf_fromzipped.o : calc_fst_deldaf_fromzipped.c
#	$(CC) $(CCFLAG) -c calc_fst_deldaf_fromzipped.c

#calc_fst_deldaf_fromzipped : calc_fst_deldaf_fromzipped.o coal_data_tped_vers.c coal_data_tped_vers.h
#	$(CC) $(CCFLAG) -o calc_fst_deldaf_fromzipped calc_fst_deldaf_fromzipped.o coal_data_tped_vers.c


	
bootstrap_freq_popstats_regions.o : bootstrap_freq_popstats_regions.c 
	$(CC) $(CCFLAG) -c bootstrap_freq_popstats_regions.c

bootstrap_freq_popstats_regions : bootstrap_freq_popstats_regions.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o bootstrap_freq_popstats_regions bootstrap_freq_popstats_regions.o coal_data_tped_vers.c


bootstrap_ld_popstats_regions.o : bootstrap_ld_popstats_regions.c 
	$(CC) $(CCFLAG) -c bootstrap_ld_popstats_regions.c

bootstrap_ld_popstats_regions : bootstrap_ld_popstats_regions.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o bootstrap_ld_popstats_regions bootstrap_ld_popstats_regions.o coal_data_tped_vers.c


bootstrap_fst_popstats_regions.o : bootstrap_fst_popstats_regions.c 
	$(CC) $(CCFLAG) -c bootstrap_fst_popstats_regions.c

bootstrap_fst_popstats_regions : bootstrap_fst_popstats_regions.o coal_data_tped_vers.c coal_data_tped_vers.h
	$(CC) $(CCFLAG) -o bootstrap_fst_popstats_regions bootstrap_fst_popstats_regions.o coal_data_tped_vers.c


clean :
	rm *.o && rm freqs_stats && rm calc_fst_deldaf && rm bootstrap_freq_popstats_regions && rm bootstrap_ld_popstats_regions && rm bootstrap_fst_popstats_regions

install :
	cp calc_fst_deldaf /usr/local/ && cp bootstrap_freq_popstats_regions /usr/local/ && cp bootstrap_ld_popstats_regions /usr/local/ && cp bootstrap_fst_popstats_regions /usr/local/

