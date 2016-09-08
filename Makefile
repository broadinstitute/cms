## top-level Makefile for cms2.0
## last updated 09.08.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

all : model composite

model :
	cd cms/model && $(MAKE)

composite :
	cd cms/combine && $(MAKE)
	
clean :
	cd cms/combine && rm *.o && cd ../model && rm *.o



