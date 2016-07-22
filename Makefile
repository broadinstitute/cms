## top-level Makefile for cms2.0
## last updated 07.14.16 	vitti@broadinstitute.org

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
	rm *.o



