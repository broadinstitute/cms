## top-level Makefile for cms2.0
## last updated 09.09.16 	vitti@broadinstitute.org

######################
## DEFINE VARIABLES ##
######################

PREFIX ?= /usr/local
BIN_DIR = /usr/bin/

CC = gcc
CCFLAG = -O0 -ggdb3 -lm -Wall 

##################
## DEFINE RULES ##
##################

export PATH

all : model composite 

model :
	cd cms/model && $(MAKE) -e install

composite :
	cd cms/combine && $(MAKE) -e install
	
clean :
	cd cms/combine && rm *.o && rm -R *.dSYM && cd ../model && rm *.o && rm -R *.dSYM
