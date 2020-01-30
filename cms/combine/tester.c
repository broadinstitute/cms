#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "cms_data.h"

//gcc -c collate_scores.c
//gcc -O0 -ggdb3 -lm -lz -Wall -o collate_scores_fromzipped collate_scores.o cms_data_zipped.c
//gcc -O0 -ggdb3 -lm -Wall -o collate_scores collate_scores.o cms_data.c

int rmdup(int *array, int length);

int main(int argc, char **argv) {
	int *positions;
	int length = 100;
	int x;
	positions = malloc(length * sizeof(int));
	positions[10] = 1;
	positions[11] = 1;
	positions[12] = 1;
	
	x = rmdup(positions, length);
	fprintf(stderr, "%d\n", x);


	fprintf(stderr, "%d\t", positions[0]);
	fprintf(stderr, "%d\t", positions[1]);
	fprintf(stderr, "%d\t", positions[2]);
	fprintf(stderr, "%d\t", positions[3]);
	fprintf(stderr, "%d\t", positions[4]);
	fprintf(stderr, "%d\t", positions[5]);
	
	return 0;
}
