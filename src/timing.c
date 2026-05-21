#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#define NDIM 3

#define RPTS 100000000

/**
Measure how long it takes to do certain operations, to compare what is the most efficient
*/


long int ns(struct timespec* start, struct timespec* end){
	return (end->tv_sec - start->tv_sec) * 1000000000 + (end->tv_nsec - start->tv_nsec);
}

double test_func(double input){
	return sin(input);
}

int main(void){
	struct timespec start, end;
	long long int timing = 0;

	int i,d;

	dsfmt_seed(time(NULL));

    /**** PUT GLOBAL SETUP CODE HERE ***********************/

	double number;
	double out;
	/*******************************************************/
	
	for (i = 0; i < RPTS; i++){
		/**** PUT REPEATING SETUP CODE HERE ********************/

		number = dsfmt_genrand() * 10.0;
		out = 0.0;

		/*******************************************************/
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		/**** PUT CODE TO TIME BETWEEN THESE COMMENT BLOCKS ****/

		//sqrt(number); // 21 ns
		number*number; // 21 ns
		//round(number);

//		for(d = 0; d < RPTS; d++){
//			out = test_func(number);
//		}
		

		/*******************************************************/
		clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		timing += ns(&start, &end);
	}
	printf("%lf\n", out);
	//timing /= RPTS;
	printf("Took %Li ns\n", timing/RPTS);

	return 0;
}

