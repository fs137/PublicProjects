#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

const double pi=3.1415926535897932384626433832795;
const double b_rad=1;

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);
double argument(double x,double y);

/** Calculates a random number between 0 and 1. */
inline double rnd() {
	return double(rand())/RAND_MAX;
}

#endif
