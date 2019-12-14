#include "common.hh"

/** Opens a file, and checks the return value to ensure that the operation
 * was successful.*/
FILE* safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) fprintf(stderr,"rmap_sim: error opening file \"%s\"",filename);
	return temp;
}

/** Function for printing fatal error messages and exiting. */
void fatal_error(const char *p,int code) {
	fprintf(stderr,"rmap_sim: %s\n",p);
	exit(code);
}

/** Calculates the argument of two-dimensional position vector. */
double argument(double x,double y) {
	const double pi=3.1415926535897932384626433832795;
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):(x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}
