#include <cstdlib>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "depl_sim_nl.hh"
#include <fstream>
#include "stdio.h"
#include <cstring>
const char filename[]="data"; 
const char pfolder[]="data/parameters";
using namespace std;
// dfname[1]: name of directory: dir -> sweep/dir
// dfname[2]: name of parameter file: pfile -> sweep/parameters/dir/pfile.par
// dfname[3,4]: frame to read from: trialname fi -> sweep/dir/[trialname]/f.[fi]
// dfname[5]: name of output trial: trialname -> sweep/dir/[trialname]
// dfname[6]: initial seed for random number generator
int main(int argc, char *dfname[]) {	
	// Create directory for output
	mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	mkdir(pfolder,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	char * buf = new char [256];
	sprintf(buf, "%s/%s/%s.par",pfolder,dfname[1],dfname[2]);
	FILE * infile = fopen(buf,"r");
	char * linebuf = new char [256];
	// NOTE: include a list of what parameters to load from 'infile' as an optional input
	// Input mode 0:
	int Npart;
	double Tsteps;
	double dt;
	double Tframe;
	double shell_rad;
	double shell_len;
	double dforce;
	double dpr;
	double pdiffuse;
	double rho_eq;
	int pbc_flag;
	double lambda_in;
	double lambda_out;
	int mode; // flag for how to parse input parameters
	int init_seed=atoi(dfname[6]); // initial seed used to set random number generators
	if(infile!=NULL) {
		// for each input parameter: (format= npart (int), shell_rad (double), dforce (double), dpr (double)), pdiffuse (double), rho_eq (double)
		int nargs = fscanf(infile,"%d\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%d\n%lf\n%lf\n%lf\n%d",&Npart,&shell_rad,&shell_len,&dforce,&dpr,&pdiffuse,&rho_eq,&pbc_flag,&Tsteps,&dt,&Tframe,&mode);
		if(nargs==12) {
			if(mode==0) {
				// pdiffuse: diffusion constant (not spherical displacement)
				// Standard parsing: the value read to pdiffuse is actually the diffusion constant:
				pdiffuse*=10;
			}
		}
		lambda_out = pdiffuse*2e-2; // since the simulation domain is extended by a buffer of width 5. 
		lambda_in = rho_eq*lambda_out*(8*(shell_len)*(shell_rad+5)*(shell_rad+5)-M_PI*2*(shell_len)*(shell_rad+2+2*dpr)*(shell_rad+2+2*dpr))
							/(4*M_PI/3*(1+rho_eq));
		printf("lambda_in: %g. lambda_out: %g\n",lambda_in,lambda_out);
		printf("Time: %g. Time per frame: %g\n",Tsteps,Tframe);
	} else {
		printf("Infile not found\n");
		delete [] buf;
		delete [] linebuf;
		return 0;
	}
	fclose(infile);
	delete [] linebuf;
	int initframe;
	if(strcmp(dfname[3],dfname[5])==0) {
		initframe = atoi(dfname[4]);
		printf("Same trial, starting at frame %d. pbc_flag=%d\n",initframe,pbc_flag);
	}	
	else {
		printf("Trial names: %s, %s",dfname[3],dfname[5]);
		initframe = 0;
	}
	// Attempt to find a time associated with the initial frame:
	sprintf(buf,"%s/%s/%s/tfile",filename,dfname[1],dfname[3]);

	// Initialize the simulation
	sprintf(buf,"%s/%s/%s",filename,dfname[1],dfname[3]);
	mkdir(buf,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	char * tfname = new char[256];
	sprintf(tfname,"%s/%s/%s/tfile",filename,dfname[1],dfname[3]);
	FILE * tfile = fopen(tfname,"r");
	double t_init=0;
	if(tfile!=NULL) {
		double tfn;
		int fn;
		while(fscanf(tfile,"%d %lf\n",&fn,&tfn)!=EOF) {
			if(fn!=initframe) {}
			else{
				t_init=tfn;
				break;
			}
		}
		fclose(tfile);
	}
	delete [] tfname;
	
	depl_sim bs(-shell_rad-5,shell_rad+5,-shell_rad-5,shell_rad+5,-shell_len,shell_len,buf,shell_rad,shell_len,dforce,dpr,pdiffuse,initframe,t_init,2.5,init_seed);
	// Add particles at coordinates from the previous frame
	sprintf(buf,"%s/%s/%s/f.%s",filename,dfname[1],dfname[3],dfname[4]);
	FILE * coordsfile;
	coordsfile = fopen(buf,"r");
	if(coordsfile!=NULL) {
		double x,y,z;
		int id;
		while(fscanf(coordsfile,"%d %lg %lg %lg",&id,&x,&y,&z)!=EOF){
			bs.put(x,y,z,pbc_flag,id);
		}
		fclose(coordsfile);
	}
	// Record the parameter file used in current traj:
	FILE * mfile;
	sprintf(buf,"%s/%s/mfile",filename,dfname[1]);
	mfile = fopen(buf,"a");
	if(mfile!=NULL){
		fprintf(mfile,"%s %s\n",dfname[3],dfname[2]);
		fclose(mfile);
	} else {
		printf("Could not open mfile for trial %s\n",dfname[5]);
	}
	// Carry out the simulation
	sprintf(buf,"%s/%s/%s",filename,dfname[1],dfname[5]);
	//bs.solve_const_n(Tsteps,Nframes,dt,lambda_in,lambda_out);
	printf("Solving over time %g and %d frames\n",Tsteps,int(Tsteps/10)+1);
	bs.solve_const_n(Tsteps,int(Tsteps/Tframe)+1,dt,lambda_in,lambda_out,pbc_flag);
	//fclose(header);
	delete [] buf;
}
