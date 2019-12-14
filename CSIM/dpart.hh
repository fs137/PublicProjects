#ifndef DPART_HH
#define DPART_HH

#include <cstdlib>
#include "gsl/gsl_rng.h"
#include <algorithm>
#include "common.hh"
#include "math.h"
// strength of the (unphysical) force used to keep particles within upper/lower z bounds
const double soft_constraint_force = 1.0;
const double pdamp=1;
const double augfactor = 1; 
// size of repulsion from surface.
const double shell_repulsion = 12; // originally: 12

struct dpart {
	/** The particle's ID. */
	int id;
	/** The x position of the particle. */
	double x;
	/** The y position of the particle. */
	double y;
	/** The z position of the particle. */
	double z;	
	/** The x force on the particle. */
	double fx;
	/** The y force on the particle. */
	double fy;
	/** The z force on the particle. */
	double fz;
	/** A 'state', mainly used for debugging purposes at this point. */
	int state;
	dpart() {}
	dpart(int id_,double x_,double y_,double z_) :
		id(id_), x(x_), y(y_), z(z_), state(0) {}
	inline void integrate(double dt) {
		x+=dt*fx;
		y+=dt*fy;
		z+=dt*fz;
	}
	inline double mass_inv() {
		return 1;
	}
	// Confine particle to cylindrical surface
	inline void clear_force_cyl_constr(double shell_dfconst, double shell_rad, double shell_len, double shell_lo,double shell_losq,double dpr,int pbc_flag) {
		double rsq=x*x+y*y;
		const double A = shell_rad+1;
		const double B = A+2*dpr;
		const double Bsq = B*B;
		if((rsq>shell_losq) && (rsq<Bsq)) {
			rsq=sqrt(rsq);
			double sf = rsq-A;
			if(sf>0) sf*=shell_dfconst*(rsq-B);
			else sf*=-shell_repulsion;
			fx=sf*x;
			fy=sf*y;
			fz=0;
		} else fx=fy=fz=0;
		if(pbc_flag == 0) {
			double zpsl = z+shell_len;
			double zmsl = z-shell_len;
			fz = soft_constraint_force*((zpsl<0)*((-zpsl>2)*(2+zpsl)-zpsl)-(z>shell_len)*((zmsl>2)*(2-zmsl)+zmsl));
		}
	}
	inline void clear_force_cyl(double shell_dfconst, double shell_rad, double shell_len, double shell_lo,double shell_losq,double dpr,int pbc_flag) {
		double rsq=x*x+y*y;
		const double A = shell_rad+1;
		const double B = A+2*dpr;
		const double Bsq=B*B;
		if(rsq>shell_losq && rsq<Bsq) {
			rsq=sqrt(rsq);
			double invr = 1./rsq;
			double sf = (rsq>A)*(-shell_dfconst*(B*invr-1.))+(rsq<A)*(shell_repulsion*(A*invr-1));
			fx=sf*x;
			fy=sf*y;
			fz=0;
		} else fx=fy=fz=0;
		if(pbc_flag == 0) {
			double zpsl = z+shell_len;
			double zmsl = z-shell_len;
			fz = soft_constraint_force*((zpsl<0)*((-zpsl>2)*(2+zpsl)-zpsl)-(z>shell_len)*((zmsl>2)*(2-zmsl)+zmsl));
		}
	}
	inline void add_force(double cfx,double cfy,double cfz) {
#pragma omp atomic
		fx+=cfx;
#pragma omp atomic
		fy+=cfy;
#pragma omp atomic
		fz+=cfz;
	}
	inline void random_kick(gsl_rng* gs,double kicksize) {
		double dx,dy,dz,rsq;
		do {
			dx=gsl_rng_uniform(gs)-0.5;
			dy=gsl_rng_uniform(gs)-0.5;
			dz=gsl_rng_uniform(gs)-0.5;
			rsq=dx*dx+dy*dy+dz*dz;
		} while(rsq>0.25);
		x+=dx*kicksize;
		y+=dy*kicksize;
		z+=dz*kicksize;
	}
};

#endif
