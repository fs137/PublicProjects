#ifndef DEPL_SIM_HH
#define DEPL_SIM_HH

#include <cstdio>
#ifndef NETWORKS_HH
#include "networks.hh"
#endif
#ifdef _OPENMP
#include "omp.h"
#endif

#include "common.hh"
#include "dpart.hh"
#include "gsl/gsl_rng.h"

/** The initial memory per each simulation block. */
const int init_region_memory=16;
const int init_desads_memory=128;
/** The maximum memory that can be allocated to a particular simulation. */
const int max_region_memory=65536;

/** The size of the hard repulsion force. */
const double rforce=20;

class depl_sim {
	public:
		/** The number of blocks in the x direction. */
		int m;
		/** The number of blocks in the y direction. */
		int n;
		/** The number of blocks in the z direction. */
		int o;
		/** The total number of blocks. */
		int mno;
		/** The minimum coordinate in the x direction. */
		const double ax;
		/** The maximum coordinate in the x direction. */
		const double bx;
		/** The minimum coordinate in the y direction. */
		const double ay;
		/** The maximum coordinate in the y direction. */
		const double by;
		/** The minimum coordinate in the z direction. */
		const double az;
		/** The maximum coordinate in the z direction. */
		const double bz;		
		/** The length of a block in the x direction. */
		double dx;
		/** The length of a block in the y direction. */
		double dy;
		/** The length of a block in the z direction. */
		double dz;
		/** The inverse length of a block in the x direction. */
		double xsp;
		/** The inverse length of a block in the y direction. */
		double ysp;
		/** The inverse length of a block in the z direction. */
		double zsp;
		/** The name of the directory in which to store the output. */
		const char *filename;
		const double shell_rad;
		const double shell_len;
		const double shell_wid;//=1.5;
		const double shell_lo;//=shell_rad-shell_wid;
		const double shell_losq;//=shell_lo*shell_lo;
		const double dforce;
		const double dpr;
		const double pdiffuse; // note: pdiffuse should correspond to a tenth of the diffusion constant.
		const double shell_dfconst;
		/** A counter for the total number of particles in the simulation, used
		 * to initialize the particles ID. (If particles are removed from
		 * the simulation, then n_bact may be higher than the actual number of particles. */
		int n_part;
		// A counter for the total number of particles added.
		int n_added;
		// A parameter used to label particles non-redundantly
		int max_id;
		/** A counter for the simulation frame to output. */
		int f_count;
		/** The simulation time. */
		double time;
		/** Time used to remap particles between frames. */
		double remap_time; 
		/** time used to compute forces between frames (diagnostic.) */
		double calculate_force_time; 
		double clear_force_time;
		// time to update positions between frames.
		double integrate_time; 
		/** squared radius separating desorbed and adsorbed regions.*/
		double Rdes_sq; 
		/** An array used to keep track of the number of particles per
		 * block during the remapping phase. */
		int * co;
		/** The total memory allocated to particles. */
		int mem;
		/** The total memory in each box. */
		int * bpa_mem;
		/** An array of particle information. */
		dpart *sa;
		/** A neighbor list. */
		Graph nl;
		/** Threshold for the neighbor list. */
		double thr_nl;
		double thrsq_nl;
		double R_nl;
		double interaction_thr;
		double interaction_thrsq;
		/** A list of particles in each box. */
		int ** bpa;
		/** A list of the box address associated with each particle. */
		int ** pba;
		/** The amount of time between consecutive neighbor list updates. */
		double tunl; // '(t)ime to (u)pdate (n)eighbor (l)ist.
		/** A timer for when the neighbor list should be updated. */
		double ttnu; // '(t)ime (t)oward (n)ext (u)pdate.
		int max_region_mem; // maximum memory over all regions
		bool haltflag;
		depl_sim(double ax_,double bx_,double ay_,double by_,double az_,double bz_,const char* filename_,
				double shell_rad_, double shell_len_,double dforce_, double dpr_,double pdiffuse_,int f_count_, double t_init_,double thr_nl_,int init_seed);
		~depl_sim();
		void step_forward(double dt,int pbc_flag);
		void step_forward_const_n(double dt,double lambda_in, double lambda_out, int pbc_flag);
		void integrate(double dt);
		double overlap_trace_full(const int & pbc_flag);
		double compute_overlap_energy(const int& s,const int& pbc_flag);
		double cyl_overlap_energy(const int& s);
		void calculate_forces(int pbc_flag);
		void remap(int pbc_flag);
		double pair_depl_energy(const double& rsq);
		double local_depl_energy(const int& s,const int& pbc_flag);
		void remap_const_n(double lambda_in, double lambda_out,double dt,int pbc_flag);
		void count_overlaps_nl(double epsilon, int& Noverlaps,int pbc_flag);
		void count_desorbed(int& Ndesorbed);
		void solve(double duration,int frames,double dt_max,int pbc_flag);
		void solve_const_n(double duration, int frames, double dt_max, double lambda_in, double lambda_out,int pbc_flag);
		void put(double  x,double y,double z,int pbc_flag,int id);
		void put(double x,double y,double z,int pbc_flag);
		void remove(int s);
		bool overlaps(double x,double y,double z,double r,int pbc_flag);
		void write(int k);
		void update_nl(int pbc_flag_);
		int total_particles();
		int total_particles_desads();
	private:
		/** A temporary character buffer used to create output
		 * filenames. */
		char *buf;
		/** A file handle to the diagnostic information file. */
		FILE *dfile;
		/** A file for times of snapshots */
		FILE* tfile;
		/** A record of ambient densities near the cylinder */
		FILE* density_record;
		/** The total number of threads available. */
		const int max_threads;
		/** The total number of threads to currently use. */
		int threads;
		/** An array of GSL random number generator states. */
		gsl_rng** gslr;
		gsl_rng * gslr_des; // GSL rng for desorbed region.
		void add_region_memory(int s);
		void add_sa_memory();
		double box_len_x;
		double box_len_y;
		double box_len_z;
		double domain_volume;
	public:
		double depletion_force(double x);
		double sph_int(double x,double L);
	private:
		void calculate_force(int s,int pbc_flag);
		void min_distance(double gx,double gy,double hx,double hy,double kx,double ky,double &ox,double &oy);
		/** Custom int function, that gives consistent stepping for
		 * negative numbers. With normal int, we have
		 * (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1). With this routine, we
		 * have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
		inline int step_int(double a) {
			return a<0?int(a)-1:int(a);
		}
		inline void int_box(double x,double y,double z,double r,int &li,int &ui,int &lj,int &uj,int &lk,int &uk,int pbc_flag) {
			double llim = x-ax-r;
			li=(llim>0?int(llim*xsp):int(llim*xsp)-1);
			double ulim = x-ax+r;
			ui=int(ulim*xsp);
			llim = y-ay-r;
			lj=(llim>0?int(llim*ysp):int(llim*ysp)-1);
			ulim = y-ay+r;
			uj=int(ulim*ysp);
			if(pbc_flag==0){
				/*li=int((x-ax-r)*xsp);if(li<0) li=0;if(li>=m) li=m-1;
				ui=int((x-ax+r)*xsp);if(ui<0) ui=0;if(ui>=m) ui=m-1;
				lj=int((y-ay-r)*ysp);if(lj<0) lj=0;if(lj>=n) lj=n-1;
				uj=int((y-ay+r)*ysp);if(uj<0) uj=0;if(uj>=n) uj=n-1;*/
				lk=int((z-az-r)*zsp);if(lk<0) lk=0;if(lk>=o) lk=o-1;
				uk=int((z-az+r)*zsp);if(uk<0) uk=0;if(uk>o) uk=o;
			} else if(pbc_flag==1) {
				llim = z-az-r;
				lk=(llim>0?int(llim*zsp):int(llim*zsp)-1);
				uk=int((z-az+r)*zsp);
			}
		}
		inline double random_kick() {
			return (2*double(rand())/RAND_MAX-1)*1e-3;
		}
		inline double dis_sq(double gx,double gy,double qx,double qy) {
			double ex=gx-qx,ey=gy-qy;
			return ex*ex+ey*ey;
		}
		inline double sqrt_approx(double rsq) {
			double eta = 0.25*rsq-1;
			double rval = -0.25*eta+1;
			rval=eta*rval+2;
			return rval;
		}
	public:
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return 0;}
#endif
};

#endif
