#include <cstdio>
#include <cstdlib>
#include "depl_sim_nl.hh"

void elementary_ecm(const int& sign,const int& a,const int& b,int& v0,int& v1) {
	if(sign==0) {
		return;
	} else if(sign==-1) {
		v0*=a;
		v1*=b;
		return;
	} else {
		v0*=b;
		v1*=a;
		return;
	}
}

void edge_crossing_state(const int& sign_x,const int& sign_y,const int& sign_z,const int& ax,const int& bx,const int& ay,const int& by,const int& az,const  int& bz,int& v0,int&v1) {
	v0=1;
	v1=1;
	elementary_ecm(sign_x,ax,bx,v0,v1);
	elementary_ecm(sign_y,ay,by,v0,v1);
	elementary_ecm(sign_z,az,bz,v0,v1);
	return;
}

// It may be faster to implement this under a single switch statement.
void standard_crossing_state_inv(const int& cs,int& wx,int& wy,int& wz) {
	// (1,2,4), (1,5,6), (1,7,9)
	// 1, (typical)
	if(cs==1) {
		wx=0;
		wy=0;
		wz=0;
		return;
	} else if(cs<10) {
		// 2,4,5,6,7,9, (non-generic) < 10
		switch(cs) {
			case 2:
				wx=-1; wy=0;wz=0;return;
			case 4:
				wx=1; wy=0;wz=0;return;
			case 5:
				wx=0;wy=-1;wz=0;return;
			case 6:
				wx=0;wy=1,wz=0;return;
			case 7:
				wx=0;wy=0;wz=-1;return;
			case 9:
				wx=0;wy=0;wz=1;return;
		}
	} else if(cs<55) {
		// 2*5,2*6,4*5,4*6,2*7,2*9,4*7,4*9,5*7,5*9,6*7,6*9, (highly non-generic) <55
		switch(cs) {
			case 10:
				wx=-1;wy=-1;wz=0;return;
			case 12:
				wx=-1;wy=1;wz=0;return;
			case 20:
				wx=1;wy=-1;wz=0;return;
			case 24:
				wx=-1;wy=1;wz=0;return;
			case 14:
				wx=-1;wy=0;wz=-1;return;
			case 18:
				wx=-1;wy=0;wz=1;return;
			case 28:
				wx=1;wy=0;wz=-1;return;
			case 36:
				wx=-1;wy=0;wz=1;return;
			case 35:
				wx=0;wy=-1;wz=-1;return;
			case 45:
				wx=0;wy=-1;wz=1;return;
			case 42:
				wx=0;wy=1;wz=-1;return;
			case 54:
				wx=0;wy=1;wz=1;return;
		} 
	} else {
		// 2*5*7, 4*5*7, 2*6*7, 4*6*7, 2*5*9, 4*5*9, 2*6*9, 4*6*9, (extremely non-generic) > 70
		switch(cs) {
			case 70:
				wx=-1;wy=-1;wz=-1;return;
			case 140:
				wx=1;wy=-1;wz=-1;return;
			case 82:
				wx=-1;wy=1;wz=-1;return;
			case 164:
				wx=1;wy=1;wz=-1;return;
			case 90:
				wx=-1;wy=-1;wz=1;return;
			case 180:
				wx=1;wy=-1;wz=1;return;
			case 108:
				wx=-1;wy=1;wz=1;return;
			case 216:
				wx=1;wy=1;wz=1;return;
		}
	}
}

/** Constructs a depl_sim class, allocating a block structure to store the
 * particles positions and opening a diagnostic output file.
 * \param[in] (ax_,bx_) the x dimensions of the simulation region.
 * \param[in] (ay_,by_) the y dimensions of the simulation region.
 * \param[in] (az_,bz_) the z dimensions of the simulation region.
 * \param[in] (m,n,o) the number of blocks to divide the simulation region into.
 * \param[in] filename the directory name in which to store the output. */
// depl_sim bs(-40,40,-40,40,-40,40,8,8,8,buf,shell_rad,shell_dfconst,dforce,dpr);
depl_sim::depl_sim(double ax_,double bx_,double ay_,double by_,double az_,double bz_,const char* filename_,
		double shell_rad_, double shell_len_, double dforce_, double dpr_, double pdiffuse_, int f_count_,double t_init_,double thr_nl_ ,int init_seed) :
	ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
	filename(filename_), shell_rad(shell_rad_),shell_len(shell_len_),shell_wid(1.5),shell_lo(shell_rad-shell_wid),shell_losq(shell_lo*shell_lo),dforce(dforce_),dpr(dpr_),pdiffuse(pdiffuse_),
	shell_dfconst(augfactor*2*M_PI*dforce*(1+dpr)*sqrt((shell_rad+dpr)/(1+shell_rad+2*dpr))),
	n_part(0), n_added(0),max_id(0),f_count(f_count_), time(t_init_),remap_time(0), calculate_force_time(0), clear_force_time(0),integrate_time(0), Rdes_sq((2+2*dpr+shell_rad)*(2+2*dpr+shell_rad)),
	mem(1024), sa(new dpart[1024]), thr_nl(thr_nl_),thrsq_nl(thr_nl_*thr_nl_), R_nl(0.5*thr_nl), pba(new int*[1024]), tunl((thr_nl-2)*(thr_nl-2)/8), ttnu(0), buf(new char[256]),
#ifdef _OPENMP
	max_threads(omp_get_max_threads()),
#else
	max_threads(1),
#endif
	gslr(new gsl_rng*[max_threads]),gslr_des(gsl_rng_alloc(gsl_rng_taus)) {	
		int i;
	gsl_rng_set(gslr_des,max_threads+1+init_seed); // seed the rng for the desorbed region with a value different from the seeds for each processor. 
	// Set up the GSL random number generators
	for(i=0;i<max_threads;i++) {
		gslr[i]=gsl_rng_alloc(gsl_rng_taus);
		gsl_rng_set(gslr[i],init_seed+i);
	}
	haltflag = false;	
	// Initialize the simulation blocks
	double L = thr_nl;
	box_len_x = bx-ax;
	box_len_y = by-ay;
	box_len_z = bz-az;
	domain_volume=box_len_x*box_len_y*box_len_z-box_len_z*(M_PI*(shell_rad+1)*(shell_rad+1));
	interaction_thr = 2+2*dpr;
	interaction_thrsq = interaction_thr*interaction_thr;
	printf("Box dimensions: %g %g %g\n",box_len_x,box_len_y,box_len_z);
	m = int(box_len_x/L);
	if(m==0) m=1;
	n = int(box_len_y/L);
	if(n==0) n=1;
	o = int(box_len_z/L);
	if(o==0) o=1;
	mno = m*n*o;
	co = new int [mno];
	bpa = new int * [mno];
	bpa_mem=new int [mno],
	dx = (bx-ax)/m;
	xsp = 1/dx;
	dy = (by-ay)/n;
	ysp = 1/dy;
	dz = (bz-az)/o;
	zsp = 1/dz;
	printf("Spacings: %g %g %g\n",xsp,ysp,zsp);
	if(xsp<0||ysp<0||zsp<0) {
		printf("bx,ax,by,ay,bz,az=%g,%g,%g,%g,%g,%g\n",bx,ax,by,ay,bz,az);
		printf("m,n,o=%d,%d,%d\n",m,n,o);
	}
	max_region_mem = init_region_memory;
	for(i=0;i<mno;i++) {
		co[i]=0;
		bpa_mem[i]=init_region_memory;
		bpa[i]=new int[init_region_memory];
	}
	// 
	// Initialize pba:
	for(i=0;i<1024;i++) pba[i]=new int[2]; 
	// Open diagnostic file
	sprintf(buf,"%s/dfile",filename);
	if (f_count_==0) {
		dfile=safe_fopen(buf,"w");
	} else {
		dfile=safe_fopen(buf,"a");
	}
	// Open time file
	sprintf(buf,"%s/tfile",filename);
	if (f_count_==0) tfile=safe_fopen(buf,"w");
	else tfile=safe_fopen(buf,"a");
	char * nbuf = new char[256];
	sprintf(nbuf,"%s/density_record",filename);
	if (f_count_==0) density_record=safe_fopen(buf,"w");
	else density_record=safe_fopen(buf,"a");
}

/** The class destructor closes the diagnostic file and frees the dynamically
 * allocated memory. */
depl_sim::~depl_sim() {
	fclose(dfile);
	fclose(tfile);
	fclose(density_record);
	for(int i=mno-1;i>=0;i--) {
		delete [] bpa[i];
	}
	for(int i=mem-1;i>=0;i--) {
		delete [] pba[i];
	}
	for(int i=max_threads-1;i>=0;i--) gsl_rng_free(gslr[i]);
	gsl_rng_free(gslr_des);
	delete [] gslr;
	delete [] buf;
	delete [] bpa;
	delete [] pba;
	delete [] bpa_mem;
	delete [] co;
	delete [] sa;
}

/** Carries out a particlesl growth simulation
 * \param[in] duration the time interval over which to simulate.
 * \param[in] frames the number of frames to store. */
void depl_sim::solve(double duration,int frames,double dt_max,int pbc_flag) {
	double t_start=time,time_interval=duration/frames,target_time,t0,t1,t2;
	int l,tb; 
	// Output the initial fields and record initial time
	if(f_count==0) {
		write(0);
		printf("# Frame %d",0);
		fprintf(dfile,"%d %d\n",0,total_particles());
		fprintf(tfile,"0 0\n");
	}

	t0=wtime();
	for(int k=1;k<=frames;k++) {

		// Compute the target time to the next output frame
		target_time=t_start+time_interval*k;l=1;

		// Calculate threads
		threads=n_part/300+1;
		if(threads>max_threads) threads=max_threads;
		remap_time=0;
		// Carry out simulation step using the regular timestep until
		// within range of the target time
		while(time+dt_max*(1+1e-8)<target_time) {
			l++;
			step_forward(dt_max,pbc_flag);
		}

		// Carry out a final simulation step, using exactly the right
		// time step to reach the target time
		step_forward(target_time-time,pbc_flag);
		
		// Output the fields
		t1=wtime();
		write(k+f_count);

		// Print diagnostic information
		t2=wtime();tb=total_particles();
		printf("# Frame %d [%d particles, %d steps, thr=%d, com=%.6g s, write=%.6g s, remap=%.6g s]\n",k+f_count,tb,l,threads,t1-t0,t2-t1,remap_time);
		fprintf(dfile,"%d %d %d %d %.8g %.8g\n",k,tb,l,threads,t1-t0,t2-t1);
		fprintf(tfile,"%d %g\n",k+f_count,target_time);
		fflush(dfile);
		t0=t2;
	}
	f_count+=frames;
	time+=duration;
}

double depl_sim::compute_overlap_energy(const int& s,const int& pbc_flag) {
	double overlap_energy = 0;
	if (pbc_flag==0) {
		dpart &b=sa[s];
		for(int j=0;j<nl.Nnbors[s];j++) {
			int sp = nl.nbors[s][j];
			if(sp<s) continue;
			double delx,dely,delz,rsq;
			dpart &bb=sa[sp];
			// compute the net force on b from bb and its images.
			delx = bb.x-b.x;
			rsq = delx*delx;
			dely=bb.y-b.y;
			rsq += dely*dely;
			delz = bb.z-b.z;
			rsq += delz*delz;
			if(rsq<4) {
				overlap_energy+=10*(4-rsq);
			}
		}
	}	
	else if	(pbc_flag==1) {
		dpart &b=sa[s];
		// determine if b can possibly interact with images of its neighbors:
		for(int i=0;i<nl.Nnbors[s];i++) {
			int sp = nl.nbors[s][i];
			if(sp<s) continue;
			// Use the edge decoration for periodic boundary conditions
			double delx,dely,delz,rsq;
			dpart &bb=sa[sp];
			delx = bb.x-b.x;
			dely = bb.y-b.y;
			delz = bb.z-b.z;
			if(nl.edge_dec[s][i]==1) {}
			else {
				int sign_x,sign_y,sign_z;
				standard_crossing_state_inv(nl.edge_dec[s][i],sign_x,sign_y,sign_z);
				delx+=sign_x*box_len_x;
				dely+=sign_y*box_len_y;
				delz+=sign_z*box_len_z;
			}
			rsq=delx*delx;
			rsq+=dely*dely;
			rsq+=delz*delz;
			if(rsq<4) {
				overlap_energy+=10*(4-rsq);
			}
		}
	} 
	return overlap_energy;
}

double depl_sim::cyl_overlap_energy(const int& s) {
	double rsq = sa[s].x*sa[s].x+sa[s].y*sa[s].y;
	double thrsq = shell_rad+1;
	thrsq*=thrsq;
	double overlap_energy = 10*(rsq<thrsq)*(thrsq-rsq);
	return overlap_energy;
}

double depl_sim::overlap_trace_full(const int& pbc_flag) {
	double E = 0;
	#pragma omp parallel for num_threads(threads)
	// check that box assignments make sense:
	for(int s=0;s<n_part;s++) {
		E+=compute_overlap_energy(s,pbc_flag);
		E+=cyl_overlap_energy(s);
	}
	return E;
}

void depl_sim::solve_const_n(double duration,int frames,double dt_max,double lambda_in, double lambda_out,int pbc_flag) {
	double t_start=time,time_interval=duration/frames,target_time,t0,t1,t2;
	int l,tb; 
	// Output the initial fields and record initial time
	if(f_count==0) {
		write(0);
		printf("# Frame %d\n",0);
		fprintf(dfile,"%d %d\n",0,total_particles());
		fprintf(tfile,"0 %g\n",t_start);
	}
	t0=wtime();
	for(int k=1;k<=frames&&!haltflag;k++) {
		// Compute the target time to the next output frame
		target_time=t_start+time_interval*k;l=1;
		printf("%g\n",target_time);
		// Calculate threads
		threads=n_part/300+1;
		if(threads>max_threads) threads=max_threads;
		calculate_force_time=0;
		integrate_time=0;
		// Carry out simulation step using the regular timestep until
		// within range of the target time
		double lambda_indt = lambda_in*dt_max;
		while(time+dt_max*(1+1e-8)<target_time&&!haltflag) {l++;step_forward_const_n(dt_max,lambda_indt,lambda_out,pbc_flag);}
		// Carry out a final simulation step, using exactly the right
		// time step to reach the target time
		step_forward_const_n(target_time-time,lambda_in*(target_time-time),lambda_out,pbc_flag);
		if(!haltflag) {
			int Ndesorbed;
			count_desorbed(Ndesorbed);
			// Output 
			t1=wtime();
			write(k+f_count);
			double overlap_indicator = overlap_trace_full(pbc_flag);
			// Print diagnostic information
			t2=wtime();tb=total_particles();
			printf("# Frame %d [%d particles, %g overlap, %d steps, thr=%d, com=%.6g s, calc_f=%.6g, integrate=%.6g]\n",k+f_count,tb,overlap_indicator,l,threads,t1-t0,calculate_force_time,integrate_time);
			fprintf(dfile,"%d %d %d %d %.8g %.8g\n",k,tb,l,threads,t1-t0,t2-t1);
			fprintf(tfile,"%d %g\n",k+f_count,target_time);
			fprintf(density_record,"%d %g\n",k+f_count,float(n_part)/domain_volume);
			fflush(dfile);
			t0=t2;
		}
	}
	f_count+=frames;
	time+=duration;
}

/** Steps the simulation forward by a given time interval, by calculating
 * forces, integrating the particles states, and remapping the particles that have
 * cross block boundaries.
 **/
void depl_sim::step_forward(double dt,int pbc_flag) {
	if(ttnu>tunl) {
		update_nl(pbc_flag);
		ttnu=0;
	}
	calculate_forces(pbc_flag);
	integrate(dt);
	remap(pbc_flag);
	time+=dt;
	ttnu+=dt;
}
void depl_sim::step_forward_const_n(double dt,double lambda_indt, double lambda_out,int pbc_flag) {
	double wt1, wt2;
	if(ttnu<tunl); 
	else{
		update_nl(pbc_flag);
		ttnu=0;
	}
	wt1 = wtime();
	calculate_forces(pbc_flag);
	wt2 = wtime();
	calculate_force_time+=wt2-wt1;
	wt1 = wt2;
	integrate(dt);
	wt2 = wtime();
	integrate_time+=wt2-wt1;
	remap_const_n(lambda_indt,lambda_out,dt,pbc_flag);
	ttnu +=dt;
	time+=dt;
}

/** Remove a particle from the simulation. */
void depl_sim::remove(int s) {
	--n_part;
	int old_sbox = pba[s][0];
	int old_sj = pba[s][1];
		co[old_sbox]-=1;
		bpa[old_sbox][old_sj]=bpa[old_sbox][co[old_sbox]]; 
		pba[bpa[old_sbox][old_sj]][1]=old_sj; 
	if(pba[bpa[old_sbox][old_sj]][0]==old_sbox) {}
	else {
		printf("Error in remove(int): mismatch between pba and bpa\n"); 
		pba[bpa[old_sbox][old_sj]][0]=old_sbox;
	}
	if(s<n_part) {
		sa[s]=sa[n_part];
		int new_sbox = pba[n_part][0];
		int new_sj = pba[n_part][1];
		pba[s][0]=new_sbox; 
		pba[s][1]=new_sj; 
		bpa[new_sbox][new_sj]=s; 
	}
	remove_v(s,nl);
}

/** Adds a particle to the simulation.*/
void depl_sim::put(double x,double y,double z,int pbc_flag,int id) {
	int nx=int((x-ax)*xsp);if(nx<0) nx=0;if(nx>=m) nx=m-1;
	int ny=int((y-ay)*ysp);if(ny<0) ny=0;if(ny>=n) ny=n-1;
	int nz=int((z-az)*zsp);if(nz<0) nz=0;if(nz>=o) nz=o-1;
	int s=nx+m*(ny+n*nz);
	if(co[s]>=bpa_mem[s]) add_region_memory(s);
	bpa[s][co[s]]=n_part; 
	if(n_part==mem) add_sa_memory(); 
	pba[n_part][0]=s;
	pba[n_part][1]=co[s]++; 
	sa[n_part]=dpart(id,x,y,z);  
	// update the neighbor list:
	add_v(nl,id);
	max_id = (max_id>id?max_id:id+1);
	n_added++;
	int li, lj, lk, ui, uj, uk;
	if(pbc_flag==0) {
		int_box(x,y,z,R_nl,li,ui,lj,uj,lk,uk,0);
		for(int ci=li;ci<ui;ci++) {
			int ci0=ci,sign_x=0;
			if (-1<ci&&ci<m); 
			else if(ci<0) {ci0+=m; sign_x=1;}
			else {ci0-=m; sign_x=-1;}
			for(int cj=lj;cj<uj;cj++) {
				int cj0=cj, sign_y=0;
				if(-1<cj and cj<n); 
				else if(cj<0) {cj0+=n; sign_y=1;}
				else {cj0-=n; sign_y=-1;}
				for(int ck=lk;ck<uk;ck++) {
					int sp = ci0+m*(cj0+n*ck);
					for(int jp=0;jp<co[sp];jp++) {
						int pv = bpa[sp][jp];
						if(pv==n_part) continue;
						dpart bb = sa[pv];
						double delsqx = bb.x-x-sign_x*box_len_x;
						delsqx*=delsqx;
						double delsqy = bb.y-y-sign_y*box_len_y;
						delsqy*=delsqy;
						double delsqz = bb.z-z;
						delsqz*=delsqz;
						if(delsqx+delsqy+delsqz<thrsq_nl) {
							add_e(nl,pv,n_part);
						}
					}
				}
			}
		}
		n_part++;
	}
	else if(pbc_flag==1) {
		int_box(x,y,z,R_nl,li,ui,lj,uj,lk,uk,1);
		for(int ci=li;ci<ui;ci++) {
			int ci0=ci,sign_x=0;
			if (-1<ci&&ci<m); 
			else if(ci<0) {ci0+=m; sign_x=1;}
			else {ci0-=m; sign_x=-1;}
			for(int cj=lj;cj<uj;cj++) {
				int cj0=cj, sign_y=0;
				if(-1<cj and cj<n); 
				else if(cj<0) {cj0+=n; sign_y=1;}
				else {cj0-=n; sign_y=-1;}
				for(int ck=lk;ck<uk;ck++) {
					int ck0=ck, sign_z=0;
					if(-1<ck&&ck<o);
					else if(ck<0) {ck0+=o; sign_z=1;}
					else{ck0-=o; sign_z=-1;}
					int sp = ci0+m*(cj0+n*ck0);
					for(int jp=0;jp<co[sp];jp++) {
						int pv = bpa[sp][jp];
						if(pv==n_part) continue;
						dpart bb = sa[pv];
						double delsqx = bb.x-x-sign_x*box_len_x;
						delsqx*=delsqx;
						double delsqy = bb.y-y-sign_y*box_len_y;
						delsqy*=delsqy;
						double delsqz = bb.z-z-sign_z*box_len_z;
						delsqz*=delsqz;
						if(delsqx+delsqy+delsqz<thrsq_nl) {
							add_e(nl,pv,n_part);
						}
					}
				}
			}
		}
		n_part++;
	}
	else {
		printf("Error: pbc_flag must be 0 or 1.\n");
	}
}

void depl_sim::put(double x,double y,double z,int  pbc_flag) {
	put(x,y,z,pbc_flag,max_id);
}

/** Calculates the total number of particles in the simulation.*/
int depl_sim::total_particles() {
	int tot=*co;
	for(int i=1;i<mno;i++) tot+=co[i];
	if(tot!=n_part) printf("Error in total_particles(): return value should agree with n_part\n");
	return tot;
}

/** Doubles the memory allocation in a simulation block. */
void depl_sim::add_region_memory(int s) {
	bpa_mem[s]<<=1;
	if(bpa_mem[s]>=max_region_memory) {
		fprintf(stderr,"Memory allocation exceeded in region %d\n",s);
		exit(1);
	}
	int* nba=new int[bpa_mem[s]];
	for(int i=0;i<co[s];i++) nba[i]=bpa[s][i];
	delete [] bpa[s];
	bpa[s]=nba;
}

void depl_sim::add_sa_memory() {
	mem <<=1;
	dpart * new_sa = new dpart[mem];
	for(int i=0;i<n_part;i++) new_sa[i]=sa[i];
	delete [] sa;
	sa = new_sa;
	int ** new_pba = new int * [mem];
	for(int i=0;i<mem;i++) new_pba[i]=new int[2];
	for(int i=0;i<n_part;i++) {
		new_pba[i][0]=pba[i][0];
		new_pba[i][1]=pba[i][1];
		delete [] pba[i];
	}
	delete [] pba;
	pba = new_pba;
}

void depl_sim::integrate(double dt) {
	const double kicksize=sqrt(pdiffuse*dt);
#pragma omp parallel for num_threads(threads)
	for(int s=0;s<n_part;s++) {
		gsl_rng* gs=gslr[omp_get_thread_num()];	
		sa[s].integrate(dt);
		sa[s].random_kick(gs,kicksize);
	}
}

/** Calculates and stores all of the forces on all of the particles. */
void depl_sim::calculate_forces(int pbc_flag) {
#pragma omp parallel for num_threads(threads)
	for(int s=0;s<n_part;s++) {
		sa[s].clear_force_cyl(shell_dfconst,shell_rad,shell_len,shell_lo,shell_losq,dpr,pbc_flag);
	}
	// Add forces due to contacts between particles.
//#pragma omp parallel for schedule(dynamic) num_threads(threads)
#pragma omp parallel for num_threads(threads)
	// check that box assignments make sense:
	for(int s=0;s<n_part;s++) {
		calculate_force(s,pbc_flag);
	}
}

/** Calculates all of the forces that are on a particle due to contacts with
 * its neighbors.
 * \param[in] s the block that the particle is in.
 * \param[in] q the index within the block. */
void depl_sim::calculate_force(int s, int pbc_flag) {
	const double rthrsq = 4*(1+dpr)*(1+dpr);
	dpart &b = sa[s];
	switch(pbc_flag) {
		case 0:
			for(int j=0;j<nl.Nnbors[s];j++) {
				int sp = nl.nbors[s][j];
				if(sp<s) continue;
				double delx,dely,delz,rsq,cfx,cfy,cfz;
				dpart &bb=sa[sp];
				delx = bb.x-b.x;
				// compute the net force on b from bb and its images.
				rsq=delx*delx;
				if(rsq>rthrsq) continue;
				dely=bb.y-b.y;
				rsq+=dely*dely;
				if(rsq>rthrsq) continue;
				delz = bb.z-b.z;
				rsq += delz*delz;
				if(rsq>rthrsq) { 
					continue;
				}
				else if(rsq<4){
					// Avoids computing a square root
					// (Also, this case is likely when particles are densely packed.)
					rsq = rforce*(4-rsq);
					cfx=rsq*delx;
					cfy=rsq*dely;
					cfz=rsq*delz;
					bb.add_force(cfx,cfy,cfz);
					b.add_force(-cfx,-cfy,-cfz);
				} else {
					rsq=sqrt_approx(rsq);
					rsq=depletion_force(rsq)/rsq;
					cfx=rsq*delx;
					cfy=rsq*dely;
					cfz=rsq*delz;
					bb.add_force(cfx,cfy,cfz);
					b.add_force(-cfx,-cfy,-cfz);
				}
			}
			break;
		case 1:
			// determine if b can possibly interact with images of the system:
			for(int i=0;i<nl.Nnbors[s];i++) {
				int sp = nl.nbors[s][i];
				if(sp<s) continue;
				double delx,dely,delz,rsq,cfx,cfy,cfz;
				double delxsq, delysq, delzsq;
				dpart &bb=sa[sp];
				// compute the net force on b from bb and its images.
				delx = bb.x-b.x;
				dely = bb.y-b.y;
				delz = bb.z-b.z;
				if(nl.edge_dec[s][i]==1) {}
				else {
					int sign_x,sign_y,sign_z;
					standard_crossing_state_inv(nl.edge_dec[s][i],sign_x,sign_y,sign_z);
					delx+=sign_x*box_len_x;
					dely+=sign_y*box_len_y;
					delz+=sign_z*box_len_z;
				}
				delxsq=delx*delx;
				delysq=dely*dely;
				delzsq=delz*delz;
				rsq=delxsq+delysq+delzsq;
				if(rsq>rthrsq) {}
				else if(rsq<4) {
					rsq = rforce*(4-rsq);
					cfx=rsq*delx;
					cfy=rsq*dely;
					cfz=rsq*delz;
					bb.add_force(cfx,cfy,cfz);
					b.add_force(-cfx,-cfy,-cfz);
				}
				else if(rsq<rthrsq) {
					rsq=sqrt_approx(rsq);
					rsq=depletion_force(rsq)/rsq;
					cfx=rsq*delx;
					cfy=rsq*dely;
					cfz=rsq*delz;
					bb.add_force(cfx,cfy,cfz);
					b.add_force(-cfx,-cfy,-cfz);
				}
			}
			break;
	}
}

/** Checks to see whether a sphere at a given position vector overlaps with any
 * existing point.*/
bool depl_sim::overlaps(double x,double y,double z,double r,int pbc_flag) {
	double rsq = r*r;
	if(pbc_flag==0) {
		int li,ui,lj,uj,lk,uk,ci,cj,ck;
		int_box(x,y,z,r,li,ui,lj,uj,lk,uk,0);
		for(ck=lk;ck<=uk;ck++) {
			for(cj=lj;cj<=uj;cj++) {
				for(ci=li;ci<=ui;ci++) {
					int cijk=ci+m*(cj+n*ck),qq;
					double delsqx,delsqy,delsqz;
					for(qq=0;qq<co[cijk];qq++) {
						int bbv = bpa[cijk][qq];
						dpart &bb=sa[bbv];
						delsqx = bb.x-x;
						delsqy = bb.y-y;
						delsqz = bb.z-z;
						delsqx*=delsqx;
						delsqy*=delsqy;
						delsqz*=delsqz;
						if(delsqx+delsqy+delsqz<rsq) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	else if(pbc_flag==1) {
		int li,ui,lj,uj,lk,uk,ci,cj,ck;
		int_box(x,y,z,r,li,ui,lj,uj,lk,uk,1);
		//printf("\t\t (int box with limits li,ui,lj,uj,lk,uk=%d,%d,%d,%d,%d,%d\n",li,ui,lj,uj,lk,uk);
		for(ck=lk;ck<=uk;ck++) {
			int sign_z = ((ck>-1)&&(ck<o)?0:(ck<0?1:-1));
			int ck0 = (sign_z==0?ck:ck+sign_z*o);;
			for(cj=lj;cj<=uj;cj++) {
				int sign_y = ((cj>-1)&&(cj<n)?0:(cj<0?1:-1));
				int cj0 = (sign_y==0?cj:cj+sign_y*n);
				for(ci=li;ci<=ui;ci++) {
					int sign_x = ((ci>-1)&&(ci<m)?0:(ci<0?1:-1));
					int ci0 = (sign_x==0?ci:ci+sign_x*m);
					int cijk=ci0+m*(cj0+n*ck0),qq;
					double delsqx,delsqy,delsqz;
					for(qq=0;qq<co[cijk];qq++) {
						int bbv = bpa[cijk][qq];
						dpart &bb=sa[bbv];
						delsqx = bb.x-x;
						delsqx = (sign_x==0?delsqx:delsqx+box_len_x*sign_x);
						delsqx*=delsqx;
						delsqy = bb.y-y;
						delsqy = (sign_y==0?delsqy:delsqy+box_len_y*sign_y);
						delsqy*=delsqy;
						delsqz = bb.z-z;
						delsqz = (sign_z==0?delsqz:delsqz+box_len_z*sign_z);
						delsqz*=delsqz;
						if(delsqx+delsqy+delsqz<rsq) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	else {
		printf("pbc_flag must be either 0 (wall) or 1 (3-torus)\n");
		return false;
	}
}

double depl_sim::depletion_force(double x) {
	const double L=1+dpr;
	x*=0.5;
	//if(x>L) return 0; // This should be handled implicitly
	return -dforce*(L+x)*(L-x);
}

double depl_sim::sph_int(double x,double L) {
	return (L-x)*(L-x)*(2*L+x);
}

class id_x {
	public:
		int id;
		int x;
		int y;
		int z;
		id_x(int id_,int x_,int y_,int z_) {
			x=x_; y=y_; z=z_; id=id_;
		}
		id_x(): id(0),x(0),y(0),z(0){}
};

/** Saves the particle positions as a text file containing integer IDs,
 * positions, lengths, and rotations.
 * \param[in] k the integer suffix to add to the filename. */
void depl_sim::write(int k) {	
	sprintf(buf,"%s/f.%d",filename,k);
	FILE *fp=safe_fopen(buf,"w");
	//FILE * fp;
	//fp = fopen(buf,"wb");
	if(fp!=NULL) {
		// Express coordinates as a pointer
		//id_x * ptr = new id_x [n_part];
		//for(int s=0;s<n_part;s++) {
		//	ptr[s].id = sa[s].id;
		//	ptr[s].x = sa[s].x; ptr[s].y = sa[s].y; ptr[s].z = sa[s].z;
		//}
		//if(fwrite(ptr,sizeof(id_x),n_part,fp)!=n_part) {
	//		fputs("Error writing data\n",stderr);
	//		return;
	//	}
		for(int s=0;s<n_part;s++) {	
			dpart b=sa[s];	
			fprintf(fp,"%d %g %g %g\n",b.id,b.x,b.y,b.z);	
		}
		fclose(fp);
		//delete [] ptr;
	}	
}

void depl_sim::update_nl(int pbc_flag) {
	// Update bpa[][], pba[][], and nl:
	// first, sort particles into bins.
	// check that nl is sensible:
	if(pbc_flag==1) {
		for(int s=0;s<n_part;s++) {
			// Compare the speed of this to the (<&&>?(<?:)) approach
			if((sa[s].x>ax)&&(sa[s].x<bx)) {}
			else sa[s].x+=box_len_x*((sa[s].x<ax)-(sa[s].x>bx));
			if((sa[s].y>ay)&&(sa[s].y<by)) {}
			else sa[s].y+=box_len_y*((sa[s].y<ay)-(sa[s].y>by));
			if((sa[s].z>az)&&(sa[s].z<bz)) {}
			else sa[s].z+=box_len_z*((sa[s].z<az)-(sa[s].z>bz));
		}
		for(int s=0;s<n_part;s++) {
			double x = sa[s].x; // sa accessed
			double y = sa[s].y;
			double z = sa[s].z;
			int ni = int((x-ax)*xsp); if(ni<0) ni+=m; if(ni>=m) ni-=m;
			int nj = int((y-ay)*ysp); if(nj<0) nj+=n; if(nj>=n) nj-=n;
			int nk = int((z-az)*zsp); if(nk<0) nk+=o; if(nk>=o) nk-=o;
			int sbpa = ni+m*(nj+n*nk);
			if((sbpa==pba[s][0]) && (s==bpa[pba[s][0]][pba[s][1]])) continue; // pba accessed
			else {
				// move particle s to its new block,
				// but first update bpa[sbpa0][j_s] and pba[bpa[sbpa0][co[sbpa0]-1]] 
				int sbpa0 = pba[s][0]; 
				int j_s = pba[s][1]; 
				co[sbpa0]-=1;
				int sp = bpa[sbpa0][co[sbpa0]]; 
				bpa[sbpa0][j_s]=sp; 
				pba[sp][1]=j_s; 
				if(pba[sp][0]==sbpa0) {}
				else {
					// If this occurs: then particle sp was somehow reassigned to a different box
					// without updating pba, or vice versa. (Usually happens when adding new structure
					// 	to lists.)
					printf("Error: particle %d should be in cell %d but is recorded in %d\n",sp,sbpa0,pba[sp][0]);
					if(sp>s) {
						printf("note: sp=%d>%d=s\n",sp,s);
					} else {
						printf("sp <= s\n");
					}
					pba[sp][0]=sbpa0;
				}
				if(bpa_mem[sbpa]>co[sbpa]) {}
				else {
					add_region_memory(sbpa);
					if(bpa_mem[sbpa]<=max_region_mem);
					else max_region_mem=bpa_mem[sbpa];
				}
				bpa[sbpa][co[sbpa]]=s;
				pba[s][0]=sbpa;
				pba[s][1]=co[sbpa];
				co[sbpa]+=1;
			}
		}	
	} else {
		for(int s=0;s<n_part;s++) {
			double x = sa[s].x; 
			double y = sa[s].y;
			double z = sa[s].z;
			int ni = int((x-ax)*xsp); if(ni<0) ni=0; if(ni>=m) ni=m-1;
			int nj = int((y-ay)*ysp); if(nj<0) nj=0; if(nj>=n) nj=n-1;
			int nk = int((z-az)*zsp); if(nk<0) nk=0; if(nk>=o) nk=o-1;
			int sbpa = ni+m*(nj+n*nk);
			if(sbpa==pba[s][0]) continue; 
			else {
				int sbpa0 = pba[s][0]; 
				int j_s = pba[s][1]; 
				int sp = bpa[sbpa0][co[sbpa0]]; 
				co[sbpa0]-=1;
				bpa[sbpa0][j_s]=sp; 
				pba[sp][1]=j_s; 
				if(pba[sp][0]==sbpa0) {}
				else {
					printf("Error: particle %d should be in cell %d but is recorded in %d\n",sp,sbpa0,pba[sp][0]);
					if(sp>s) {
						printf("note: sp>s\n");
					} else {
						printf("sp <= s\n");
					}
					pba[sp][0]=sbpa0;
				}
				if(bpa_mem[sbpa]>co[sbpa]) {}
				else{
					add_region_memory(sbpa);
					if(bpa_mem[sbpa]>max_region_mem) max_region_mem=bpa_mem[sbpa];
				}
				bpa[sbpa][co[sbpa]]=s;
				pba[s][0]=sbpa;
				pba[s][1]=co[sbpa];
				co[sbpa]+=1;
			}
		}	
	}	
	// use the box lists to simplify construction of neighbor lists:
	for(int s=0;s<n_part;s++) nl.Nnbors[s]=0;
	for(int s=0;s<n_part;s++) {
		double x = sa[s].x;
		double y = sa[s].y;
		double z = sa[s].z;
		int li, ui, lj, uj, lk, uk;
		if(pbc_flag==1) {
			int_box(x,y,z,thr_nl,li,ui,lj,uj,lk,uk,1);
			for(int ci=li;ci<=ui;ci++) {
				int sign_x=(ci>=m)-(ci<0);
				int ci0=(sign_x==0?ci:ci-sign_x*m);
				double x_offset = x-sign_x*box_len_x;
				for(int cj=lj;cj<=uj;cj++) {
					int sign_y=(cj>=n)-(cj<0);
					int cj0=(sign_y==0?cj:cj-sign_y*n);
					double y_offset = y-sign_y*box_len_y;
					for(int ck=lk;ck<=uk;ck++) {
						int sign_z=(ck>=o)-(ck<0);
						int ck0=(sign_z==0?ck:ck-sign_z*o);
						int sbpa = ci0+m*(cj0+n*ck0);
						double z_offset = z-sign_z*box_len_z;
						for(int q=0;q<co[sbpa];q++) {
							int sp = bpa[sbpa][q];
							if(sp<=s) continue;
							double xp = sa[sp].x;
							double yp = sa[sp].y;
							double zp = sa[sp].z;
							double delx = x_offset-xp;
							double dely = y_offset-yp;
							double delz = z_offset-zp;
							if(delx*delx+dely*dely+delz*delz<thrsq_nl) {
								// Consider checking if the edge exists in nl first
								// 	Do this by skipping particles in nl
								// 	and removing edges from nl that are 
								// 	no longer salient.
								// if not, check if the particles are interacting.
								// 	(if so, increment N_aberrant)
								add_e(nl,s,sp);
								int nis = nl.Nnbors[s]-1;
								int nisp = nl.Nnbors[sp]-1;
								edge_crossing_state(sign_x,sign_y,sign_z,2,4,5,6,7,9,nl.edge_dec[s][nis],nl.edge_dec[sp][nisp]);
							}
						}
					}
				}
			}	
		}
		else if(pbc_flag==0) {
			int_box(x,y,z,thr_nl,li,ui,lj,uj,lk,uk,0);
			for(int ci=li;ci<=ui;ci++) {
				for(int cj=lj;cj<=uj;cj++) {
					for(int ck=lk;ck<=uk;ck++) {
						int sbpa = ci+m*(cj+n*ck);
						for(int q=0;q<co[sbpa];q++) {
							int sp = bpa[sbpa][q];
							if(sp<=s) continue;
							double xp = sa[sp].x;
							double yp = sa[sp].y;
							double zp = sa[sp].z;
							double delx = x-xp;
							double dely = y-yp;
							double delz = z-zp;
							if(delx*delx+dely*dely+delz*delz<thrsq_nl) {
								add_e(nl,s,sp);
							}
						}
					}
				}
			}
		}
	}
}

double depl_sim::pair_depl_energy(const double& rsq) {
	double eta = 0.25*rsq-1;
	double dr = eta*(0.5+0.125*eta); 
	double f1 = dpr-dr; 
	double f2 = 3+2*dpr+dr;
	return -0.667*dforce*f1*f1*f2+(eta<0)*(-eta)*rforce;
}

// Compute the interaction energy associated with a single particle
double depl_sim::local_depl_energy(const int& s,const int& pbc_flag) {
	double E = 0;
	double x = sa[s].x;
	double y = sa[s].y;
	double z = sa[s].z;
	switch(pbc_flag) {
		case 0:
			for(int i=0;i<nl.Nnbors[s];i++) {
				int si = nl.nbors[s][i];
				double delx = sa[si].x-x;
				double rsq = delx*delx;
				if(rsq>interaction_thrsq) continue;
				double dely = sa[si].y-y;
				rsq+=dely*dely;
				if(rsq>interaction_thrsq) continue;
				double delz = sa[si].z-z;
				rsq+=delz*delz;
				if(rsq>interaction_thrsq) continue;
				else {
					E+=pair_depl_energy(rsq);
				}
			}
			break;
		case 1:
			for(int i=0;i<nl.Nnbors[s];i++) {
				int si=nl.nbors[s][i];
				double delx=sa[si].x-x;
				double dely=sa[si].y-y;
				double delz=sa[si].z-z;
				if(nl.edge_dec[s][i]==1) {}
				else {
					int sign_x,sign_y,sign_z;
					standard_crossing_state_inv(nl.edge_dec[s][i],sign_x,sign_y,sign_z);
					delx+=sign_x*box_len_x;
					dely+=sign_y*box_len_y;
					delz+=sign_z*box_len_z;
				}
				double rsq=delx*delx+dely*dely+delz*delz;
				if(rsq>interaction_thrsq) continue;
				else {
					E+=pair_depl_energy(rsq);
				}
			}
			break;
	}
	return E;
}

// checks if particles have moved from one block to another, and also if a particle has moved across the boundary between the desorbed and adsorbed regions. If a particle is in the desorbed region, remap also deletes the particle with a small probability, and adds new particles to a free site in the desorbed region at a rate that depends on the free volume.
void depl_sim::remap_const_n(double lambda_indt,double lambda_out,double dt,int pbc_flag) {
	int s=0;
	double invdt = 1/dt;
	switch(pbc_flag) {
		case 0:
			while(s<n_part) {
				double x = sa[s].x;
				double y = sa[s].y;
				if(x*x+y*y>Rdes_sq) {
					double test = gsl_rng_uniform(gslr_des);
					if(test*invdt<lambda_out) {
						// Compute the energy needed to remove the particle
						double ecost = local_depl_energy(s,0);
						ecost = exp(ecost);
						test = gsl_rng_uniform(gslr_des);
						if(test<ecost) {
							remove(s);
						}
						// skip the remaining steps without incrementing s
						continue;
					}
				}
				s++;
			}
			break;
		case 1:	
			while(s<n_part) {
				double x = sa[s].x;
				double y = sa[s].y;
				if(x*x+y*y>Rdes_sq) {
					double test = gsl_rng_uniform(gslr_des);
					if(test*invdt<lambda_out) {
						double ecost = local_depl_energy(s,pbc_flag);
						ecost = exp(ecost);
						test = gsl_rng_uniform(gslr_des);
						if(test<ecost) {
							remove(s);
						}
						continue;
					}
				}
				s++;
			}
			break;
	}
	// Test if an attempt will be made to add a new particle to the desorbed region. 
	double test = gsl_rng_uniform(gslr_des);
	if(test<lambda_indt) {
		// Attempt away!
		double xx,yy,zz;
		xx = yy = zz = 0;
		int Ntry = 0;
		while(xx*xx+yy*yy<Rdes_sq&&Ntry<1e6) {
			xx = ax+box_len_x*gsl_rng_uniform(gslr_des);
			yy = ay+box_len_y*gsl_rng_uniform(gslr_des);
			zz = az+box_len_z*gsl_rng_uniform(gslr_des);
		}
		if(Ntry<1e6) {
			if(!overlaps(xx,yy,zz,2,pbc_flag)){
				put(xx,yy,zz,pbc_flag);
			}
		}
	}
}

/** Scans all of the points in each block, and remaps if they have moved from
 * one block to another. */
void depl_sim::remap(int pbc_flag) {
	int s=0;
	while(s<n_part) {
		double x = sa[s].x;
		double y = sa[s].y;
		double z = sa[s].z;
		if(pbc_flag==0) {
			if(x>bx||x<ax||y>by||y<ay||z>bz||z<az) {
				remove(s);
				continue;
			}
		}
		if(pbc_flag==1) {
			if(x>ax && x<bx) {}
			else if(x>bx) {
				x-=box_len_x;
			}
			else if(x<ax) {
				x+=box_len_x;
			}
			if(y>ay && y<by) {}
			else if(y>by) {
				y-=box_len_y;
			}
			else if(y<ay) {
				y+=box_len_y;
			}
			if(z>az && z<bz) {}
			else if(z>bz) {
				z-=box_len_z;
			}
			else if(z<az) {
				z+=box_len_z;
			}
			sa[s].x=x;
			sa[s].y=y;
			sa[s].z=z;
		}
		s++;
	}	
}

void depl_sim::count_desorbed(int & Ndesorbed) {
	Ndesorbed = 0;
	for(int s=0;s<mno;s++) {
		for(int ell=0;ell<co[s];ell++) {
			dpart b = sa[bpa[s][ell]];
			if(b.x*b.x+b.y*b.y>Rdes_sq) Ndesorbed++;
		}
	}
}


void depl_sim::count_overlaps_nl(double epsilon, int & Noverlaps,int pbc_flag) {
	Noverlaps = 0;
	double sepsq = 1-epsilon;
	sepsq*=4*sepsq;
	if(pbc_flag==1) {
		for(int i=0;i<n_part;i++) {
			double x = sa[i].x;
			double y = sa[i].y;
			double z = sa[i].z;
			for(int j=0;j<nl.Nnbors[i];j++) {
				int ij = nl.nbors[i][j];
				double xx = sa[ij].x;
				double yy = sa[ij].y;
				double zz = sa[ij].z;
				double delxsq0 = x-xx;
				double delxsq1 = delxsq0+box_len_x;
				double delxsq2 = delxsq0-box_len_x;
				delxsq0*=delxsq0;
				delxsq1*=delxsq1;
				delxsq2*=delxsq2;
				if(delxsq0>delxsq1) delxsq0=delxsq1;
				if(delxsq0>delxsq2) delxsq0=delxsq2;
				double delysq0 = y-yy;
				double delysq1 = delysq0+box_len_y;
				double delysq2 = delysq0-box_len_y;
				delysq0*=delysq0;
				delysq1*=delysq1;
				delysq2*=delysq2;
				if(delysq0>delysq1) delysq0=delysq1;
				if(delysq0>delysq2) delysq0=delysq2;
				double delzsq0 = z-zz;
				double delzsq1 = delzsq0+box_len_z;
				double delzsq2 = delzsq0-box_len_z;
				delzsq0*=delzsq0;
				delzsq1*=delzsq1;
				delzsq2*=delzsq2;
				if(delzsq0>delxsq1) delzsq0=delzsq1;
				if(delzsq0>delxsq2) delzsq0=delzsq2;
				double delsq = delxsq0+delysq0+delzsq0;
				if(delsq<sepsq) {
					Noverlaps++;
				}
			}
		}
	}
	else if(pbc_flag==0) {
		for(int i=0;i<n_part;i++) {
			double x = sa[i].x;
			double y = sa[i].y;
			double z = sa[i].z;
			for(int j=0;j<nl.Nnbors[i];j++) {
				int ij = nl.nbors[i][j];
				double xx = sa[ij].x;
				double yy = sa[ij].y;
				double zz = sa[ij].z;
				double delsqx, delsqy,delsqz;
				delsqx = x-xx;
				delsqx*=delsqx;
				delsqy = y-yy;
				delsqy*=delsqy;
				delsqz = z-zz;
				delsqz*=delsqz;
				double delsq = delsqx+delsqy+delsqz;
				if(delsq<sepsq) {
					Noverlaps++;
				}
			}
		}
	}
}
