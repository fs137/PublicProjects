from init_cond import *
Nsamples = 15
# Generate a set of trials spanning a wide range of orientations at a fixed cylinder radius
R_drop = 1 # This should really be computed from nucleation theory
	   # Also, instead of selecting a fixed shape, consider choosing random droplet shapes and 
	   # 	averaging each over the lattice orientation angle
R_cyl = 2.5 # Leave a reasonable amount of room for growth
depl_str = 1.0 # This could be tuned to avoid random secondary nucleation
dpr = 0.05
pdiffuse = 0.05
rho_eq = 0.008
T_max = 2.6e5
dt = 0.001
T_frame = 260
thetas = np.linspace(0,np.pi/6,Nsamples)
cyl_pars = [R_cyl,depl_str,dpr,pdiffuse,rho_eq,T_max,dt,T_frame]
sim_name = "droplet1"
mkdir_s(sim_name)
for j in range(Nsamples):
	print "Processing sample "+str(j)
	# [R_cyl,depl_str,dpr,pdiffuse,rho_eq,T_max,dt,T_frame]
	drop_pars = [R_drop,thetas[j]]
	droplet_ic(cyl_pars,drop_pars,sim_name,trial_name="t"+str(j))
