from init_cond import *

# generate_ring_ic(sim_pars,ring_shape_pars,sim_name,trial_name)

R_cyl = 2.0
extra_len = 20
depl_str = 0.8
dpr = 0.05
pdiffuse = 0.05
rho_eq = 0.008
T_max = 5e6
dt = 5e-3
T_frame=1e3

sim_name = "ring3"
ex_dict = get_closed_crystals(R_cyl+1,curve_list,mode='dict')
ex_keys = ex_dict.keys()
# [R_cyl,depl_str,dpr,pdiffuse,rho_eq,T_max,dt,T_frame]=sim_pars
sim_pars = [R_cyl,extra_len,depl_str,dpr,pdiffuse,rho_eq,T_max,dt,T_frame]
for i in range(len(ex_keys)):
	ci = ex_keys[i]
	for j in range(len(ex_dict[ci])):
		ring_shape_pars = [minimal_ring,ex_dict[ci][j][0],ex_dict[ci][j][1],ex_dict[ci][j][2],R_cyl+1]
		trial_name = "t"+str(i)
		generate_ring_ic(sim_pars,ring_shape_pars,sim_name,trial_name)
