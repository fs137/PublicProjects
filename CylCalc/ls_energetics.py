# Define curves for the line slip energy per unit height
# Each line slip site is associated with an energy cost, typically fixed
# 	(except near exceptional points.)
# If the cost is 'epsilon', and the (flattened) line slip generator is vls,
# 	then the cost per unit height is 'epsilon/vls[1]'
# 	'epsilon' is computed from the interaction potential
from regen import *

# Load interpolation data and reconstruct guess curves
#print "Initializing"
#print "Loading interpolation data"
#[Rsamples,tsamples,v2Ds,curve_features,curve_featuresc]=load_interp_data()
#print "Constructing guess curves"
#curve_list=make_guess_curves(Rsamples,tsamples,v2Ds,curve_features,curve_featuresc)

def ls_cc_theta_check():
	v2dR_curves = {}
	for j in range(len(curve_list)):
		if curve_list[j][4]==1:
			[na,nb,ocase,icase]=curve_list[j][:4]
			# Define the exact v2dR(theta) function:"
			v2dR_exact = lambda theta: scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](theta),(na,nb,theta,icase)).x
			v2dR_curves[j]=v2dR_exact
	# Try computing v2dR_exact:
	print "Testing v2dR curves"
	for j in v2dR_curves:
		[na,nb,ocase,icase]=curve_list[j][:4]
		tmin = curve_list[j][6]
		tmax = curve_list[j][7]
		tsamples = np.linspace(tmin,tmax,100)
		for jj in range(100):
			v2dR_curves[j](tsamples[jj])
		
def V_depl_cyl_min(depl_str,dpr,R_cyl):
	return -np.pi*depl_str*(1+dpr)*np.sqrt((R_cyl+dpr)/(R_cyl+1+2*dpr))*dpr**2

def V_depl_min(depl_str,dpr):
	return -2./3.*(depl_str)*(dpr**2)*(3+2*dpr)

def V_depl(r,depl_str,dpr,tol=1e-5):
	x = 0.5*r
	C = 1+dpr
	return -2./3.*depl_str*((C-x)**2)*(2*C+x)*int((x<C)and(x>1-tol))

def plot_V_depl(depl_str=1.0,dpr=0.05,k_core=40,Nsamples=1000):
	V_min = V_depl_min(depl_str,dpr)
	r_samples = np.linspace(0,2+3*dpr,Nsamples)
	V_samples = np.zeros(Nsamples)
	for j in range(Nsamples):
		if r_samples[j]<2:
			V_samples[j]= 0.5*k_core*(2-r_samples[j])**2+V_min
		else:
			V_samples[j]=V_depl(r_samples[j],depl_str,dpr)
	plt.plot(r_samples,V_samples)
	plt.plot([2,2],[3*V_min,-3*V_min],'--r')
	plt.axis([2-3*dpr,2+3*dpr,1.1*V_min,-3*V_min])
	#plt.text(2-2*dpr,-2.5*V_min,"AO depl. potential, $C_{depl}="+str(depl_str)+"$, $\delta="+str(dpr)+"$, $k_{core}="+str(k_core)+"$",fontsize=28)
	plt.xlabel("$\ell$ (units of $r_{core}$)",fontsize=16)
	plt.ylabel("Energy (units of $k_BT n_S r_{core}^3$)",fontsize=16)
	plt.show()

def V_Morse(r,sigma,epsilon):
	x = 0.5*r
	expterm = np.exp(-sigma*(x-1))
	rval = epsilon*expterm*(expterm-2)
	return rval

# CHECK THIS (THEORY)
# Energy associated with a single line slip particle (ignoring the mirror image)
def get_ls_energy(v2d,R,depl_str,dpr,V_min):
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	v02d = v2d[:2]+v2d[4:]
	v12d = v2d[4:]-v2d[:2]
	v0 = x3d(v02d,R)
	v1 = x3d(v12d,R)
	d0sq = np.dot(v0,v0)
	d1sq = np.dot(v1,v1)
	Esample = 2.5*V_min # Energy ignoring the missing bonds
	if d1sq<rthrsq:
		Esample+=0.5*V_depl(np.sqrt(d1sq),depl_str,dpr) 
	if d0sq<rthrsq:
		Esample+=0.5*V_depl(np.sqrt(d0sq),depl_str,dpr)
	return Esample

def get_ls_tension(v2d,R,depl_str,dpr,V_min):
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	v02d = v2d[:2]+v2d[4:]
	v12d = v2d[4:]-v2d[:2]
	v0 = x3d(v02d,R)
	v1 = x3d(v12d,R)
	d0sq = np.dot(v0,v0)
	d1sq = np.dot(v1,v1)
	Esample = -V_min # Energy cost of the missing bond (typical, when dpr is small)
	if d1sq<rthrsq:
		Esample+=V_depl(np.sqrt(d1sq),depl_str,dpr) 
	if d0sq<rthrsq:
		Esample+=V_depl(np.sqrt(d0sq),depl_str,dpr)
	Esample/=abs(v2d[1])
	return Esample

# CHECK THIS (THEORY)
def get_ls_energy_theta(v2dR,theta,depl_str,dpr,V_min):
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	v02d = v2dR[:2]+v2dR[2:4]
	v12d = v2dR[2:4]-v2dR[:2]
	v0 = x3d(v02d,v2dR[4])
	v1 = x3d(v12d,v2dR[4])
	d0sq = np.dot(v0,v0)
	d1sq = np.dot(v1,v1)
	gamma = 2.5*V_min # From five nearest neighbors
	if d1sq>rthrsq:
		pass
	else:
		gamma+=0.5*V_depl(np.sqrt(d1sq),depl_str,dpr)
	if d0sq>rthrsq:
		pass
	else:
		gamma+=0.5*V_depl(np.sqrt(d0sq),depl_str,dpr)
	return gamma

def get_ls_tension_theta(v2dR,theta,depl_str,dpr,V_min):
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	v02d = v2dR[:2]+v2dR[2:4]
	v12d = v2dR[2:4]-v2dR[:2]
	v0 = x3d(v02d,v2dR[4])
	v1 = x3d(v12d,v2dR[4])
	d0sq = np.dot(v0,v0)
	d1sq = np.dot(v1,v1)
	gamma = -V_min
	if d1sq<rthrsq:
		gamma+=V_depl(np.sqrt(d1sq),depl_str,dpr)
	if d0sq<rthrsq:
		gamma+=V_depl(np.sqrt(d0sq),depl_str,dpr)
	gamma/=abs(v2dR[1])
	return gamma

# CHECK THIS (THEORY)
# Problem: in the gapped or widened line slip case, the energy
#	tends to lie lower than the maximal packing bound. 
#	This could be from adding an extra particle in each layer,
#		or from miscalculating the width between adjacent layers
#		
# Energy per unit height of a crystal with line slip defect
def get_domain_tension(v2d,ns,R,depl_str,dpr,V_min,V_min_cyl,tol=1e-5):
	# Compute three tensions, associated with each type of line slip defect
	[na,nb]=ns
	va = v2d[:2]
	vb = v2d[2:4]
	vls = v2d[4:]
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	# NOTE: shouldn't matter if nb is zero or not
	if nb!=0:
		gamma_domain = (3*abs(nb)+2)*V_min
		# NOTE: choice of sign depends on relative dot products of va, vls
		v02d = v2d[:2]+v2d[4:]
		v12d = v2d[4:]-v2d[:2]
		v0 = x3d(v02d,R)-np.array([R,0,0])
		v1 = x3d(v12d,R)-np.array([R,0,0])
		d0sq = np.dot(v0,v0)
		d1sq = np.dot(v1,v1)
		if d0sq<rthrsq:
			gamma_domain+=V_depl(np.sqrt(d0sq),depl_str,dpr)
		if d1sq<rthrsq:
			gamma_domain+=V_depl(np.sqrt(d1sq),depl_str,dpr)
		gamma_domain += (abs(nb)+1)*V_min_cyl # Contribution from interaction with the cylinder
		gamma_domain/=abs(v2d[1])
	else:
		anom = 0
		rthrsq = 2+2*dpr
		rthrsq*=rthrsq
		v02d = v2d[:2]+v2d[4:]
		v12d = v2d[4:]-v2d[:2]
		v0 = x3d(v02d,R)-np.array([R,0,0])
		v1 = x3d(v12d,R)-np.array([R,0,0])
		d0sq = np.dot(v0,v0)
		d1sq = np.dot(v1,v1)
		if d0sq<rthrsq:
			anom+=V_depl(np.sqrt(d0sq),depl_str,dpr)
		if d1sq<rthrsq:
			anom+=V_depl(np.sqrt(d1sq),depl_str,dpr)
		gamma_domain = (2*V_min+V_min_cyl+anom)*abs(na)
		# Determine which vector from v2dR corresponds to the line slip
		gamma_domain/=abs(v2d[5])
	return gamma_domain

def get_domain_tension_theta(v2dR,ns_oi,theta,depl_str,dpr,V_min,V_min_cyl,tol=1e-5):
	[na,nb,ocase,icase]=ns_oi
	v2d=get_v2d_of_v2dR(v2dR,theta,[ocase,icase])
	R=v2dR[4]
	va = v2d[:2]
	vb = v2d[2:4]
	vls = v2d[4:]
	rthrsq = 2+2*dpr
	rthrsq*=rthrsq
	if nb!=0:
		gamma_domain = (3*abs(nb)+2)*V_min
		# Include possible contribution from nearby line slip sites
		rthrsq = 2+2*dpr
		rthrsq*=rthrsq
		v02d = v2dR[:2]+v2dR[2:4]
		v12d = v2dR[2:4]-v2dR[:2]
		v0 = x3d(v02d,v2dR[4])-np.array([v2dR[4],0,0])
		v1 = x3d(v12d,v2dR[4])-np.array([v2dR[4],0,0])
		d0sq = np.dot(v0,v0)
		d1sq = np.dot(v1,v1)
		if d0sq<rthrsq:
			gamma_domain+=V_depl(np.sqrt(d0sq),depl_str,dpr)
		if d1sq<rthrsq:
			gamma_domain+=V_depl(np.sqrt(d1sq),depl_str,dpr)
		gamma_domain+=(abs(nb)+1)*V_min_cyl
		gamma_domain/=abs(v2dR[1])
	else:
		# Every site belongs to a line slip defect
		anom = 0
		rthrsq = 2+2*dpr
		rthrsq*=rthrsq
		v02d = v2dR[:2]+v2dR[2:4]
		v12d = v2dR[2:4]-v2dR[:2]
		v0 = x3d(v02d,v2dR[4])-np.array([v2dR[4],0,0])
		v1 = x3d(v12d,v2dR[4])-np.array([v2dR[4],0,0])
		d0sq = np.dot(v0,v0)
		d1sq = np.dot(v1,v1)
		if d0sq<rthrsq:
			anom+=V_depl(np.sqrt(d0sq),depl_str,dpr)
		if d1sq<rthrsq:
			anom+=V_depl(np.sqrt(d1sq),depl_str,dpr)
		gamma_domain = (2*V_min+V_min_cyl+anom)*abs(na)
		# Determine which vector from v2dR corresponds to the line slip
		gamma_domain/=abs(v2dR[3])
	return gamma_domain

# CHECK THIS
# Note on usage:
#	Nsamples = number of samples per curve (independent of the length of each curve)
def compute_domain_tensions(curve_list,depl_str,dpr,Nsamples,CR_INDEX=""):
	# The tension per unit height of a domain is also greater if particles are 
	# packed more densely.
	mkdir_s("CR_TENSION"+CR_INDEX)
	CR_PFNAME = "CR_TENSION"+CR_INDEX+"/PFILE"
	np.savetxt(CR_PFNAME,np.array([depl_str,dpr]))
	rthrsq = (2+2*dpr)
	rthrsq*=rthrsq
	V_min = V_depl_min(depl_str,dpr)
	for j in range(len(curve_list)):
		[na,nb,ocase,icase,branched]=curve_list[j][:5]
		Esamples = np.zeros(Nsamples)
		Rsamples = np.zeros(Nsamples)
		csamples = np.zeros(Nsamples,dtype=int)
		print "Computing for curve "+str(j)+" ("+str(curve_list[j][:5])+")"
		if curve_list[j][4]==0:
			# curve parametrized by radius.
			R0 = curve_list[j][6]
			R1 = curve_list[j][7]
			Rsamples = np.linspace(R0,R1,Nsamples) # This may require adjustment for extrapolation
			tsamples = np.zeros(Nsamples)
			for jj in range(Nsamples):
				V_min_cyl = V_depl_cyl_min(depl_str,dpr,R_cyl=Rsamples[jj]-1)
				# Compute the line slip vector, and the associated energy
				v2d_guess = curve_list[j][5](Rsamples[jj]) # initial guess
				v2dr = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,Rsamples[jj],[ocase,icase]))
				if v2dr.success:
					v2d = v2dr.x
					v0 = get_v0_of_vabls(v2d,ocase,icase) # 'outer case', 'inner case' for line slip branches
					tsamples[jj] = np.arccos(0.5*v0[1])
					gammas =  get_domain_tension(v2d,[na,nb],Rsamples[jj],depl_str,dpr,V_min,V_min_cyl)
					Esamples[jj] =gammas
		elif curve_list[j][4]==1:
			# parametrized by theta: solve for both R(theta) and E(theta)
			theta0 = curve_list[j][6]
			theta1 = curve_list[j][7]
			tsamples = np.linspace(theta0,theta1,Nsamples)
			for jj in range(Nsamples):
				v2dR_guess = curve_list[j][5](tsamples[jj])
				v2dRr = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,tsamples[jj],icase))
				if v2dRr.success:
					v2dR = v2dRr.x
					V_min_cyl = V_depl_cyl_min(depl_str,dpr,R_cyl=v2dR[4]-1)
					Rsamples[jj] = v2dR[4]
					v2d = get_v2d_of_v2dR(v2dR,tsamples[jj],[ocase,icase])
					Esamples[jj] = get_domain_tension_theta(v2dR,[na,nb,ocase,icase],tsamples[jj],depl_str,dpr,V_min,V_min_cyl)
					#gamma=get_domain_tension(v2d,[na,nb],Rsamples[jj],depl_str,dpr,V_min,V_min_cyl)
					#Esamples[jj]=gamma
					#csamples[jj]=gammas[-1]
		max_energy = np.max(Esamples)
		print "Maximum energy for curve "+str(j)+" is "+str(max_energy)
		if max_energy>0:
			print str(na)
		plot_data = np.zeros([Nsamples,3])
		plot_data[:,0] = Rsamples
		plot_data[:,1] = tsamples
		plot_data[:,2] = Esamples
		np.savetxt("CR_TENSION"+str(CR_INDEX)+"/RTE_curve"+str(j)+".dat",plot_data)
		#np.savetxt("CR_TENSION"+str(CR_INDEX)+"/cases"+str(j)+".dat",csamples)

def depl_ERcurve0(s,pars):
	[na,nb,ocase,icase,depl_str,dpr,V_min,j]=pars
	# parametrized by radius
	v2d = scipy.optimize.root(ls_cconstraint,curve_list[j][5](s),(na,nb,s,[ocase,icase])).x
	V_min_cyl = V_depl_cyl_min(depl_str,dpr,s) # Replace by generic parameter function
	return get_domain_tension(v2d,[na,nb],s,depl_str,dpr,V_min,V_min_cyl)

def depl_ERcurve1(s,pars):
	# parametrized by theta
	[na,nb,ocase,icase,depl_str,dpr,V_min,j]=pars
	v2dR = scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](s),(na,nb,s,icase)).x
	V_min_cyl = V_depl_cyl_min(depl_str,dpr,v2dR[4])
	v2d = get_v2d_of_v2dR(v2dR,s,[ocase,icase])
	E_depl = get_domain_tension(v2d,[na,nb],v2dR[4],depl_str,dpr,V_min,V_min_cyl)
	return np.array([v2dR[4],E_depl])

def N_excess_ls_nbors(v2d,R):
	v02d = v2d[:2]+v2d[4:]
	v12d = v2d[4:]-v2d[:2]
	v0 = x3d(v02d,R)-np.array([R,0,0])
	v1 = x3d(v12d,R)-np.array([R,0,0])
	d0sq = np.dot(v0,v0)
	d1sq = np.dot(v1,v1)
	return int(d0sq<rthrsq)+int(d1sq<rthrsq)

def solve_new_contact_transitions(CR_INDEX,tol=1e-4):
	# For each curve, determine the points at which a new contact appears/vanishes across
	#	the line slip defect
	# Open CR_INDEX parameters:
	CR_PFNAME = "CR_TENSION"+CR_INDEX+"/PFILE"
	try:
		cr_pars = np.loadtxt(CR_PFNAME)
		depl_str = cr_pars[0]
		dpr = cr_pars[1]
	except:
		depl_str = 1.0
		dpr = 0.05
	rthr = 2+2*dpr
	rthrsq = rthr*rthr
	V_min = V_depl_min(depl_str,dpr)
	nc_segments = []
	for j in range(len(curve_list)):
		[na,nb,ocase,icase]=curve_list[j][:4]
		# Solve for the radii and orientations at which the line slip defect develops or loses a 
		#	transverse interaction 
		# Starting from a commensurate lattice, measure the intervals over which the line slip defect loses or gains a contact
		ncs = [curve_list[j][6]]
		Nsamples = int((curve_list[j][7]-curve_list[j][6])/tol)+1
		samples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
		# NOTE: initially there should be one excess contact
		N_excess0 = 1
		if curve_list[j][4]==0:
			for jj in range(1,Nsamples):
				v2d_guess = curve_list[j][5](samples[jj])
				v2d = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,samples[jj],[ocase,icase])).x
				# Compute the number of excess contacts
				N_excess1 = N_excess_ls_nbors(v2d,samples[jj])
				if N_excess1!=N_excess0:
					ave_excess = 0.5*(N_excess1+N_excess0)
					# Solve for the exact crossing
					v2d_gen = lambda R: scipy.optimize.root(ls_cconstraint,curve_list[j][5](R),(na,nb,R,[ocase,icase])).x
					crossing_finder = lambda R: N_excess_ls_nbors(v2d_gen(R),R)-ave_excess
					cr_exact=scipy.optimize.brentq(crossing_finder,samples[jj-1],samples[jj])
					ncs.append(cr_exact)
					# Update the number of excess transverse interactions
					N_excess0 = N_excess1
		else:
			for jj in range(1,Nsamples):
				v2dR_guess = curve_list[j][5](samples[jj])
				v2dR=scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,samples[jj],icase)).x
				v2d = get_v2d_of_v2dR(v2dR,samples[jj],[ocase,icase])
				N_excess1 = N_excess_ls_nbors(v2d,v2dR[4])
				if N_excess1!=N_excess0:
					ave_excess = 0.5*(N_excess1+N_excess0)
					# Solve for the exact crossing
					v2dR_gen = lambda t: scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](t),(na,nb,t,icase)).x
					N_excess_gen = lambda v2dR,t: N_excess_ls_nbors(get_v2d_of_v2dR(v2dR,t,[ocase,icase]),v2dR[4])
					crossing_finder = lambda t: N_excess_gen(v2dR_gen(t),t)-ave_excess
					cr_exact=scipy.optimize.brentq(crossing_finder,samples[jj-1],samples[jj])
					ncs.append(cr_exact)
					# Update the number of excess transverse interactions
					N_excess0 = N_excess1
		nc_segments.append(ncs)
	# NOTE: the modulus of continuity can be studied independently within each segment.
	return nc_segments
	
# NOTE: avoid trivial intersections, near commensurate crystal orientations.
def compute_dt_intersections(CR_INDEX="",tol=1e-6,Rmax=1e99,printout=None):
	# Load RTE data and define approximate curves
	# Determine which can feasibly intersect
	# Use the approximate curves to define exact curves
	# For each pair of exact curves with feasible intersection,
	#	-define a joint function of two variables to R^2 (E1(s1)-E2(s2),R1(s1)-R2(s2))
	# 	-solve for roots within the bounds of each curve.
	# Rmax: bound on domain over which curve_list is exhaustive (optional.)
	RE_bnds = []
	feasible = graph()
	# Load the parameters for CR_INDEX:
	# 	(use default values if a parameter file doesn't exist)
	try:
		CR_PFNAME = "CR_TENSION"+str(CR_INDEX)+"/PFILE"
		CR_PARS = np.loadtxt(CR_PFNAME)
		depl_str = CR_PARS[0]
		dpr = CR_PARS[1]
	except:
		depl_str = 1.0
		dpr = 0.05
	V_min = V_depl_min(depl_str,dpr)
	for j in range(len(curve_list)):
		# NOTE: there is no need for interpolation in for this function
		RTEdat = np.loadtxt("CR_TENSION"+str(CR_INDEX)+"/RTE_curve"+str(j)+".dat")
		RE_bnds.append([[np.min(RTEdat[:,0]),np.max(RTEdat[:,0])],[np.min(RTEdat[:,2]),np.max(RTEdat[:,2])]])
		feasible.add_node(j)
		# Compare bounds in R and E space with existing nodes in 'feasible'
		# (when there are many curves to compare, consider sorting curves into bounding boxes first)
		for jj in range(j):
			# Check if jj is dual to j:
			if jj==curve_list[j][16]:
				continue
			if box_overlap(RE_bnds[jj],RE_bnds[j]):
				feasible.add_edge(j,jj)
	# Now that guess curves have been defined, determine the exact (R,E) curves
	#	Define v2d(R) or v2dR(theta) exactly
	#	Define V_min_cyl(R) or V_min_cyl(theta) exactly
	#	Define E(R) or ER(theta) 
	ER_curves_ex = []
	intersecting = graph()
	for j in range(len(curve_list)):
		[na,nb,ocase,icase]=curve_list[j][:4]	
		isect_pars = [na,nb,ocase,icase,depl_str,dpr,V_min,j]
		ER_curves_ex.append([isect_pars,curve_list[j][4]])
		intersecting.add_node(j)
		existence_tol = np.sqrt(tol)
		print "Checking intersections with curve "+str(j)
		for nnj in feasible.nodes[j].nbors:
			# Look for an intersection between curves j and nnj, or determine if one doesn't exist
			# (seek existence first)
			# NOTE: save time by keeping track of the duals to each new edge
			if (nnj<j):
				# First check if intersections for the mirror image of node nnj has been defined for node j
				nnj_ = curve_list[nnj][16]
				if nnj_ in intersecting.nodes[j].nbors:
					# Add nnj as a neighbor and define roots directly from intersecting.edges[(nnj_,j)]
					roots = []
					roots_ = intersecting.edges[(nnj_,j)]
					for jj in range(len(roots_)):
						Rjj = roots_[jj][0]
						# Define theta associated with the complementary curve nnj:
						if curve_list[nnj][4]==0:
							v2d_exact = scipy.optimize.root(ls_cconstraint,curve_list[nnj][5](Rjj),(na,nb,Rjj,[ocase,icase])).x
							v0 = get_v0_of_vabls(v2d_exact,ocase,icase)
							theta = np.arcsin(0.5*v0[1])
							roots.append((Rjj,theta,roots_[jj][2]))
						else:
							# Need to solve for the complementary theta at Rjj, with an appropriate guess.
							# (e.g. try to solve curve_list[nnj][5](theta)=R for the initial guess, and then try to refine
							#	it with an auxiliary solver)
							v2dR_exact = lambda theta: scipy.optimize.root(ls_cconstraint_theta,curve_list[nnj][5](theta),(na,nb,theta,icase)).x
							R_criterion = lambda theta: v2dR_exact(theta)[4]-Rjj
							theta_Rjj1 = scipy.optimize.brentq(R_criterion,curve_list[nnj][6],curve_list[nnj][-3])
							theta_Rjj2 = scipy.optimize.brentq(R_criterion,curve_list[nnj][-3],curve_list[nnj][7])
							ERjj1 = depl_ERcurve1(theta_Rjj1,[na,nb,ocase,icase,depl_str,dpr,V_min,nnj])
							ERjj2 = depl_ERcurve1(theta_Rjj2,[na,nb,ocase,icase,depl_str,dpr,V_min,nnj])
							diff1 = (ERjj1-roots_[jj][2])**2
							diff2 = (ERjj2-roots_[jj][2])**2
							if diff1<diff2:
								roots.append((Rjj,theta_Rjj1,roots_[jj][2]))
							else:
								roots.append((Rjj,theta_Rjj2,roots_[jj][2]))
							
					intersecting.add_edge(j,nnj)
					intersecting.edges[(nnj,j)]=roots
					continue
				# Existence portion: find the point of closest approach between the two curves,
				# or compute differences between the curves along a complete set of shared parametrizations
				# Solve for the exact root (if any). If a minimization approach is used, seek a root if the 
				#	minimum distance is below a certain threshold.
				if (curve_list[j][4]==0) and (curve_list[nnj][4]==0):
					# Find the common domain
					# Look for roots of Ej(R)-Ennj(R) within the shared domain
					common_domain = [max(curve_list[j][6],curve_list[nnj][6]),min(curve_list[j][7],curve_list[nnj][7])]
					diffEjEnnj = lambda R: depl_ERcurve0(R,ER_curves_ex[j][0])-depl_ERcurve0(R,ER_curves_ex[nnj][0])
					# Seek existence of a root up to a specified tolerance
					Nivals = int((common_domain[1]-common_domain[0])/existence_tol)+1
					Rsamples = np.linspace(common_domain[0],common_domain[1],Nivals)
					crossings = []
					E0 = diffEjEnnj(Rsamples[0])
					for jj in range(1,Nivals):
						E1 = diffEjEnnj(Rsamples[jj])
						if np.sign(E0)!=np.sign(E1):
							crossings.append(jj-1)
						E0 = E1
					# Define a list of exact roots associated with crossings
					# (each crossing is automatically associated with a unique root)
					roots = []
					for jji in range(len(crossings)):
						jj0 = crossings[jji]
						R0 = Rsamples[jj0]
						R1 = Rsamples[jj0+1]
						R_ = scipy.optimize.brentq(diffEjEnnj,R0,R1)
						# Include the associated angle (from v0_of_v2d)
						v2d_guess = curve_list[j][5](R_)
						v2d_exact = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,R_,[ocase,icase])).x
						v0 = get_v0_of_vabls(v2d_exact,ocase,icase)
						theta = np.arcsin(0.5*v0[1])
						roots.append((R_,theta,depl_ERcurve0(R_,ER_curves_ex[j][0])))
					if len(roots)>0:
						intersecting.add_edge(j,nnj)
						intersecting.edges[(nnj,j)]=roots
						
				elif curve_list[j][4]==0 and curve_list[nnj][4]==1:
					print "Tricky case (1st type)"
					[Rj0,Rj1]= curve_list[j][6:8]
					[t0,t1]=curve_list[nnj][6:8]
					guess_tol = np.sqrt(tol) # Can tune prefactor
					NRsamples = int((Rj1-Rj0)/guess_tol)+1
					Ntsamples = int((t1-t0)/guess_tol)+1
					Rsamples = np.linspace(Rj0,Rj1,NRsamples)
					tsamples = np.linspace(t0,t1,Ntsamples)
					# Define a two dimensional function of R_j and theta_nnj to (E_j(R_j)-E_nnj(theta_nnj),R_nnj(theta_nnj)-R_j)
					# Solve for roots with Newton's method, or with an appropriate 2D generalization.
					isect_test_fun = lambda Rtheta: depl_ERcurve1(Rtheta[1],ER_curves_ex[nnj][0])-np.array([Rtheta[0],depl_ERcurve0(Rtheta[0],ER_curves_ex[j][0])])
					# Construct a list of initial guesses for roots, based on the norm of isect_test_fun and eventually on an appropriate modulus of continuity
					init_guesses = []
					# Choose samples according to the relevant modulus of continuity
					for jj in range(NRsamples):
						for jjj in range(Ntsamples):
							guess = isect_test_fun(Rsamples[jj],tsamples[jjj])
							if np.linalg.norm(guess)<guess_tol:
								init_guesses.append([np.array([Rsamples[jj],tsamples[jjj]]),np.linalg.norm(guess)])
					# Seek a root near each guess: start
					roots = []
					for jj in range(len(init_guesses)):
						Nattempts = max(int(10./init_guesses[jj][1])+1,10)
						for jjj in range(Nattempts):
							guess = init_guesses[jj][0]+np.random.normal(0,guess_tol,2)
							Rt_result = scipy.optimize.root(isect_test_fun,guess,tol=tol)
							if Rt_result.success:
								# Check that Rt_result.x is actually in the domain of the curve
								if (curve_list[j][6]<Rt_result.x[0]<curve_list[j][7]) and (curve_list[nnj][6]<Rt_result.x[1]<curve_list[nnj][7]):
									pass
								else:
									break
								# An exact intersection has been found
								# (NOTE: the root may have been found previously)
								for j4 in range(len(roots)):
									delRt = np.array([roots[j4][0],roots[j4][1]])-Rt_result.x
									if np.dot(delRt,delRt)>tol:
										pass
									else:
										break
								Rt_root = Rt_result.x
								roots.append((Rt_root[0],Rt_root[1],depl_ERcurve0(Rt_root[0],ER_curves_ex[j][0])))
							break
					if len(roots)>0:
						intersecting.add_edge(j,nnj)
						intersecting.edges[(nnj,j)]=roots
				elif curve_list[j][4]==1 and curve_list[nnj][4]==0:
					print "Tricky case (2nd type)"
					# NOTE: partition the domain into blocks, and make several attempts within each. 
					#	(this can be tricky, because of rapid changes in modulus of continuity near commensurate lattice sites.)
					# 	Instead of making an exhorbitant number of attempts, try to adjust the domain partitions so that a few
					#	attempts in each will produce a root, if any. Also, attempt minimization from within a given block in advance 
					#	of root finding
					[Rj0,Rj1]= curve_list[nnj][6:8]
					[t0,t1]=curve_list[j][6:8]
					guess_tol = 0.2*np.sqrt(np.sqrt(tol)) # Can tune prefactor
					NRsamples = int((Rj1-Rj0)/guess_tol)+1
					Ntsamples = int((t1-t0)/guess_tol)+1
					Rsamples = np.linspace(Rj0,Rj1,NRsamples)
					tsamples = np.linspace(t0,t1,Ntsamples)
					# Define a two dimensional function of R_j and theta_nnj to (E_j(R_j)-E_nnj(theta_nnj),R_nnj(theta_nnj)-R_j)
					# Solve for roots with Newton's method 
					isect_test_fun = lambda Rtheta: depl_ERcurve1(Rtheta[1],ER_curves_ex[j][0])-np.array([Rtheta[0],depl_ERcurve0(Rtheta[0],ER_curves_ex[nnj][0])])
					isect_skep_fun = lambda Rtheta: np.linalg.norm(isect_test_fun(Rtheta))
					# Construct a list of initial guesses for roots, based on the norm of isect_test_fun and an appropriate modulus of continuity
					init_guesses = []
					# Choose samples according to the relevant modulus of continuity
					print "Choosing guesses (from "+str(NRsamples*Ntsamples)+" sample points)"
					for jj in range(NRsamples):
						for jjj in range(Ntsamples):
							guess = isect_test_fun(np.array([Rsamples[jj],tsamples[jjj]]))
							if np.linalg.norm(guess)<guess_tol:
								init_guesses.append([np.array([Rsamples[jj],tsamples[jjj]]),np.linalg.norm(guess)])
					# Seek a root near each guess:
					print "Finding roots near each guess (out of "+str(len(init_guesses))+" total)"
					roots = []
					for jj in range(len(init_guesses)):
						Nattempts = max(int(10./init_guesses[jj][1])+1,10)
						for jjj in range(Nattempts):
							guess = init_guesses[jj][0]+np.random.normal(0,guess_tol,2)
							Rt_result = scipy.optimize.root(isect_test_fun,guess)
							if Rt_result.success:
								# An exact intersection has been found
								if (curve_list[nnj][6]<Rt_result.x[0]<curve_list[nnj][7]) and (curve_list[j][6]<Rt_result.x[1]<curve_list[j][7]):
									pass
								else:
									break
								for j4 in range(len(roots)):
									delRt = np.array([roots[j4][0],roots[j4][1]])-Rt_result.x
									if np.dot(delRt,delRt)>guess_tol:
										pass
									else:
										break
								Rt_root = Rt_result.x
								roots.append([Rt_root[0],Rt_root[1],depl_ERcurve1(Rt_root[1],ER_curves_ex[j][0])])
								break
					if len(roots)>0:
						print "(Adding a new edge between curves "+str([j,nnj])+" associated with "+str(len(roots))+" intersections)"
						intersecting.add_edge(j,nnj)
						intersecting.edges[(nnj,j)]=roots
				elif curve_list[j][4]==1 and curve_list[nnj][4]==1:
					# Similar to first case, except in theta space instead of R space
					# (solve for R as well though)
					# Find the common domain
					# Look for roots of Ej(R)-Ennj(R) within the shared domain
					common_domain = [max(curve_list[j][6],curve_list[nnj][6]),min(curve_list[j][7],curve_list[nnj][7])]
					diffEjEnnj = lambda t: depl_ERcurve1(t,ER_curves_ex[j][0])-depl_ERcurve1(t,ER_curves_ex[nnj][0])
					# Seek existence of a root up to a specified tolerance
					Nivals = int((common_domain[1]-common_domain[0])/existence_tol)+1
					tsamples = np.linspace(common_domain[0],common_domain[1],Nivals)
					crossings = []
					E0 = diffEjEnnj(tsamples[0])
					for jj in range(1,Nivals):
						E1 = diffEjEnnj(tsamples[jj])
						if np.sign(E0)!=np.sign(E1):
							crossings.append(jj-1)
						E0 = E1
					# Define a list of exact roots associated with crossings
					# (each crossing is automatically associated with a unique root)
					roots = []
					for jji in range(len(crossings)):
						jj0 = crossings[jji]
						t0 = tsamples[jj0]
						t1 = tsamples[jj0+1]
						t_ = scipy.optimize.brentq(diffEjEnnj,t0,t1)
						v2dR_guess = curve_list[j][5](t_) 
						v2dR_exact = scipy.optimize.root(ls_cconstraint_theta,v2d_guess,(na,nb,t_,icase))
						roots.append((v2dR_exact[4],t_,depl_ERcurve1(t_,ER_curves_ex[j][0])))
					if len(roots)>0:
						intersecting.add_edge(j,nnj)
						intersecting.edges[(nnj,j)]=roots
	# Determine the ground state curve segments, starting at the smallest possible cylinder radius
	# (i) Find the curve(s) that attain the smallest cylinder radii. 
	# (ii) Starting from these, find the next intersection (e.g. by ordering edges by the min and max of their associated lists of roots)
	# (iii) Define the gound state curve by a list of intervals (between adjacent roots) together with the curve index(ices) associated with
	#	each. Note: if just one curve index is provided, the other can be found from curve_list[curve_index][16].
	print "Determining ground state"
	ground_states = []
	first_curve = 0
	min_R = 1e99
	for j in range(len(curve_list)):
		if curve_list[j][4]==0:
			R_cand = curve_list[j][6]
		else:
			R_cand = min(curve_list[j][10],curve_list[j][11])
		if min_R>R_cand:
			first_curve=j
			min_R = R_cand
	# Find the next intersection on the first branch
	R0 = min_R
	curve0 = first_curve
	ground_states.append([R0,curve0])
	# Indices: first_curve, and first_curve_p (which may equal first_curve)
	exhausted = False
	while not exhausted:
		# Find the next intersection on current branch
		R1 = 1e99
		exhausted = True
		for nnj in intersecting.nodes[curve0].nbors:
			for root in intersecting.edges[(min(nnj,curve0),max(nnj,curve0))]:
				if (root>R0) and (root<R1):
					R1=root
					curve1 = nnj
					exhausted=False
		if not exhausted:
			ground_states.append([R1,curve1])
			R0 = R1
			curve0=curve1
	# Add the terminal radius of the last curve:
	last_curve = ground_states[-1][1]
	last_radius=max(curve_list[last_curve][10:14])
	ground_states.append([last_radius,last_curve])
	if printout!=None:
		# Display the ground state curve over each interval up to Rmax
		for j in range(len(ground_states)-1):
			print "Displaying "+nth(j)+" segment of ground state"
			curve_j = ground_states[j][1]
			[na,nb,ocase,icase]=curve_list[curve_j][:4]
			pars = [na,nb,ocase,icase,depl_str,dpr,V_min,curve_j]
			R0 = ground_states[j][0]
			R1 = ground_states[j+1][0]
			t0 = ground_states[j][1]
			t1 = ground_states[j+1][1]
			Rsamples = np.linspace(R0,R1,50)
			Esamples = np.zeros(50)
			if curve_list[curve_j][4]==0:
				for jj in range(50):
					Esamples[jj]=depl_ERcurve0(Rsamples[jj],pars)
			else:
				for jj in range(50):
					Esamples[jj]=depl_ERcurve1(tsamples[jj],pars)
			plt.plot(Rsamples,Esamples,s=0.5)
		plt.savefig(printout[0],dpi=printout[1])
		plt.show()
		plt.clf()
	return ground_states

compute_dt_intersections(printout=["test_dt_landscape.png",1028])

# Note: the orientation energy also depends on the density of points along the z axis
# Make a plot of the line slip tension as a function of radius for each curve:
def compute_ls_tensions(depl_str,dpr,Nsamples):
	# Create directory if it does not already exist:
	mkdir_s("LS_TENSION")
	rthrsq = (2+2*dpr)
	rthrsq*=rthrsq
	V_min = V_depl_min(depl_str,dpr)
	for j in range(len(curve_list)):
		plot_data = np.zeros([2,Nsamples])
		# Universal curve parameters:
		[na,nb,ocase,icase,branched]=curve_list[j][:5]
		Esamples = np.zeros(Nsamples)
		Rsamples = np.zeros(Nsamples)
		print "Computing for curve "+str(j)
		if curve_list[j][4]==0:
			# curve parametrized by radius.
			R0 = curve_list[j][6]
			R1 = curve_list[j][7]
			Rsamples = np.linspace(R0,R1,Nsamples) # This may require adjustment for extrapolation
			for jj in range(Nsamples):
				# Compute the line slip vector, and the associated energy
				v2d_guess = curve_list[j][5](Rsamples[jj]) # initial guess
				v2dr = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,Rsamples[jj],[ocase,icase]))
				if v2dr.success:
					v2d = v2dr.x
					# Energy depends on either vls+va or vls-va (whichever is shorter)
					v02d = v2d[:2]+v2d[4:]
					v12d = v2d[4:]-v2d[:2]
					v0 = x3d(v02d,Rsamples[jj])
					v1 = x3d(v12d,Rsamples[jj])
					d0sq = np.dot(v0,v0)
					d1sq = np.dot(v1,v1)
					Esamples[jj] = -V_min
					if d1sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d1sq),depl_str,dpr)
					if d0sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d0sq),depl_str,dpr)
					Esamples[jj]/=abs(v2d[1])
		elif curve_list[j][4]==1:
			# parametrized by theta: solve for both R(theta) and E(theta)
			theta0 = curve_list[j][6]
			theta1 = curve_list[j][7]
			tsamples = np.linspace(theta0,theta1,Nsamples)
			for jj in range(Nsamples):
				v2dR_guess = curve_list[j][5](tsamples[jj])
				v2dRr = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,tsamples[jj],icase))
				if v2dRr.success:
					v2dR = v2dRr.x
					Rsamples[jj] = np.copy(v2dR[4])
					v02d = v2dR[:2]+v2dR[2:4]
					v12d = v2dR[2:4]-v2dR[:2]
					v0 = x3d(v02d,Rsamples[jj])
					v1 = x3d(v12d,Rsamples[jj])
					d0sq = np.dot(v0,v0)
					d1sq = np.dot(v1,v1)
					Esamples[jj] = -V_min
					if d1sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d1sq),depl_str,dpr)
					if d0sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d0sq),depl_str,dpr)
					Esamples[jj]/=abs(v2dR[1])
		plot_data[0,:] = np.copy(Rsamples)
		plot_data[1,:] = np.copy(Esamples)
		np.savetxt("LS_TENSION/E_curve_"+str(j)+".dat",plot_data[:,:])
	#plt.savefig("ls_energies.png")
	#plt.clf()

def plot_ls_tensions():
	print "Generating line slip energy plot"
	for j in range(len(curve_list)):
		dfile = "LS_TENSION/E_curve_"+str(j)+".dat"
		Edat = np.loadtxt(dfile)
		plt.plot(Edat[0,:],Edat[1,:],linewidth=0.1)
	plt.xlabel("Radius")
	plt.ylabel("l.s. tension per unit height")
	plt.savefig("ls_tensions.png",dpi=1024)
	plt.clf()

def E_packing_bnd(R,V_min,V_min_cyl):
	return (np.pi*R)*(3*V_min+V_min_cyl)/(np.sqrt(3))

def plot_cr_tensions(depl_str,dpr,image_name="cr_tensions_lb.png",CR_INDEX="",subtract_lb=False,text=True,curve_emphasis=[],gs_bounds=None):
	V_min = V_depl_min(depl_str,dpr)
	cases = {}
	print "Generating plot of full crystal tension for curve list of length "+str(len(curve_list))
	for j in range(len(curve_list)):
		[na,nb,ocase,icase]=curve_list[j][:4]
		dfile = "CR_TENSION"+str(CR_INDEX)+"/RTE_curve"+str(j)+".dat"
		RTEdat=np.loadtxt(dfile)
		E_lb_samples = np.zeros(len(RTEdat))
		if subtract_lb:
			for jj in range(len(RTEdat)):
				V_min_cyl = V_depl_cyl_min(depl_str,dpr,RTEdat[jj,0]-1)
				E_lb_samples[jj]=E_packing_bnd(RTEdat[jj,0],V_min,V_min_cyl)
		if j not in curve_emphasis:
			plt.plot(RTEdat[:,0],RTEdat[:,2]-E_lb_samples[:],'k',linewidth=0.3)
		else:
			plt.plot(RTEdat[:,0],RTEdat[:,2]-E_lb_samples[:],'m',linewidth=1.5)
	# Plot the energies at points where fractional vacancies have equivalent geometries:
	E_R_equiv = {}
	for j in range(len(curve_list)):
		if j in R_equiv:
			[na,nb,ocase,icase]=curve_list[j][:4]
			case = [ocase,icase]
			if curve_list[j][4]==0:
				V_min_cyl = V_depl_cyl_min(depl_str,dpr,R_equiv[j]-1)
				# Compute the curve energy at radius R: get_domain_tension(v2d,nb,R,depl_str,dpr,V_min,V_min_cyl)
				v2d_R_equiv = scipy.optimize.root(ls_cconstraint,curve_list[j][5](R_equiv[j]),(na,nb,R_equiv[j],case)).x
				E_val = get_domain_tension(v2d_R_equiv,[na,nb],R_equiv[j],depl_str,dpr,V_min,V_min_cyl)
				if subtract_lb:
					E_val-=E_packing_bnd(R_equiv[j],V_min,V_min_cyl)
				R_val = R_equiv[j]
				E_R_equiv[j]=[R_val,E_val]
			elif curve_list[j][4]==1:
				V_min_cyl = V_depl_cyl_min(depl_str,dpr,R_equiv[j][0]-1)
				v2dR_R_equiv = scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](R_equiv[j][1]),(na,nb,R_equiv[j][1],icase)).x
				v2d_R_equiv = get_v2d_of_v2dR(v2dR_R_equiv,R_equiv[j][1],[ocase,icase])
				E_val = get_domain_tension(v2d_R_equiv,[na,nb],R_equiv[j][0],depl_str,dpr,V_min,V_min_cyl)
				if subtract_lb:
					E_val-=E_packing_bnd(R_equiv[j][0],V_min,V_min_cyl)
				R_val = R_equiv[j][0]
				E_R_equiv[j]=[R_val,E_val]
	E_R_vals = np.array(E_R_equiv.values())
	plt.scatter(E_R_vals[:,0],E_R_vals[:,1],s=0.2,marker='.',edgecolors='none')
	# Include the ground state (if flagged)
	if gs_bounds!=None:
		gs_Rsamples = np.linspace(gs_bounds[0],gs_bounds[1],int(500*(gs_bounds[1]-gs_bounds[0]))+1)
		gs_Esamples = np.zeros(len(gs_Rsamples))
		for j in range(len(gs_Rsamples)):
			V_min_cyl = V_depl_cyl_min(depl_str,dpr,gs_Rsamples[j])
			print "Processing sample "+str(j)
			ex_list = get_closed_crystals(gs_Rsamples[j],curve_list)
			min_E = 1e99
			min_jj = -1
			# Expect two nearly identical minima, with transposed commensurate indices
			for jj in range(len(ex_list)):
				v2d=ex_list[jj][0]
				#get_domain_tension(v2d,ns,R,depl_str,dpr,V_min,V_min_cyl,tol=1e-5)		
				E_jj = get_domain_tension(v2d,ex_list[jj][1:3],gs_Rsamples[j],depl_str,dpr,V_min,V_min_cyl)
				if E_jj<min_E:
					min_E=E_jj
					min_jj=jj
			gs_Esamples[j]=min_E
			if subtract_lb:
				gs_Esamples[j]-=E_packing_bnd(gs_Rsamples[j],V_min,V_min_cyl)
		plt.plot(gs_Rsamples,gs_Esamples,'r',linewidth=0.5)
		plt.axis([gs_bounds[0]-0.01,gs_bounds[1]+0.01,-0.001,0.07])
	if text:
		for j in range(len(curve_list)):
			if j in R_equiv:
				[na,nb,ocase,icase]=curve_list[j][:4]
				if curve_list[j][4]==1:
					if (na,nb) not in cases:
						plt.text(E_R_equiv[j][0],E_R_equiv[j][1],str(j),fontsize=0.25,color='b')
						#plt.text(RTEdat[-1,0],RTEdat[-1,2],str([na,nb]),fontsize=1,color='b')
						cases[(na,nb)]=True
				else:
					if (na,nb) not in cases:
						plt.text(E_R_equiv[j][0],E_R_equiv[j][1],str(j),fontsize=0.25,color='r')
						#plt.text(RTEdat[-1,0],RTEdat[-1,2],str([na,nb]),fontsize=1,color='r')
						cases[(na,nb)]=True
	# Include a lower bound assuming h.c.p. volume filling:
	if not subtract_lb:
		E_lower_bnd = np.zeros([100,2])
		for i in range(100):
			R_i = 8.*i/100
			V_min = V_depl_min(depl_str,dpr)
			V_min_cyl = V_depl_cyl_min(depl_str,dpr,R_i)
			E_lower_bnd[i,0] = R_i
			E_lower_bnd[i,1] = E_packing_bnd(R_i,V_min,V_min_cyl)
		plt.plot(E_lower_bnd[:,0],E_lower_bnd[:,1],linewidth=0.15)
	plt.xlabel("Cylinder radius (units of $r_{core}$)")
	if subtract_lb:
		plt.ylabel("energy/unit length minus lin. lower bnd. (units of $n_Sk_BTr_{core}^2$)")
	else:
		plt.ylabel("energy per unit length (units of $n_Sk_BTr_{core}^2$)")
	plt.savefig(image_name,dpi=1024)
	plt.clf()

def ref_and_rev(ind):
	rval1 = (ind[2],ind[3],ind[0],ind[1])
	rval2 = (ind[1],ind[0],ind[3],ind[2])
	rval3 = (ind[3],ind[2],ind[1],ind[0])
	rval = [tuple(ind),rval1,rval2,rval3]
	return rval


# CHECK THIS
# Note on usage: each (integer valued) mode corresponds to a different parametrization
# 	Mode: 	Type:
#	0	curve1(R), curve2(R)
#	1 	curve1(R), curve2(theta)
#	2 	curve1(theta), curve2(theta)
def check_isect(curve1,curve2,Ecurve1,Ecurve2,Nsamples,pars,mode=0):
	[depl_str,dpr,V_min] = pars
	isect = False # Return
	[na1,nb1] = curve1[:2]
	case1 = curve1[2:4]
	[na2,nb2] = curve2[:2]
	case2 = curve2[2:4]
	Esamples = -1e99*np.ones(Nsamples)
	signs = np.zeros(Nsamples,dtype=int)
	if mode==0:
		Rmin1 = min(curve1[6:8])
		Rmax1 = max(curve1[6:8])
		Rmin2 = min(curve2[6:8])
		Rmax2 = max(curve2[6:8])
		# Globla min and max:
		Rmin = max(Rmin1,Rmin2)
		Rmax = min(Rmax1,Rmax2)
		# both parametrized by R
		Rsamples = np.linspace(Rmin,Rmax,Nsamples)
		# Look for sign changes in ERcurves[j][0]-ERcurves[jj][0]
		for ell in range(Nsamples):
			Esamples[ell] = Ecurve1[0](Rsamples[ell],na1,nb1,case1)-Ecurve2[0](Rsamples[ell],na2,nb2,case2)
			signs[ell] = (Esamples[ell]>-1e50)*np.sign(Esamples[ell])
		for ell in range(Nsamples-1):
			if signs[ell]!=signs[ell+1]:
				isect = True
				break
	elif mode==1:
		# curve 1 parametrized by R, curve 2 by theta
		Rmin1 = min(curve1[6:8])
		Rmax1 = max(curve1[6:8])
		tmin2 = min(curve2[6:8])
		tmax2 = max(curve2[6:8])
		tsamples = np.linspace(tmin2,tmax2,Nsamples)
		# Here we have E1(R), E2(theta), R2(theta), and we wish to solve E1(R)==E2(theta), R==R2(theta)
		# Essentially, we define a 1D function E1(R2(theta))==E2(theta) and look for its roots with bracketing
		# Auxiliary function
		aux_deltaE = lambda v2dR,theta: Ecurve1[0](v2dR[4],na1,nb1,case1)-get_domain_tension_theta(v2dR,nb,theta,depl_str,dpr,V_min,V_depl_cyl_min(depl_str,dpr,v2dR[4]-1))
		for ell in range(Nsamples):
#			v2dR = Ecurve2[0](tsamples[ell],na2,nb2,case2[1])
			v2dR = Ecurve2[0](tsamples[ell],na2,nb2,case2[1])
			if v2dR[4]>Rmin and v2dR[4]<Rmax:
				Esamples[ell] = aux_deltaE(v2dR,tsamples[ell])
				signs[ell] = (Esamples[ell]>-1e50)*np.sign(Esamples[ell])
		# Check for sign changes of deltaE on the relevant domain
		for ell in range(Nsamples-1):
			if signs[ell]>0 and signs[ell+1]<0 or signs[ell]<0 and signs[ell+1]>0:
				isect = True
				break
		return isect
	elif mode==2:
		# both parametrized by theta
		R11 = curve1[10]
		R21 = curve1[11]
		R2vm1 = curve1[13]
		R12 = curve2[10]
		R22 = curve2[11]
		R2vm2 = curve2[13]
		t11 = curve1[6]
		t21 = curve1[7]
		t2v1 = curve1[14]
		t12 = curve2[6]
		t22 = curve2[7]
		t2v2 = curve2[14]
		t1samples = np.linspace(t11,t21,Nsamples)
		# Look for existence of solutions to R1 = R2, E1 = E2
#		Rdiff = lambda R1,t2: R1-Ecurve2[0](t2,na2,nb2,case2[1])[4]
		Rdiff = lambda R1,t2: R1-Ecurve2[0](t2)[4]
		# Solve for theta2(theta1) through R1(theta1)=R2(theta2)
		# 	R1 and R2 are not necessarily monotonic. 
		# First, for each theta1, find points theta2a,theta2b where Rdiff(theta1,theta2) changes sign
		theta2_p = lambda R1: scipy.optimize.brentq(lambda x: Rdiff(R1,x),t21,t2v2)
		theta2_m = lambda R1: scipy.optimize.brentq(lambda x: Rdiff(R1,x),t2v2,t22)
		Esamples_m = -1e99*np.ones(Nsamples)
		signs_m = np.zeros(Nsamples,dtype=int)
		# Check each possible branch to see where the curves might intersect
		for ell in range(Nsamples):
			v2dR1 = Ecurve1[0](t1samples[ell],na1,nb1,case1[1])
			V_min_cyl1 = V_depl_cyl_min(depl_str,dpr,v2dR1[4]-1)
			if v2dR1[4]>R12 and v2dR1[4]<R2vm2:
				t2 = theta2_p(t1samples[ell])
#				v2dR2 = Ecurve2[0](t2,na2,nb2,case2[1])
				v2dR2 = Ecurve2[0](t2,na2,nb2,case2[1])
				Esamples[ell] = get_domain_tension_theta(v2dR1,nb,theta,V_min,V_min_cyl1)-get_domain_tension_theta(v2dR2,nb,t2,V_min,V_depl_cyl_min(depl_str,dpr,v2dR2[4]-1))
				signs[ell] = np.sign(Esamples[ell])
			if v2dR1[4]>R22 and v2dR1[4]<R2vm2:
				t2 = theta2_m(t1samples[ell])
#				v2dR2 = Ecurve2[0](t2,na2,nb2,case2[1])
				v2dR2 = Ecurve2[0](t2,na2,nb2,case2[1])
				V_min_cyl2 = V_depl_cyl_min(depl_str,dpr,v2dR2[4]-1)
				Esamples_m[ell] = get_domain_tension_theta(v2dR1,nb,theta,V_min,)-get_domain_tension_theta(v2dR2,nb,t2,V_min,V_min_cyl2)
				signs_m[ell] = np.sign(Esamples_m[ell])
		# Check for sign changes in each branch:
		for ell in range(Nsamples-1):
			if signs[ell]>0 and signs[ell+1]<0 or signs[ell]<0 and signs[ell+1]>0:
				isect = True
				break
			if signs_m[ell]>0 and signs_m[ell+1]<0 or signs_m[ell]<0 and signs_m[ell+1]>0:
				isect = True
				break
	return isect

# CHECK THIS
# update_curve_prox(curve_prox,curve_list,ERcurves,j)
def update_curve_prox(curve_prox,curve_list,ERcurves,j,pars,Nsamples=10):
	depl_str = pars[0]
	dpr = pars[1]
	V_min = pars[2]
	if j in curve_prox:
		pass
	else:
		curve_prox[j] = []
	[naj,nbj] = curve_list[j][:2]
	casej = curve_list[j][2:4]
	if curve_list[j][4]==0:
		Rmin = min(curve_list[j][6:8])
		Rmax = max(curve_list[j][6:8])
		for jj in range(j+1,len(curve_list)):
			# Check if the R-ranges overlap
			Rminjj = min(curve_list[jj][10:12])
			Rmaxjj = max(curve_list[jj][10],curve_list[jj][11],curve_list[jj][13])
			if Rminjj>Rmax or Rmin > Rmaxjj:
				continue
			else:
				[najj,nbjj] = curve_list[jj][:2]
				casejj = curve_list[jj][2:4]
				if curve_list[jj][4]==0:
					isect = check_isect(curve_list[j],curve_list[jj],ERcurves[j],ERcurves[jj],Nsamples,pars,0)
					if isect:
						curve_prox[j].append(jj)
						if jj in curve_prox:
							curve_prox[jj].append(j)
						else:
							curve_prox[jj] = [j]
				elif curve_list[jj][4]==1:
					isect = check_isect(curve_list[j],curve_list[jj],ERcurves[j],ERcurves[jj],Nsamples,pars,1)
					if isect:
						curve_prox[j].append(jj)
						if jj in curve_prox:
							curve_prox[jj].append(j)
						else:
							curve_prox[jj]=[j]
	elif curve_list[j][4]==1:
		Rmin = min(curve_list[j][10:12])
		Rmax = max(curve_list[j][10],curve_list[j][11],curve_list[j][12])
		for jj in range(j+1,len(curve_list)):
			Rminjj = min(curve_list[jj][10:12])
			Rmaxjj = max(curve_list[jj][10],curve_list[jj][11],curve_list[jj][13])
			if Rminjj>Rmax or Rmin > Rmaxjj:
				continue
			else:
				if curve_list[jj][4]==0:
					isect = check_isect(curve_list[jj],curve_list[j],ERcurves[jj],ERcurves[j],Nsamples,pars,1)
					if isect:
						curve_prox[j].append(jj)
						if jj in curve_prox:
							curve_prox[jj].append(j)
						else:
							curve_prox[jj]=[j]
				else:
					isect = check_isect(curve_list[j],curve_list[jj],ERcurves[j],ERcurves[jj],Nsamples,pars,2)
					if isect:
						curve_prox[j].append(jj)
						if jj in curve_prox:
							curve_prox[jj].append(j)
						else:
							curve_prox[jj]=[j]
	return


# CHECK THIS
def polish_roots(curve_prox,curve_list,ERcurves,pars,j,jj):
	[naj,nbj,ocasej,icasej] = curve_list[j][:4]	
	[najj,nbjj,ocasejj,icasejj] = curve_list[jj][:4]
	if (j,jj) in curve_prox:
		if curve_list[j][4]==0 and curve_list[jj][4]==0:
			Rminj = min(curve_list[j][6:8])
			Rmaxj = max(curve_list[j][6:8])
			# Enumerate intersections between curves j and jj:
			Rminjj = min(curve_list[jj][6:8])
			Rmaxjj = max(curve_list[jj][6:8])
			Rmin = max(Rminj,Rminjj)
			Rmax = min(Rmaxj,Rmaxjj)
			roots = custom_bracket(lambda R: ERcurves[j][0](R,naj,nbj,[ocasej,icasej])-ERcurves[j][0](R,najj,nbjj,[ocasejj,icasejj]),Rmin,Rmax)
		elif curve_list[j][4]==0 and curve_list[jj][4]==1:
			# Define similar auxiliary functions of theta2:
			tminjj = min(curve_list[jj][6:8])
			tmaxjj = max(curve_list[jj][6:8])
			[t1jj,t2jj]= curve_list[jj][6:8]
			t2vmjj = curve_list[jj][14]
			# Check for endpoints in theta2:
			[R1jj,R2jj] = curve_list[jj][10:12]
			R2vmjj = curve_list[jj][13]
			Rminjj = min(R1jj,R2jj)
			Rmaxjj = max(R1jj,R2jj,R2vmjj)
			# Refinement: maintain a list of intervals over which to parametrize Rj by thetajj 
			theta_root = -1e99*np.ones(4)
			if R1jj<Rminj:
				# Upper part of curve overshadows [Rminj,(some endpoint1)]
				# Use the interval R1jj,R2vmjj to determine where Rjj(theta)==Rminj (upper branch)
#				Rjj = lambda theta: ERcurves[jj][0](theta,najj,nbjj,icasejj)[4]-Rminj
				Rjj = lambda theta: ERcurves[jj][0](theta)[4]-Rminj
				theta_root[0] = scipy.optimize.brentq(Rjj,t1jj,t2vmjj)
			elif R1jj<Rmaxj:
				theta_root[0] = t1jj
			if R2jj<Rminj:
				# Lower part of curve overshadows [Rminj,(some endpoint2)]
				# Same, with interval R2jj,R2vmjj to determine where Rjj(theta)==Rminj (lower branch)
#				Rjj = lambda theta: ERcurves[jj][0](theta,najj,nbjj,icasejj)[4]-Rminj
				Rjj = lambda theta: ERcurves[jj][0](theta)[4]-Rminj
				theta_root[1] = scipy.optimize.brentq(Rjj,t2jj,t2vmjj)
			elif R2jj<Rmaxj:
				theta_root[1] = t2jj
			if Rmaxjj>Rmaxj:
				# Solve for values of theta where Rjj(theta)=Rmaxj
				if R1jj<Rmaxj:
					# Upper part of curve overshadows [(some initial point),Rmaxj]
#					Rjj = lambda theta: ERcurves[jj][0](theta,najj,nbjj,icasejj)[4]-Rmaxj
					Rjj = lambda theta: ERcurves[jj][0](theta)[4]-Rmaxj
					theta_root[2] = scipy.optimize.brentq(Rjj,t1jj,t2vmjj)
					# Note: this case implies theta_root[0] is defined
					# 	If this condition is not met, the upper curve is not relevant
				if R2jj<Rmaxj:
					# Lower part of curve overshadows [(some initial point),Rmaxj]
#					Rjj = lambda theta: ERcurves[jj][0](theta,najj,nbjj,icasejj)[4]-Rmaxj
					Rjj = lambda theta: ERcurves[jj][0](theta)[4]-Rmaxj
					theta_root[3] = scipy.optimize.brentq(Rjj,t2jj,t2vmjj)
					# Note: this case implies theta_root[1] is defined
					#	Condition not met implies the lower branch is not relevant
			else:
				# Domain consists of a single interval (use theta_root[0] and theta_root[2] to describe the endpoints)
				theta_root[2] = np.copy(theta_root[1])
				theta_root[3] = theta_root[1] = -1e99
			# Note: there should be at most two intervals: [theta_root[0],theta_root[2]] and possibly [theta_root[1],theta_root[3]]	
			theta_ivals = []
			if theta_root[0]>-1e50 and theta_root[2]>-1e50:
				theta_ivals.append([theta_root[0],theta_root[2]])
			if theta_root[1]>-1e50 and theta_root[3]>-1e50:
				theta_ivals.append([theta_root[1],theta_root[3]])
			# Define functions to compute the deltaE(theta) root(s)
			R_of_theta = lambda theta: ERcurves[jj][0](theta,najj,nbjj,icasejj)[4] # CHECK SYNTAX
			# Use theta_root for endpoints in a root-finding algorithm for the difference of energies. (v2dR,nb,theta,V_min)
			dE = lambda theta: get_domain_tension_theta(ERcurves[jj][0](theta,najj,nbjj,icasejj),nbjj,theta,depl_str,dpr,V_min)-ERcurves[j][0](R_of_theta(theta),naj,nbj,[ocasej,icasej]) # CHECK SYNTAX
			# Find roots of dE within each interval:
			if len(theta_ivals)==1:
				roots = custom_bracket(dE,theta_ivals[0][0],theta_ivals[0][1])
			elif len(theta_ivals)==2:
				roots1 = custom_bracket(dE,theta_ivals[0][0],theta_ivals[0][1])
				roots2 = custom_bracket(dE,theta_ivals[1][0],theta_ivals[1][1])
				roots = roots1+roots2
		elif curve_list[j][4]==1 and curve_list[jj][4]==1:
			[t1j,t2j] = curve_list[j][6:8]
			[t1jj,t2jj] = curve_list[jj][6:8]
			# Solve for thetajj in terms of thetaj
			v2dRjj = lambda theta: scipy.optimize.root(ls_cconstraint_theta,curve_list[jj][5](theta),(najj,nbjj,theta,icasejj)).x # Syntax:
			v2dRj = lambda theta: scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](theta),(naj,nbj,theta,icasej)).x # Syntax
			RjmRjj = lambda thetaj,thetajj: v2dRj(thetaj)[4]-v2dRjj(thetajj)[4]
			EjmEjj = lambda thetaj,thetajj: get_domain_tension(v2dRj(thetaj),nbj,depl_str,dpr,V_min,V_depl_cyl_min(depl_str,dpr,v2dRj(thetaj)[4]-1))-get_domain_tension(v2dRjj(thetajj),nbjj,depl_str,dpr,V_min,V_depl_cyl_min(depl_str,dpr,v2dRjj(thetajj)[4]-1)) #CHECK SYNTAX
			vf = lambda thetas: np.array([RjmRjj(thetas[0],thetas[1]),EjmEjj(thetas[0],thetas[1])])
			roots = custom_2d_solver(vf,[t1j,t2j],[t1jj,t2jj])
		else:
			print "Error in polish_roots(curve_prox,curve_list,ERcurves,pars,j,jj):\nUsage:\n The index of the curve parametrized by R should always come first (i.e. be passed in the 'j' slot.)" 
			quit()
	else:
		print "Error in polish_roots(curve_prox,curve_list,ERcurves,pars,j,jj):\n curves "+str([j,jj])+" are not proximal"
		quit()
	return roots
		
def refine_intersections(curve_prox,curve_list,ERcurves,pars):
	isect_pts = {} # Structure to be returned
	[depl_str,dpr,V_min] = pars
	for j in range(len(curve_list)):
		if curve_list[j][4]==0:
		# Parametrized by R
			Rminj = min(curve_list[j][6:8])
			Rmaxj = max(curve_list[j][6:8])
			[naj,nbj,ocasej,icasej] = curve_list[j][:4]
			try:
				for jj in curve_prox[j]:
					if (jj,j) in isect_pts:
						continue
					else:
						[najj,nbjj,ocasejj,icasejj]= curve_list[j][:4]
						isect_pts[(j,jj)] = polish_roots(curve_prox,curve_list,ERcurves,pars,j,jj)
			except:
				print "Warning: curve_prox does not contain curve "+str(j)
				pass
	return isect_pts

def compute_gs_tensions(depl_str,dpr,Rmin,Rmax,Nsamples,tol=1e-6):
	V_min = V_depl_min(depl_str,dpr)
	Rs = np.linspace(Rmin,Rmax,Nsamples)
	Egs = np.zeros(Nsamples)
	tpts = []
	tRs = []
	tEgs = []
	tcurves = []
	for j in range(Nsamples):
		V_min_cyl = V_depl_cyl_min(depl_str,dpr,Rs[j])
		ex_list = get_closed_crystals(Rs[j],curve_list)
		Egs[j] = 1e99
#		pcurves = [(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1)]
		min_jj = -1
		for jj in range(len(ex_list)):
			v2d = ex_list[jj][0]
			nb = ex_list[jj][2]
			gamma_jj = get_domain_tension(v2d,nb,Rs[j],depl_str,dpr,V_min,V_min_cyl)
			if gamma_jj<Egs[j]:
				Egs[j] = gamma_jj
				min_jj = jj
	valid = np.where(Egs<1e50)
	print "Generating full energy plot"
	for j in range(len(curve_list)):
		dfile = "LS_TENSION/E_curve_"+str(j)+".dat"
		Edat = np.loadtxt(dfile)
		plt.plot(Edat[0,:],Edat[1,:],'k',linewidth=0.5)
	plt.plot(Rs[valid],Egs[valid],'c',linewidth=.7)
	plt.xlabel("Radius")
	plt.title("Ground state line slip tension")
	plt.xlabel("Radius")
	plt.ylabel("$\lambda$ crystal")
	plt.savefig("cr_tensions_gs.png",dpi=1024)
	
def compute_full_gs_tensions(curve_list,depl_str,dpr,Rmin,Rmax,Nsamples,tol=1e-6):
	V_min = V_depl_min(depl_str,dpr)
	Rs =np.linspace(Rmin,Rmax,Nsamples) 
	Egs = np.zeros(Nsamples)
	tpts = []
	tRs = []
	tEgs = []
	tcurves = []
#	current_curves = [(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1)] # Expect two curves with comparable energy, when n0!=n5
	for j in range(Nsamples):
		V_min_cyl = V_depl_cyl_min(depl_str,dpr,Rs[j])
		ex_list = get_closed_crystals(Rs[j],curve_list)
		Egs[j] = 1e99
#		pcurves = [(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1),(-1,-1,-1,-1)]
		min_jj = -1
		for jj in range(len(ex_list)):
			v2d = ex_list[jj][0]
			[na,nb,ocase,icase]=ex_list[jj][1:5]
			gamma_jj = get_domain_tension(v2d,[na,nb],Rs[j],depl_str,dpr,V_min,V_min_cyl)
			if gamma_jj<Egs[j]:
				Egs[j] = gamma_jj
				min_jj = jj
	valid = np.where(Egs<1e50)
	tpts = []
	tRs = []
	tEgs = []
	tcurves = []
	ERcurves = [] # A list of lists of curves ('true' curves and auxiliary curves)
	N_ecurves = [] # Number of curves in each entry of Ecurves
	print "Generating full energy plot"
	Nset = 0
	for j in range(len(curve_list)):
		[na,nb,ocase,icase] = curve_list[j][:4]
		dfile = "CR_TENSION/RTE_curve"+str(j)+".dat"
		RTEdat = np.loadtxt(dfile)
		plt.plot(RTEdat[0,:],RTEdat[2,:],'k',linewidth=0.5)
		ERcurves.append([])
		N_ecurves.append(0)
		# NOTE: Use these curves to define a R^2->R^2 function for each pair of possibly relevant curves
		# 	(next step: first check for intersection, then once likely, look for a solution)
		if curve_list[j][4]==0:
			# Curve already parametrized by 'R'
			E_of_R = lambda R,na,nb,case: get_domain_tension(scipy.optimize.root(ls_cconstraint,curve_list[j][5](R),(na,nb,R,case)).x,nb,R,depl_str,dpr,V_min,V_depl_cyl_min(depl_str,dpr,R-1))
			ERcurves[j].append(E_of_R) 
			N_ecurves[j]=1
		elif curve_list[j][4]==1:
			if Nset==0:
				print "Defining first v2dR_of_theta curve"
			Nset+=1
			# Auxiliary curves needed: must parametrize both R(theta) and E(theta)Major Moment
			# Check that curve_list[j][5] does indeed return a 5 dimensional vector:
			test_vector = curve_list[j][5](0.5*(curve_list[j][6]+curve_list[j][7]))
			if len(test_vector)!=5:
				print "Error in compute_full_gs_tensions: curve parametrized by theta returns a vector with dimension "+str(len(test_vector))
			v2dR_of_theta = lambda theta,na,nb,icase: scipy.optimize.root(ls_cconstraint_theta,curve_list[j][5](theta),(na,nb,theta,icase)).x
			ERcurves[j].append(v2dR_of_theta)
			if Nset==1:
				print "(Done)"
	print "Number of v2dR_of_theta curves set: "+str(Nset)
	# Determine which curves are proximal (in [E,R] space)
	curve_prox = {}
	for j in range(len(curve_list)):
		update_curve_prox(curve_prox,curve_list,ERcurves,j,[depl_str,dpr,V_min])
	# Compute intersections between nearby curves
	isect_pts = refine_intersections(curve_prox,curve_list,ERcurves,[depl_str,dpr,V_min])
	# Note on usage: isect_pts is a dictionary that gives a list of points in (E,R) space for every pair of intersecting curves
	# 		that is, if curves (i,j) intersect at points (E1,R1), (E2,R2), ..., (EN,RN), then isect_pts[(i,j)]=[(E1,R1),...,(EN,RN)] 
	#for j in range(len(curve_list)):
	#	dfile = "LS_TENSION/E_curve_"+str(j)+".dat"
	#	Edat = np.loadtxt(dfile)
		#plt.text(Edat[0,len(Edat[0,:])/3],Edat[1,len(Edat[1,:])/3],str(curve_list[j][15]),fontsize=2)
	plt.plot(Rs[valid],Egs[valid],'c',linewidth=.7)
	plt.scatter(tRs,tEgs,s=2,c="y",edgecolors='none')
	plt.xlabel("Radius")
	plt.title("Ground state line slip tension")
	plt.xlabel("Radius")
	plt.ylabel("$\gamma$")
	plt.savefig("crystal_tensions_gs.png",dpi=1024)
	plt.clf()


def plot_theta_gs_R(R_min,R_max,depl_str=1.0,dpr=0.05,Nsamples=None,Ncsamples=100,postfix="",curve_emphasis=[]):
	if Nsamples==None:
		Nsamples = 200*(int(R_max-R_min)+1)
	Rsamples_gs = np.linspace(R_min,R_max,Nsamples)
	tsamples_gs = np.zeros(Nsamples)
	Rtsamples_gs = []
	V_min = V_depl_min(depl_str,dpr)
	print "Computing ground state energies in interval "+str([R_min,R_max])
	for j in range(len(Rsamples_gs)):
		print str(j)+" of "+str(len(Rsamples_gs))
		V_min_cyl = V_depl_cyl_min(depl_str,dpr,Rsamples_gs[j])
		ex_list=get_closed_crystals(Rsamples_gs[j],curve_list)
		min_E = 1e99
		min_jj = -1
		# Expect two nearly identical minima, with transposed commensurate indices
		for jj in range(len(ex_list)):
			v2d=ex_list[jj][0]
			#get_domain_tension(v2d,ns,R,depl_str,dpr,V_min,V_min_cyl,tol=1e-5)		
			E_jj = get_domain_tension(v2d,ex_list[jj][1:3],Rsamples_gs[j],depl_str,dpr,V_min,V_min_cyl)
			if E_jj<min_E:
				min_E=E_jj
				min_jj=jj
		[n0,n5,n0p,n5p]=ex_list[min_jj][-1]
		# Look for the dual orientation (with reversed n0,n5)
		dual_jj = -1
		for jj in range(len(ex_list)):
			if jj==min_jj:
				continue
			[n0jj,n5jj,n0jjp,n5jjp]=ex_list[jj][-1]
			if ((n0jj==n5)and (n5jj==n0)) and ((n0jjp==n5p) and (n5jjp==n0p)):
				dual_jj=jj
		# Define theta(s):
		v0 = get_v0_of_vabls(ex_list[min_jj][0],ex_list[min_jj][3],ex_list[min_jj][4])
		tsamples_gs[j]=np.arcsin(v0[1]*0.5)
		Rtsamples_gs.append([Rsamples_gs[j],tsamples_gs[j]])
		if dual_jj!=-1:
			v0p = get_v0_of_vabls(ex_list[dual_jj][0],ex_list[dual_jj][3],ex_list[dual_jj][4])
			thetap = np.arcsin(v0p[1]*0.5)
			Rtsamples_gs.append([Rsamples_gs[j],thetap])
	Rtsamples_gs = np.array(Rtsamples_gs)
	print "Displaying curves"
	for j in range(len(curve_list)):
		[na,nb,ocase,icase]=curve_list[j][:4]
		if curve_list[j][4]==0:
			# Parametrized by R
			Rsamples = np.linspace(curve_list[j][6],curve_list[j][7],Ncsamples)
			tsamples = np.zeros(Ncsamples)
			for jj in range(Ncsamples):
				v2d_guess = curve_list[j][5](Rsamples[jj])
				v2dr = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,Rsamples[jj],[ocase,icase]))
				if v2dr.success:
					v2d = v2dr.x
					v0 = get_v0_of_vabls(v2d,ocase,icase)
					tsamples[jj] = np.arcsin(0.5*v0[1])
				else:
					print "Error computing theta at R="+str(Rsamples[jj])
		else:
			# Parametrized by theta
			Rsamples = np.zeros(Ncsamples)
			tsamples = np.linspace(curve_list[j][6],curve_list[j][7],Ncsamples)
			for jj in range(Ncsamples):
				v2dR_guess = curve_list[j][5](tsamples[jj])
				v2dRr = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,tsamples[jj],icase))
				if v2dRr.success:
					v2dR = v2dRr.x
					Rsamples[jj] = np.copy(v2dR[4])
					#v02d = v2dR[:2]+v2dR[2:4]
					#v12d = v2dR[2:4]-v2dR[:2]
					#v0 = x3d(v02d,Rsamples[jj])
					#v1 = x3d(v12d,Rsamples[jj])
				else:
					print "Error computing R(theta) at R="+str(tsamples[jj])
		if j not in curve_emphasis:
			plt.plot(Rsamples,tsamples,c='k',linewidth=0.3)
		else:
			plt.plot(Rsamples,tsamples,c='m',linewidth=1.0)
	print "Overlaying the ground state curve"
	plt.scatter(Rtsamples_gs[:,0],Rtsamples_gs[:,1],s=1.5,c='r')
	plt.xlabel("Cylinder radius")
	plt.ylabel("Crystal orientation")
	bounds = np.array([[R_min,0],[R_max,np.pi/3]])
	plt.plot([R_min,R_min],[0,np.pi/3],c='b',linewidth=0.2)
	plt.plot([R_max,R_max],[0,np.pi/3],c='b',linewidth=0.2)
	plt.savefig("ground_state_theta_R_"+postfix+".png",dpi=1024)
	pars = np.array([depl_str,dpr])
	np.savetxt("ground_state_theta_R_"+postfix+".par",pars)
	plt.show()

def plot_curves_energies(depl_str,dpr,Nsamples):
	# Make a 3D line plot in (R,theta,gamma) space
	V_min = V_depl_min(depl_str,dpr)
	for j in range(len(curve_list)):
		print "Plotting curve "+str(j)
		if curve_list[j][4]==0:
			# Parametrized by R
			Rsamples = np.linspace(curve_list[j][10],curve_list[j][11],Nsamples)
			tsamples = np.zeros(Nsamples)
			Esamples = np.zeros(Nsamples)
			for jj in range(Nsamples):
				v2d_guess = curve_list[j][5](Rsamples[jj])
				v2dr = scipy.optimize.root(ls_cconstraint,v2d_guess)
				if v2dr.success:
					v2d = v2dr.x
					tsamples[jj] = np.arcsin(0.5*v2d[1])
					Esamples[jj] = get_ls_tension(v2d,Rsamples[jj],depl_str,dpr,V_min)
				else:
					tsamples[jj] = 1e99
					Esamples[jj] = 1e99
		else:
			# Parametrized by theta
			Rsamples = np.zeros(Nsamples)
			tsamples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
			Esamples = np.zeros(Nsamples)
			for jj in range(Nsamples):
				v2dR_guess = curve_list[j][5](tsamples[jj])
				v2dRr = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess)
				if v2dRr.success:
					v2dR = v2dRr.x
					Rsamples[jj] = v2dR[4]
					v02d = v2dR[:2]+v2dR[2:4]
					v12d = v2dR[2:4]-v2dR[:2]
					v0 = x3d(v02d,Rsamples[jj])
					v1 = x3d(v12d,Rsamples[jj])
					d0sq = np.dot(v0,v0)
					d1sq = np.dot(v1,v1)
					Esamples[jj] = -V_min
					if d1sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d1sq),depl_str,dpr)
					if d0sq<rthrsq:
						Esamples[jj]+=V_depl(np.sqrt(d0sq),depl_str,dpr)
					Esamples[jj]/=abs(v2dR[1])
				else:
					Esamples[jj] = 1e99
					Rsamples[jj] = 1e99
		# Plot the current curve segment (where it is defined)
		valid = np.where(Esamples<1e50)
		fig = plt.figure()
		#ax = fig.add_subplot(111,projection='3d')
		#Axes3D.plot(Rsamples[valid],tsamples[valid],Esamples[valid])
	plt.show()

Rmin = 1e99
Rmax = 0
for j in range(len(curve_list)):
	if min(curve_list[j][10:12])<Rmin:
		Rmin = min(curve_list[j][10:12])
	if max(curve_list[j][10:12])>Rmax:
		Rmax = max(curve_list[j][10:12])
print [Rmin,Rmax]	
#print "Computing line slip tensions (per unit height)"
#compute_ls_tensions(1.0,0.05,200)
#print "Computing domain tensions per unit height"
#compute_domain_tensions(curve_list,1.0,0.05,500)

#print "Plotting tension curves"
#plot_ls_tensions()
#print "Plotting full tension"
#plot_cr_tensions(1.0,0.05)

#print "Computing ground state curve"
#print "Computing ground state tensions"
#compute_full_gs_tensions(curve_list,1.0,0.05,Rmin,Rmax,1024)
