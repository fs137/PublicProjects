# AUTHOR: Henry Wilkin
import numpy as np
import scipy.optimize
import scipy.interpolate
from matplotlib import pyplot as plt
from enumeration import *
import time
import os

hrt3 = 0.5*np.sqrt(3)

def x3d(x2d,R):
	Rinv = 1./R
	rval = np.zeros(3)
	rval[2]=x2d[1]
	rval[0]=R*np.cos(x2d[0]*Rinv)
	rval[1]=R*np.sin(x2d[0]*Rinv)
	return rval

def v3d(x2d1,x2d2,R):
	x3d1 = x3d(x2d1,R)
	x3d2 = x3d(x2d2,R)
	rval = x3d2-x3d1
	return rval

def R_criterion_v2d(R,v):
	v3d_ = x3d(v,R)
	v3d_[0]-=R
	return np.dot(v3d_,v3d_)-4

def get_R_from_v2d(v2d):
	# Use the representative vector with largest azimuthal component
	azimuthal = np.array([v2d[0],v2d[2],v2d[4]])
	azimuthal*=azimuthal
	j = np.argmax(azimuthal)
	j*=2
	R_root = scipy.optimize.root(R_criterion_v2d,1.0,(v2d[j:j+2]))
	if R_root.success:
		return R_root.x
	else:
		print "Failed to obtain valid cylinder radius from "+str(v2d[j:j+2])
		return None

def x2d_of_x3d(x3d_):
	R = np.linalg.norm(x3d_[:2])
	theta = np.arccos(x3d_[0]/R)
	theta=theta if theta>0 else -theta
	return np.array([R*theta,x3d_[2]])

# idx datatype:
idx = np.dtype([('id','i8'),('x','f8'),('y','f8'),('z','f8')])
# (Redundant)
def load_interp_data():
	try:
		Rsamples = np.loadtxt("LSDATA/ls_Rsamples.dat")
		tsamples = np.loadtxt("LSDATA/ls_thetas.dat")
		v2Ds = np.loadtxt("LSDATA/ls_v2Ds.dat")
		curve_features = np.loadtxt("LSDATA/ls_curvef.dat",dtype=int)
		curve_featuresc = np.loadtxt("LSDATA/ls_curvefc.dat")
	except:
		print "Error: one or more data files couldn't be found. Run enumeration code first."
		quit()
	return Rsamples,tsamples,v2Ds,curve_features,curve_featuresc

def circle(r,npts):
	dtheta = 2*np.pi/npts
	generator = np.exp(dtheta*1j)
	xs = np.zeros([npts+1,2])
	cpt = r
	for j in range(npts+1):
		xs[j,0] = cpt.real
		xs[j,1] = cpt.imag
		cpt*=generator
	return xs

# Function to evaluate curves exactly
def exact_v2dR_theta(curve_entry,theta):
	# First make sure that theta is within the required bounds:
	if min(curve_entry[6:8])<=theta<=max(curve_entry[6:8]):
		v2dR_guess = curve_entry[5](theta)
		[na,nb,ocase,icase] = curve_entry[:4]
		v2dRr = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,theta,icase))
		if v2dRr.success:
			return v2dRr.x
		else:
			print "Error (exact_v2dR_theta): failed to converge from guess"
	else:
		print "Error (exact_v2dR_theta): value of theta out of bounds"
		quit()

def exact_v2d_R(curve_entry,R):
	if min(curve_entry[6:8])<=R<=max(curve_entry[6:8]):
		v2d_guess = curve_entry[5](R)
		[na,nb,ocase,icase]=curve_entry[:4]
		v2dr = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,R,[ocase,icase]))
		if v2dr.success:
			return v2dr.x
		else:
			print "Error (exact_v2d_R): failed to converge from guess"
	else:
		print "Error (exact_v2d_R): value of R out of bounds"
		quit()

# Obtain an array of exact v2d samples parametrized by theta or R:
def exact_v2d_samples(curve_entry,Nsamples):
	samples = []
	if curve_entry[4]==1:
		tsamples = np.linspace(curve_entry[6],curve_entry[7],Nsamples)
		for i in range(Nsamples):
			v2dR = exact_v2dR_theta(curve_entry,tsamples[i])
			v2d = get_v2d_of_v2dR(v2dR,tsamples[i],curve_entry[2:4])
			samples.append([v2dR[4],v2d,curve_entry[:2]])
	else:
		Rsamples = np.linspace(curve_entry[6],curve_entry[7],Nsamples)
		for i in range(Nsamples):
			v2d = exact_v2d_R(curve_entry,Rsamples[i])
			samples.append([Rsamples[i],v2d,curve_entry[:2]])
	return samples

# Function to define an initial set of functions used to provide initial guesses
# for root finding
#	NOTE: we would like to use the resulting list of curves 
#	to construct lattices with various orientations at a 
#	specific radius R. Hence, it would be useful to include
#	the minimum and mmaximum R for each curve, together with
#	'multivalued domains' where applicable: the domain over 
#	which a given curve includes two orientations for a given
#	radius (case 3(ab))
def make_guess_curves(Rsamples,tsamples,v2Ds,curve_features,curve_featuresc):
	# Determine the curve type and limits from curve_features:
	# Reconstruct the list of curves using the established limits:
	curve_list = []
	for i in range(len(curve_features)):
		imin = curve_features[i,4]
		imax = curve_features[i,5]
		if imax-imin<2:
			print "Problem encountered at curve "+str(i)+", case="+str([ocase,icase])
		na = curve_features[i,0]
		nb = curve_features[i,1]
		[n0,n5,n0p,n5p] = curve_features[i,7:11]
		ip = curve_features[i,11]
		ocase = curve_features[i,2]
		icase = curve_features[i,3]
		branch_state = curve_features[i,6]
		# Avoid redundant curves:
		R_pt1 = curve_featuresc[i,2]
		R_pt2 = curve_featuresc[i,3]
		R2valmin = curve_featuresc[i,4]
		R2valmax = curve_featuresc[i,5]
		theta_vm = curve_featuresc[i,6]
		#print [ocase,icase,na,nb,imax-imin]
		#plt.scatter(Rsamples[imin:imax,0],tsamples[imin:imax])
		#plt.show()
		#plt.clf()
		# Check if Rsamples[imin:imax,0] is ordered:
		radii = np.copy(Rsamples[imin:imax,0])
		latvec = np.copy(v2Ds[imin:imax,:])
		thetas = np.copy(tsamples[imin:imax])
		if curve_features[i,6]==0:
			if radii[-1]>radii[0]:
				pass
			else:
				radii= radii[::-1]
				latvec = latvec[::-1,:]
			curve = scipy.interpolate.interp1d(radii,latvec,axis=0,kind='cubic',fill_value="extrapolate")
			curve_list.append([na,nb,ocase,icase,0,curve,Rsamples[imin,0],Rsamples[imax-1,0],
				Rsamples[imin,0],Rsamples[imax-1,0],Rsamples[imin,0],Rsamples[imax-1,0],-10,-10,-10,[n0,n5,n0p,n5p],ip])
		if curve_features[i,6]==1:
			# Only parametrize by theta if it is necessary
			if thetas[-1]>thetas[0]:
				pass
			else:
				thetas = thetas[::-1]
				latvec = latvec[::-1,:]
				radii = radii[::-1]
			v2DRs = np.zeros([imax-imin,5])
			# first generate v2DRs from v2Ds:
			for ii in range(imax-imin):
				v2DRs[ii,:2] = np.copy(latvec[ii,:2])
				v2DRs[ii,2:4] = np.copy(latvec[ii,4:6])
				v2DRs[ii,4] = radii[ii]
			# interpolate over tsamples[imin:imax]:
			if thetas[0]<thetas[-1]:
				pass
			else:
				print "Error: endpoints of theta interval are reversed"
			curve = scipy.interpolate.interp1d(thetas,v2DRs,axis=0,kind='cubic',fill_value="extrapolate")
			curve_list.append([na,nb,ocase,icase,1,curve,tsamples[imin],tsamples[imax-1], # 0 through 7
				curve_featuresc[i,0],curve_featuresc[i,1],R_pt1,R_pt2,R2valmin, # 8 through 12
				R2valmax,theta_vm,[n0,n5,n0p,n5p],ip]) # 13 through 15
			# NOTE: special care is needed near 'critical points', where v2D(R) becomes multivalued
			# Since endpoints are determined by commensurate crystals, this amounts to estimating  
			# the maximum of R(theta). 
	return curve_list

# Automatically load curve data and construct guess curves:
print "Loading data for interpolation"
[Rsamples,tsamples,v2Ds,curve_features,curve_featuresc]=load_interp_data()
print "Generating approximate curves for guesses"
curve_list = make_guess_curves(Rsamples,tsamples,v2Ds,curve_features,curve_featuresc) 
print "(done)"
# Determine nmax:
nmax_global = 0
for j in range(len(curve_list)):
	[n0,n5,n0p,n5p] = curve_list[j][15]
	nmax_global = max(n0,n5,nmax_global)
	
def make_ccls_plot(mode=0,Nsamples = 100,pfix="",curve_emphasis=[],labels=False):
	marked = {}
	for j in range(len(curve_list)):
		if j in marked:
			continue
		[na,nb,ocase,icase] = curve_list[j][:4]
		j_dual = curve_list[j][-1]
	 	marked[j_dual]=True
		case = [ocase,icase]
		if mode==0:
			# (R,theta) plot
			if curve_list[j][4]==0:
				Rsamples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
				tsamples = np.zeros(Nsamples)
				for jj in range(Nsamples):
					v2d_guess = curve_list[j][5](Rsamples[jj])
					v2d_exact = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,Rsamples[jj],case)).x 
					v0 = get_v0_of_vabls(v2d_exact,ocase,icase)
					tsamples[jj] = np.arcsin(0.5*v0[1])
			elif curve_list[j][4]==1:
				Rsamples = np.zeros(Nsamples)
				tsamples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
				for jj in range(Nsamples):
					v2dR_guess = curve_list[j][5](tsamples[jj])
					v2dR_exact = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,tsamples[jj],icase)).x
					Rsamples[jj] = np.copy(v2dR_exact[4])
			# Plot Rsamples.tsamples
			if j not in curve_emphasis:
				plt.plot(Rsamples,tsamples,linewidth=0.4,c='b')
			else:
				plt.plot(Rsamples,tsamples,linewidth=1.0,c='r')
			if labels:
				if j_dual>-1:
					plt.text(Rsamples[Nsamples/4],tsamples[Nsamples/4],str(j)+" and "+str(j_dual),fontsize=2)
				else:
					plt.text(Rsamples[Nsamples/4],tsamples[Nsamples/4],str(j),fontsize=2)
			
	plt.savefig("lscc_curve_map_nmax"+str(nmax_global)+pfix+".png",dpi=1024)
	plt.clf()

# Plot the angles made by the line slip generator
# (Assume vector has positive 
def plot_ls_angle(Nsamples=100):
	for j in range(len(curve_list)):
		[na,nb,ocase,icase] = curve_list[j][:4]
		case = [ocase,icase]
		# (R,theta) plot
		if curve_list[j][4]==0:
			Rsamples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
			tsamples = np.zeros(Nsamples)
			for jj in range(Nsamples):
				v2d_guess = curve_list[j][5](Rsamples[jj])
				v2d_exact = scipy.optimize.root(ls_cconstraint,v2d_guess,(na,nb,Rsamples[jj],case)).x 
				
				tsamples[jj] = np.arcsin(0.5*v2d_exact[1])
				if v2d_exact[0]<0:
					tsamples[jj]*=-1
				if tsamples[jj]<0:
					tsamples[jj]+=np.pi
				
		elif curve_list[j][4]==1:
			Rsamples = np.zeros(Nsamples)
			t0samples = np.linspace(curve_list[j][6],curve_list[j][7],Nsamples)
			tsamples = np.zeros(Nsamples)
			for jj in range(Nsamples):
				v2dR_guess = curve_list[j][5](t0samples[jj])
				v2dR_exact = scipy.optimize.root(ls_cconstraint_theta,v2dR_guess,(na,nb,t0samples[jj],icase)).x
				Rsamples[jj] = np.copy(v2dR_exact[4])
				tsamples[jj] = np.arcsin(0.5*v2dR_exact[1])
				if v2dR_exact[0]<0:
					tsamples[jj]*=-1
				if tsamples[jj]<0:
					tsamples[jj]+=np.pi
		# Plot Rsamples.tsamples
		plt.plot(Rsamples,tsamples,linewidth=0.7)
		plt.text(Rsamples[Nsamples/4],tsamples[Nsamples/4],"curve "+str(j),fontsize=3)
	plt.savefig("ls_angles_nmax"+str(nmax_global)+".png",dpi=1024)
	plt.clf()

# NOTE: 8.11.2019 
#	Using 'dict' mode produces fewer orientations than 'list' mode. (This is because when v2dR is parametrized by theta, 
#	there can be multiple orientations associated with the given curve index.)
# Find all closed crystals with radius R:
#	(consider storing crystal orientations with a dictionary structure instead)
def get_closed_crystals(R,curve_list,mode='list',curve_index=False):
	if mode=='list':
		rel_curves = [] # a list of 'relevant' curves
		measured = {}
		for i in range(len(curve_list)):
			# Check if the associated 'edge' has already been accounted for:
			ci = tuple(curve_list[i][-2])
			reverse_curve_indices=tuple([ci[2],ci[3],ci[0],ci[1]])
			if reverse_curve_indices in measured:
				continue
			measured[ci]=True
			measured[reverse_curve_indices]=True
			if curve_list[i][4]==0:
				R1 = curve_list[i][10]
				R2 = curve_list[i][11]
				Rmin = min(R1,R2)
				Rmax = max(R1,R2)
				if Rmin<R<Rmax:
					rel_curves.append(i)
			elif curve_list[i][4]==1:
				Rmin = min(curve_list[i][10],curve_list[i][11])
				Rmax = curve_list[i][13]
				if Rmin<R<Rmax:
					rel_curves.append(i)
					# rel_curves.append(curve_list[i])
				elif max(curve_list[i][10],curve_list[i][11])>Rmax:
					print "Error: R2valmax should equal or exceed the radius of the two endpoints"
					quit()
		# Loop over relevant curves, and solve for the *exact* lattice generators (possibly for each branch.)
		ex_list = []
		for ri in range(len(rel_curves)):
			i = rel_curves[ri]
			#na = rel_curves[i][0]
			#nb = rel_curves[i][1]
			#[n0,n5,n0p,n5p] = rel_curves[i][15]
			#ocase = rel_curves[i][2]
			#icase = rel_curves[i][3]
			#branched = rel_curves[i][4]
			na = curve_list[i][0]
			nb = curve_list[i][1]
			[n0,n5,n0p,n5p] = curve_list[i][15]
			ocase = curve_list[i][2]
			icase = curve_list[i][3]
			branched = curve_list[i][4]
			# Check if the current curve is associated with two branches:
			if branched==1:
				# Check if R is in the multibranched domain:
				#theta_min = rel_curves[i][6]
				#theta_max = rel_curves[i][7]
				#theta_mid = rel_curves[i][14]
				theta_min = curve_list[i][6]
				theta_max = curve_list[i][7]
				theta_mid = curve_list[i][14]
				if curve_list[i][12]<R<curve_list[i][13]:
					# Solve for orientations in both branches (use rel_curves[i][14])
					theta_min = curve_list[i][6]
					theta_max = curve_list[i][7]
					theta_mid = curve_list[i][14]
					if theta_min>theta_mid or theta_max<theta_min or theta_max<theta_mid:
						print "Error (regen.py): bounds "+str([theta_min,theta_mid,theta_max])+" are not ordered"
						quit()
					# Find the maximum of 'curve' in the interval theta_min, theta_max,
					# (and check if it exceeds R)
					# First define inline function returning just the radius:
					temp_Rtcurve = lambda theta : curve_list[i][5](theta)[4]-R
					t_root1 = scipy.optimize.brentq(temp_Rtcurve,theta_min+1e-9,theta_mid)
					t_root2 = scipy.optimize.brentq(temp_Rtcurve,theta_mid,theta_max-1e-9)
					v2dR1 = curve_list[i][5](t_root1)
					v2dR2 = curve_list[i][5](t_root2)
					v2d1 = np.zeros(6)
					v2d2 = np.zeros(6)
					v2d1[:2] = np.copy(v2dR1[:2])
					v2d1[4:6] = np.copy(v2dR1[2:4])
					projxy = 2*np.cos(t_root1)
					v02d1 = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root1)])
					projxy = 2*np.cos(t_root2)
					v02d2 = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root2)])
					if curve_list[i][2:4]==[2,0]:
						v2d1[2:4] = v02d1
						v2d2[2:4] = v02d2
						
					elif curve_list[i][2:4]==[2,1]:
						v2d1[2:4] = -v2d1[:2]+v02d1
						v2d2[2:4] = -v2d2[:2]+v02d2
					else:
						print "Error: branches only occur in the third external case"
						quit()
					# Refine the solutions with ls_cconstraint
					v2d1r = scipy.optimize.root(ls_cconstraint,v2d1,(na,nb,R,[ocase,icase]))
					v2d2r = scipy.optimize.root(ls_cconstraint,v2d2,(na,nb,R,[ocase,icase]))
					if v2d1r.success:
						if curve_index:
							entry = [v2d1r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p],i]
						else:
							entry = [v2d1r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]]
						ex_list.append(entry)
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,1])
					# append the two orientations to ex_list:
					if v2d2r.success:
						if curve_index:
							entry = [v2d2r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p],i]
						else:
							entry = [v2d2r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]]
						ex_list.append(entry)
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,2])
				else:
					# Solve for the single orientation at angle theta
					temp_Rtcurve = lambda theta : curve_list[i][5](theta)[4]-R
					t_root = scipy.optimize.brentq(temp_Rtcurve,theta_min,theta_max-1e-9)
					v2dR = curve_list[i][5](t_root)
					v2d = np.zeros(6)
					v2d[:2] = np.copy(v2dR[:2])
					v2d[4:6] = np.copy(v2dR[2:4])
					projxy = 2*np.cos(t_root)
					v02d = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root)])
					if curve_list[i][2:4]==[2,0]:
						v2d[2:4] = v02d
					elif curve_list[i][2:4]==[2,1]:
						v2d[2:4] = -v2d[:2]+v02d
					else:
						print "Error: branches only occur in the third external case"
						quit()
					# Refine the solutions with ls_cconstraint
					v2dr = scipy.optimize.root(ls_cconstraint,v2d,(na,nb,R,[ocase,icase]))
					if v2dr.success:
						if curve_index:
							entry = [v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p],i]
						else:
							entry = [v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]]
						ex_list.append(entry)
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,1])

			else:
				# use curve_list[i][5] to provide an initial guess for v2d:
				# Parametrized by R
				v2d = curve_list[i][5](R)
				# refine the initial guess using ls_cconstraint:
				v2dr = scipy.optimize.root(ls_cconstraint,v2d,(na,nb,R,[ocase,icase]))
				if v2dr.success:
					if curve_index:
						entry = [v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p],i]
					else:
						entry = [v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]]
					ex_list.append(entry)
				else:
					print "Failed to find orientation at [na,nb,ocase,icase]="+str([na,nb,ocase,icase])
	elif mode=='dict':
		rel_curves = {} # a list of 'relevant' curves 
				# 	(curves that could feasibly attain the given cylinder radius.)
		measured = {}
		for i in range(len(curve_list)):
			# Check if the associated 'edge' has already been accounted for:
			ci = tuple(curve_list[i][-2])
			reverse_curve_indices=tuple([ci[2],ci[3],ci[0],ci[1]])
			if reverse_curve_indices in measured:
				continue
			measured[ci]=True
			measured[reverse_curve_indices]=True
			if curve_list[i][4]==0:
				R1 = curve_list[i][10]
				R2 = curve_list[i][11]
				Rmin = min(R1,R2)
				Rmax = max(R1,R2)
				if Rmin<R<Rmax:
					rel_curves[i]=curve_list[i]
			elif curve_list[i][4]==1:
				Rmin = min(curve_list[i][10],curve_list[i][11])
				Rmax = curve_list[i][13]
				if Rmin<R<Rmax:
					rel_curves[i]=curve_list[i]
				elif max(curve_list[i][10],curve_list[i][11])>Rmax:
					print "Error: R2valmax should equal or exceed the radius of the two endpoints"
					quit()
		# Loop over relevant curves, and solve for the *exact* lattice generators (possibly for each branch.)
		ex_list = {}
		for i in rel_curves:
			na = rel_curves[i][0]
			nb = rel_curves[i][1]
			[n0,n5,n0p,n5p] = rel_curves[i][15]
			ocase = rel_curves[i][2]
			icase = rel_curves[i][3]
			branched = rel_curves[i][4]
			
			# Check if the current curve is associated with two branches:
			if rel_curves[i][4]==1:
				# Check if R is in the multibranched domain:
				theta_min = rel_curves[i][6]
				theta_max = rel_curves[i][7]
				theta_mid = rel_curves[i][14]
				if rel_curves[i][12]<R<rel_curves[i][13]:
					# Solve for orientations in both branches (use rel_curves[i][14])
					theta_min = rel_curves[i][6]
					theta_max = rel_curves[i][7]
					theta_mid = rel_curves[i][14]
					if theta_min>theta_mid or theta_max<theta_min or theta_max<theta_mid:
						print "Error (regen.py): bounds "+str([theta_min,theta_mid,theta_max])+" are not ordered"
						quit()
					# Find the maximum of 'curve' in the interval theta_min, theta_max,
					# (and check if it exceeds R)
					# First define inline function returning just the radius:
					temp_Rtcurve = lambda theta : rel_curves[i][5](theta)[4]-R
					t_root1 = scipy.optimize.brentq(temp_Rtcurve,theta_min+1e-9,theta_mid)
					t_root2 = scipy.optimize.brentq(temp_Rtcurve,theta_mid,theta_max-1e-9)
					v2dR1 = rel_curves[i][5](t_root1)
					v2dR2 = rel_curves[i][5](t_root2)
					v2d1 = np.zeros(6)
					v2d2 = np.zeros(6)
					v2d1[:2] = np.copy(v2dR1[:2])
					v2d1[4:6] = np.copy(v2dR1[2:4])
					projxy = 2*np.cos(t_root1)
					v02d1 = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root1)])
					projxy = 2*np.cos(t_root2)
					v02d2 = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root2)])
					if rel_curves[i][2:4]==[2,0]:
						v2d1[2:4] = v02d1
						v2d2[2:4] = v02d2
						
					elif rel_curves[i][2:4]==[2,1]:
						v2d1[2:4] = -v2d1[:2]+v02d1
						v2d2[2:4] = -v2d2[:2]+v02d2
					else:
						print "Error: branches only occur in the third external case"
						quit()
					# Refine the solutions with ls_cconstraint
					v2d1r = scipy.optimize.root(ls_cconstraint,v2d1,(na,nb,R,[ocase,icase]))
					v2d2r = scipy.optimize.root(ls_cconstraint,v2d2,(na,nb,R,[ocase,icase]))
					if v2d1r.success:
						add_dl(ex_list,i,[v2d1r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]])
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,1])
					# append the two orientations to ex_list:
					if v2d2r.success:
						add_dl(ex_list,i,[v2d2r.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]])
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,2])
				else:
					# Solve for the single orientation at angle theta
					temp_Rtcurve = lambda theta : rel_curves[i][5](theta)[4]-R
					t_root = scipy.optimize.brentq(temp_Rtcurve,theta_min,theta_max-1e-9)
					v2dR = rel_curves[i][5](t_root)
					v2d = np.zeros(6)
					v2d[:2] = np.copy(v2dR[:2])
					v2d[4:6] = np.copy(v2dR[2:4])
					projxy = 2*np.cos(t_root)
					v02d = np.array([2*R*np.arcsin(projxy/(2*R)),2.*np.sin(t_root)])
					if rel_curves[i][2:4]==[2,0]:
						v2d[2:4] = v02d
					elif rel_curves[i][2:4]==[2,1]:
						v2d[2:4] = -v2d[:2]+v02d
					else:
						print "Error: branches only occur in the third external case"
						quit()
					# Refine the solutions with ls_cconstraint
					v2dr = scipy.optimize.root(ls_cconstraint,v2d,(na,nb,R,[ocase,icase]))
					if v2dr.success:
						add_dl(ex_list,i,[v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]])
					else:
						print "Failed to find orientation at [na,nb,ocase,icase,branch]="+str([na,nb,ocase,icase,1])
					
			else:
				# use rel_curves[i][5] to provide an initial guess for v2d:
				# Parametrized by R
				v2d = rel_curves[i][5](R)
				# refine the initial guess using ls_cconstraint:
				v2dr = scipy.optimize.root(ls_cconstraint,v2d,(na,nb,R,[ocase,icase]))
				if v2dr.success:
					add_dl(ex_list,i,[v2dr.x,na,nb,ocase,icase,[n0,n5,n0p,n5p]])
				else:
					print "Failed to find orientation at [na,nb,ocase,icase]="+str([na,nb,ocase,icase])

	return ex_list

def positive_vectors(v2d):
	rval = np.copy(v2d)
	vadvb = np.dot(v2d[:2],v2d[2:4])
	if vadvb>0:
		vc = v2d[:2]-v2d[2:4]
	else:
		vc = v2d[:2]+v2d[2:4]
	if np.dot(rval[4:],rval[:2])<0:
		rval[4:]*=-1
	if np.dot(rval[2:4],rval[:2])<0:
		rval[2:4]*=-1
	if np.dot(rval[2:4],rval[4:])<0:
		rval[2:4]=vc
	return rval

def gls_distances(R,curve_list):
	ex_list = get_closed_crystals(R,curve_list)
	gls_dists = []
	for j in range(len(ex_list)):
		v2d = positive_vectors(ex_list[j][0])
		vc = v2d[:2]-v2d[2:4]
		gls1_0 = v2d[4:]+v2d[:2] 
		gls1_gen = vc
		gls1_0-=2*vc
		gls2_0 = v2d[4:]-v2d[:2]
		gls2_gen = v2d[2:4]
		gls2_0-=2*v2d[2:4]
		dist3d = lambda dv1: np.linalg.norm(x3d(dv1,R))
		gls_seps = np.zeros(5)
		for j in range(5):
			gls_seps[j]=dist3d(gls1_0+j*gls1_gen)
		gls1_sep = np.min(gls_seps)
		for j in range(5):
			gls_seps[j]=dist3d(gls2_0+j*gls2_gen)
		gls2_sep = np.min(gls_seps)
		dv1 = gls1_sep
		dv2 = gls2_sep
		if dv1<dv2:
			gls_dists.append([dv1,dv2])
		else:
			gls_dists.append([dv2,dv1])
	return gls_dists
# Plot curves with the list of exact orientations at radius 2.5:
def plot_curves_and_lattices(curve_list,ex_list):
	for i in range(len(curve_list)):
		case123 = curve_list[i][2]
		caseab = curve_list[i][3]
		if curve_list[i][4]==0:
			# use curve_list[i][5] to generate points on the curve:
			R1 = curve_list[i][6]
			R2 = curve_list[i][7]
			Rmin = min(R1,R2)
			Rmax = max(R1,R2)
			#radii = np.linspace(Rmin+0.01*(Rmax-Rmin),Rmax-0.01*(Rmax-Rmin),100)
			radii = np.linspace(Rmin,Rmax,100)
			thetas = np.zeros(100)
			for ii in range(100):
				v2d=curve_list[i][5](radii[ii])
				v02d = get_v0_of_vabls(v2d,case123,caseab)
				thetas[ii] = np.arcsin(0.5*v02d[1])
			plt.plot(radii,thetas)
		elif curve_list[i][4]==1:
			#curve_list.append([na,nb,ocase,icase,1,curve,tsamples[imin],tsamples[imax-1],
			#	curve_featuresc[i,0],curve_featuresc[i,1],curve_featuresc[i,2],curve_featuresc[i,3],curve_featuresc[i,4],
			#	curve_featuresc[i,5],curve_featuresc[i,6]])
			t1 = curve_list[i][6]
			t2 = curve_list[i][7]
			tmin = min(t1,t2)
			tmax = max(t1,t2)
			#thetas = np.linspace(tmin+0.01*(tmax-tmin),tmax-0.01*(tmax-tmin),100)
			thetas = np.linspace(tmin,tmax,100)
			radii = np.zeros(100)
			for ii in range(100):
				v2dr = curve_list[i][5](thetas[ii])
				radii[ii] = v2dr[4]
			plt.plot(radii,thetas)
	for i in range(len(ex_list)):
		v02d = get_v0_of_vabls(ex_list[i][0],ex_list[i][3],ex_list[i][4])
		# check if v02d[0]<0:
		if v02d[0]<0:
			print "Error: v02d has wrong orientation"
		theta = np.arcsin(0.5*v02d[1])
		plt.scatter([2.5],[theta])
	plt.savefig("ex_lattice_test.png")
	plt.clf()

# Functions to define elementary defect properties: distances associated with the 
#	bond broken by the line slip, and the separation across a gapped line slip defect.

# Get a standard set of vectors, with similar local properites independent of the 'parity' 
# 	of the crystal
def standard_vectors(v2d,pars=None):
	if pars==None:
		pars_defined = False
	else:
		pars_defined = True
		[na,nb]=pars
	# If na, nb are given, also find values nua, nu2 such that nua*ua+nu2*u2+u3 forms a closed loop
	# 	given ua = va, u2 = (ca*va+cb*vb), and u3 = cls*vls, where the 'c' coeffs are all signs, 
	#	we have 
	#		nua*va+nu2*(ca*va+cb*vb)+c3*cls*vls = na*va+nb*vb+vls, or
	#		nua+nu2*ca = na,   nu2*cb = nb,  and  cls = c3 
	#	or,
	#		nu2 = nb*cb,   nua = na-ca*cb*nb, 
	# v2d[:2]: va, v2d[2:4]: vb, v2d[4:]: vls
	ua = np.copy(v2d[:2])
	# Choose a down-ward pointing line slip generator if possible
	#aux_phase = 1
	#if ua[1]>0:
	#	ua*=-1
	#	aux_phase = -1
	ub = np.copy(v2d[2:4])
	vls = v2d[4:]
	vlsdua = np.dot(vls,ua)
	cls = 2*int(vlsdua>=0)-1
	u3 = vls*cls
	# Define the two vectors u1,u2 associated with gapped line slip generators,
	#	as well as the vector u3 associated with crossing the line slip
	# NOTE: u1 is ostensibly the fast vector, and u2 the slow one
	uadub = np.dot(ua,ub)
	c2b = 2*int(uadub<0)-1
	c2a = 1
	c1b = -c2b
	c1a = 0
	u1 = c1b*ub
	u2 = ua+c2b*ub
	# nua*va+nu2*(va-vb) = na*va+nb*vb => nua= na-nu2, nu2 =-nb
	# Sort vectors u1 and u2 according to their alignment with u3:
	# 	(the generator associated with the fast 'kink' has a positive 
	#	dot product with u3. [conj.])
	if np.dot(u3,u1)<0:
		# Exchange the fast and slow vectors
		aux_u = u1
		u1 = u2
		u2 = aux_u
		[c2a,c2b] = [c1a,c1b] # (c1a and c1b are not required at this stage)
	if not pars_defined:
		return [ua,u1,u2,u3]
	else:
		# Compute nua, nu2 from na, nb, c2a, c2b, cls:
		nu2 = nb*c2b
		nua = na-c2a*nu2 # Check this
		return [ua,u1,u2,u3,nua,nu2,cls] 
		

def get_defect_geos(v2d,R,case):
	[ua,u1,u2,u3] = standard_vectors(v2d)
	# Compute geometric properties of features in crystals
	# Gaps associated with the line slip defect:
	gv1 = ua-u3
	gv2 = ua+u3
	g1 = np.linalg.norm(gv1)
	g2 = np.linalg.norm(gv2)
	# Gaps associated with 'gapped' line slip defects:
	lsgv1 = u3-ua-u1
	lsgv2 = -u3-ua-u2
	lsg1 = np.linalg.norm(lsgv1)
	lsg2 = np.linalg.norm(lsgv2)
	# Compute the channel width associated with each g.l.s. generator:
	cwv1 = u3+u1
	cwv2 = u3-u2
	cw1 = np.linalg.norm(cw1)
	cw2 = np.linalg.norm(cw2)
	# Compute the expected distance between adjacent f.v. minima:
	fvv1 = u3-u1
	fvv2 = u3+u2
	fv1 = np.linalg.norm(fvv1)
	fv2 = np.linalg.norm(fvv2)
	# Maintain a single edge reference list:
	basic_len_list = [g1,g2,lsg1,lsg2,fv1,fv2,cw1,cw2]
	basic_v_list = [gv1,gv2,lsgv1,lsgv2,fvv1,fvv2,cwv1,cwv2]
	return basic_len_list, basic_v_list

def get_defect_geos(ex_list,R):
	# Get a list of orientations at radius R:
	for i in range(len(ex_list)):
		v2d = ex_list[i][0]
		case = ex_list[i][3:]
		[bll, bvl] = get_defect_geos(v2d,R,case)
		ex_list[i].append(bll)
		ex_list[i].append(bvl)

# Generate points on a lattice defined by v2d at radius R:
def generate_pts(v2d,na,nb,R,z0,z1):
	# Starting at [0,z0], generate points along a crystal with line-slip.
	# generators: va, vb. (vls implicit)
	C = 2*np.pi*R
	height = abs(z1-z0)
	zmin = min(z0,z1)
	zmax = max(z0,z1)
        if v2d[1]>=-epsilon:
	    x_list = [np.array([0,zmin])] 
        else:
            x_list = [np.array([0,zmax])]
	Ncols = abs(nb)+1
	# a 'buffer' associated with corners in the crystal:
	corner_buffer = abs(na)+abs(nb)+2
	Nrows = 2*(int(height/abs(v2d[1]))+2*corner_buffer)
	x_array = np.zeros([Nrows,Ncols,2])
	x_array[0,0,:] = x_list[0]
	occ = np.zeros([Nrows,Ncols],dtype=int)
	# toroidal boundary conditions are used implicitly in numpy conventions,
	# though not with the 'twist' that is relevant in this case. 
	# If the line-slip is originally situated along, e.g., the zeroth column, then 
	# it is `twisted' from that column in the wrapped domain.
        llimR = -corner_buffer
	ulimR = Nrows-corner_buffer
	llimC = 0
	ulimC = Ncols
	for i in range(llimR,ulimR):
		for j in range(llimC,ulimC):
			# compute the height at site (i,j):
			x_array[i,j,:] = x_array[0,0,:]+v2d[:2]*i+v2d[2:4]*j
			if zmin<x_array[i,j,1]<zmax:
				occ[i,j] = 1
				x_list.append(x_array[i,j,:])		
	return [x_list,x_array,occ]

def display_fv_pts(new_xs,C):
	for j in range(len(new_xs)):
		while new_xs[j,0]<0:
			new_xs[j,0]+=C
		while new_xs[j,0]>C:
			new_xs[j,0]-=C
	plt.scatter(new_xs[:,0],new_xs[:,1],c='b')
	plt.axes().set_aspect('equal','datalim')
	plt.show()



def compute_twist(v2d,R,Na,N1,N2):
	# First, define the vectors associated with the fast and slow line slip defects
	[ua,u1,u2,u3] = standard_vectors(v2d)
	# Compute the angle swept out by the line slip defect from top to bottom:
	ls_angle = (Na+1)*ua[0]+N1*u1[0]+N2*u2[0] # One extra segment along ua to close the loop
	ls_angle/=R
	# The 'twist' is then given by the angular deficit: 
	twist_ = 2*np.pi-(ls_angle%2*np.pi) # extra rotation needed to ensure closed loop within xy plane as well
	return twist_

def generate_layers(std_vecs,Nlayers):
	# Generate layers along u2, stacked along ua
	[ua,u1,u2,u3,nua,nu2,nu3] = std_vecs
	# Compute the number of rows needed to form a line slip defect with Nlayers contacts
	# nua*ua+nu2*u2+u3-[C,0] = [0,0]
	Nrows = Nlayers
	Ncols = abs(nu2)+1
	# Construct the pristine array
	x_array = np.zeros([Nrows,Ncols,2])
	occ = np.ones([Nrows,Ncols],dtype='int')
	for i in range(Nrows):
		for j in range(Ncols):	
			x_array[i,j,:] = i*ua+j*u2
	try:
		print occ.shape
	except:
		print type(occ)
		print occ
		quit()
	return [x_array,occ]

	

def images(index,v_shift,h_shift,dist=3):
	ext = 2*dist+1
	index_im = np.zeros([ext,ext,2],dtype=int)
	for i in range(ext):
		v_sign = i-dist
		for j in range(ext):
			h_sign = j-dist
			index_im[i,j,:] = index+v_sign*v_shift+h_sign*h_shift
	return index_im

def index_norm(index,R_lim,C_lim):
	return np.dot(index,index) 

sq_disps = [[0,0],[1,0],[0,1],[1,1],[0,-1],[-1,0],[-1,-1],[1,-1],[-1,1]]

def least_norm(index,v_shift,h_shift,R_lim,C_lim,dist=4):
	A = np.array([[v_shift[0],h_shift[0]],[v_shift[1],h_shift[1]]])
	rcoeffs = np.linalg.solve(A,index)
	icoeffs = np.zeros(2,dtype=int)
	icoeffs[0]=np.round(rcoeffs[0])
	icoeffs[1]=np.round(rcoeffs[1])
	#ln_index = -icoeffs[0]*v_shift-icoeffs[1]*h_shift+index
	#icoeffs = np.zeros([9,2],dtype=int)
	#icoeffs[0,0] = np.round(rcoeffs[0])
	#icoeffs[0,1] = np.round(rcoeffs[1])
	#for i in range(1,9):
	#	icoeffs[i,0] = icoeffs[0,0]+sq_disps[i][0]
	#	icoeffs[i,1] = icoeffs[0,1]+sq_disps[i][1]
	#cands = index-np.outer(icoeffs[:,0],v_shift)-np.outer(icoeffs[:,1],h_shift)
	cands = images(index-icoeffs[0]*v_shift-icoeffs[1]*h_shift,v_shift,h_shift,dist)
	ln_index = np.copy(index)
	ln = 1e99
	ln_set = False
	# Enumerate feasible candidates
	fcands = []
	for i in range(2*dist+1):
		for j in range(2*dist+1):
			if cands[i,j,0]>=0 and cands[i,j,1]>=0:
				fcands.append(cands[i,j,:])
	# Choose a standard candidate with lowest norm
	for i in range(len(fcands)):
		cand_ln = index_norm(fcands[i],R_lim,C_lim)
		if cand_ln<ln:
			ln_index = np.copy(fcands[i])
			ln = cand_ln
			ln_set = True
	if ln_set:
		# Find candidates with the least norm
		ln_cands = []
		for i in range(len(fcands)):
			if index_norm(fcands[i],R_lim,C_lim)==ln:
				ln_cands.append(fcands[i])
		if len(ln_cands)>1:
			print "Exceptional case"
			# Choose the candidate with smallest 0-component
			ln_0val = ln_index[0]
			for i in range(len(ln_cands)):
				if ln_cands[i][0]<ln_0val:
					ln_0val = ln_cands[i][0]
					ln_index = ln_cands[i]
		return ln_index	
	else:
		print "Error in least_norm. Couldn't find positive candidate"
		print index
		print v_shift
		print h_shift
		print rcoeffs
		print fcands
		quit()


signs = np.array([[(1,1),(1,0),(1,-1)],[(0,1),(0,0),(0,-1)],[(-1,1),(-1,0),(-1,-1)]],dtype=int)

def cr_images(index,h_shift,v_shift):
	index_im = []
	for i in range(3):
		index_im.append([])
		for j in range(3):
			rep = tuple(np.array(index)+signs[i,j,0]*h_shift+signs[i,j,1]*v_shift)
			index_im[i].append(rep)
	return index_im



class crystal_site:
	def __init__(self,index):
		# An index (or address)
		self.index = index
		# A list of neighbors
		self.nbors = []
		# A set of local data
		self.local = {}
		# A function used to update local data upon adding neighbors
		# (returns 0, but should also act on itself)
		self.local_update = lambda csite,new_csite: 0 
	def add_nbor(self,new_csite):
		index = new_csite.index
		if index not in self.nbors:
			self.nbors.append(index)
			self.local_update(self,new_csite)
	def add_hex_nbhd(self):
		index_im = np.zeros([6,2],dtype=int)
		index0 = np.array(self.index)
		for i in range(6): 
			index_im[i,:] = np.copy(index0)
		index_im[0,0]+=1
		index_im[1,0]+=1
		index_im[1,1]-=1
		index_im[2,1]-=1
		index_im[3,0]-=1
		index_im[4,0]-=1
		index_im[4,1]+=1
		index_im[5,1]+=1
		self.local['hnbhd'] = []
		for i in range(6):
			self.local['hnbhd'].append(tuple(index_im[i,:]))
		
	
# Class for discrete crystal topology:
class crystal_top:
	def __init__(self):
		self.bulk = {}
		self.bdry = {}
		self.occ = {}
	def add_node_bulk(self,index):
		# Create a neighbor list for each site
		self.bulk[index] = crystal_site(index)
		self.bulk[index].local['type'] = 'bulk'
		self.occ[index] = 1
	def add_node_bdry(self,index):
		self.bdry[index] = crystal_site(index)
		self.bdry[index].local['type'] = 'bdry'
		self.occ[index]=1
		if index not in self.bulk:
			# Include 'index' as a bulk site as well (for now)
			self.bulk[index] = crystal_site(index)
			self.bulk[index].local['type'] = 'bulk'	
	def add_edge(self,index1,index2):
		if index1 in self.bulk and index2 in self.bulk:
			self.bulk[index1].add_nbor(self.bulk[index2])
			self.bulk[index2].add_nbor(self.bulk[index1])
		if index1 in self.bdry and index2 in self.bulk:
			self.bdry[index1].add_nbor(self.bulk[index2])
		if index2 in self.bdry and index1 in self.bulk:
			self.bdry[index2].add_nbor(self.bulk[index1])
		if index1 in self.bdry and index2 in self.bdry:
			self.bdry[index1].add_nbor(self.bdry[index2])
			self.bdry[index2].add_nbor(self.bdry[index1])

	def display(self,gen1,gen2,h_shift,v_shift):
		# Define basic circle:
		circ = circle(1,20)
		pos = {}
		for index in self.bulk:
			new_x = gen1*index[0]+gen2*index[1]
			pos[index]=new_x
		for index in self.bdry:
			pos[index]=gen1*index[0]+gen2*index[1]
		xs = []
		for index in pos:
			xs.append(pos[index])
		xs = np.array(xs)
		# Plot xs
		plt.scatter(xs[:,0],xs[:,1])
		# plot edges 
		for index in self.bulk:
			plt.plot(pos[index][0]+circ[:,0],pos[index][1]+circ[:,1],'g')
			plt.plot(pos[index][0]-h_shift[0]+circ[:,0],pos[index][1]-h_shift[1]+circ[:,1],'--m')
			plt.plot(pos[index][0]+h_shift[0]+circ[:,0],pos[index][1]+h_shift[1]+circ[:,1],'--c')
			plt.plot(pos[index][0]+v_shift[0]+circ[:,0],pos[index][1]+v_shift[1]+circ[:,1],'--y')
			plt.plot(pos[index][0]-v_shift[0]+circ[:,0],pos[index][1]-v_shift[1]+circ[:,1],'--b')
			for j in range(len(self.bulk[index].nbors)):
				index_j = self.bulk[index].nbors[j]
				plt.plot([pos[index][0],pos[index_j][0]],[pos[index][1],pos[index_j][1]],'--k')
		plt.axes().set_aspect('equal','datalim')
		plt.show()

def expand_crystal_topology(cr_top,h_shift,v_shift):
	# Define a dictionary of occupied sites (including images)
	occ = {}
	for index in cr_top.occ:
		occ[index] = cr_top.occ[index]
	# Include images of points in cr_top
	for index in cr_top.occ:
		index_im = cr_images(index,h_shift,v_shift)
		for i in range(3):
			for j in range(3):
				occ[index_im[i][j]]=occ[index]
	# Maintain a list of boundary sites
	bdry = []
	for index in cr_top.bdry:
		bdry.append(index)
	exp_possible = True
	while exp_possible:
		exp_possible= False
		new_bdry = []
		# For each element in bdry, add unvisited neighbors to new_bdry
		for index in bdry:
			# Check this!
			cr_top.bulk[index].add_hex_nbhd()
			for i in range(6):
				index_i = cr_top.bulk[index].local['hnbhd'][i]
				if not index_i in occ:
					exp_possible = True
					# Update the list of occupied sites:
					index_i_im = cr_images(index_i,h_shift,v_shift)
					for ii in range(3):
						for jj in range(3):
							occ[index_i_im[ii][jj]]=1
					# Add index_i to cr_top.bulk
					cr_top.add_node_bulk(index_i)
					# Add index_i to new_bdry
					new_bdry.append(index_i)
		# Update bdry:
		bdry = new_bdry
	# Define edges
	for index in cr_top.bulk:
		cr_top.bulk[index].add_hex_nbhd()
		for j in range(6):
			index_j = cr_top.bulk[index].local['hnbhd'][j]
			if index_j in cr_top.occ:
				if cr_top.occ[index_j]==1:
					cr_top.add_edge(index,index_j)
	# Define the boundary component
	for index in cr_top.bulk:
		if len(cr_top.bulk[index].nbors)<6:
			cr_top.add_node_bdry(index)

					
def self_intersection(cr_top,h_shift,v_shift,dist=1):
	# Check if ls_seq intersects with any of its nearest images
	offsets = []
	ext = 2*dist+1
	for i in range(ext):
		factor_h = i-dist
		for j in range(ext):
			factor_v = j-dist
			offsets.append(factor_h*h_shift+factor_v*v_shift)
	for i in range(len(offsets)):
		if offsets[i][0]==0 and offsets[i][1]==0:
			continue
		for index in cr_top.bulk:
			arr_index = np.array(index)+offsets[i]
			if tuple(arr_index) in cr_top.bulk:
				return True
	return False

ls_incr = [np.array([1,0],dtype=int),np.array([1,-1],dtype=int),np.array([0,1],dtype=int)]

def crystal_top_test(R,test_seq_):
	# Note: Extend test_seq_ if it does not accommodate the underlying crystal shape
	ex_list = get_closed_crystals(R,curve_list)
	for j in range(len(ex_list)):
		# Obtain standard vectors and loop indices for orientation
		std_vecs = standard_vectors(ex_list[j][0],[ex_list[j][1],ex_list[j][2]])
		if std_vecs[4]==0 or std_vecs[5]==0:
			continue
		[ua,u1,u2,u3,nua,nu2,nu3] = std_vecs[:7]
		print [nu2,nua]	
		h_shift= np.array([np.sign(nu2)*nua,abs(nu2)+1],dtype=int)
		h_shift_v = h_shift[0]*ua+(h_shift[1]-1)*u2+np.sign(nu2*nu3)*u3
		self_intersect = True
		test_seq__=test_seq_
		v_shift_v = np.zeros(2)
		v_shift = np.zeros(2,dtype=int)
		for i in range(len(test_seq__)):
			v_shift+=ls_incr[test_seq__[i]]
			v_shift_v+=ls_incr[test_seq__[i]][0]*ua+ls_incr[test_seq__[i]][1]*u2 
		# Check if v_shift requires augmentation
		# Look for self-intersections
		test_seq = test_seq__
		# Define the crystal domain and topology that will accommodate the wrapped line slip defect
		# Add the 'left' and 'right' sides of the line slip to crystal_top
		# (also, compute the index shift associated with the line slip configuration)
		cr_top = crystal_top()
		node = np.zeros(2,dtype=int)
		node_seq = []
		pnode = np.zeros(2,dtype=int)
		nnode = np.copy(h_shift)
		nnode[1]-=1
		nnode_seq = [tuple(nnode)]
		nnode_set = 1
		for i in range(len(test_seq)):
			node_seq.append(tuple(node))
			# Add 'node' to cr_top:
			cr_top.add_node_bdry(tuple(node))
			# Add neighbors of 'node' from across the boundary to cr_top
			if test_seq[i]==0:
				if test_seq[i-1]<2:
					nnode[0]+=1
					nnode_seq.append(tuple(nnode))
					nnode_set+=1
				else:
					# keep nnode and nnode_seq unchanged
					pass
			elif test_seq[i]==1:
				if test_seq[i-1]<2 and test_seq[i-2]<2:
					# Set value of 'occ' to 0 at current site:
					cr_top.occ[tuple(nnode)]=0
					# modify nnode_seq[i-1] 
					nnode[1]-=1
					nnode_seq[nnode_set-1]=tuple(nnode)
					nnode[0]+=1
					nnode_seq.append(tuple(nnode))
					nnode_set+=1
				else:
					if test_seq[i-1]==2:
						# This case should not occur (it is contractible)
						print "Warning: line slip configuration is contractible"
						# Remove the last element of nnode_seq
						nnode_set-=1
						nnode[1]-=1
					if test_seq[i-2]==2:
						# A sheared vacancy is formed along the line slip.
						cr_top.occ[tuple(nnode)]=0
						nnode[0]+=1
						nnode[1]-=1
						nnode_seq[nnode_set-1]=tuple(nnode)
			elif test_seq[i]==2:
				if test_seq[i-1]<2:
					nnode[0]+=1
					# Add the vertical bank
					nnode_seq.append(tuple(nnode))
					nnode[1]+=1
					# Add horizontal bank
					nnode_seq.append(tuple(nnode))
					nnode_set+=2
				else:
					nnode[1]+=1
					nnode_seq.append(tuple(nnode))
					nnode_set+=1
			# increment 'node' according to test_seq:
			node+=ls_incr[test_seq[i]]
			cr_top.add_edge(tuple(pnode),tuple(node))
			pnode = np.copy(node)
		#cr_top.add_edge((0,0),tuple(node_seq[-1]))
		#cr_top.display(ua,u2,h_shift_v,v_shift_v)
		self_intersect = self_intersection(cr_top,h_shift,v_shift)
		if self_intersect:
			print "(Skipping)"
			continue
		for i in range(len(node_seq)):
			cr_top.add_edge(tuple(node_seq[i-1]),tuple(node_seq[i]))
		# Add edges along nnode_seq
		for i in range(len(nnode_seq)-1):
			cr_top.add_node_bdry(tuple(nnode_seq[i]))
		for i in range(1,len(nnode_seq)-1):
			cr_top.add_edge(tuple(nnode_seq[i]),tuple(nnode_seq[i-1]))
		cr_top.add_edge(tuple(nnode_seq[0]),tuple(nnode_seq[-2]))
		# Display current state of cr_top:
		cr_top.display(ua,u2,h_shift_v,v_shift_v)
		print "Extending to full domain"	
		# Expand cr_top to full interior domain, accounting for periodicity
		expand_crystal_topology(cr_top,h_shift,v_shift)
		cr_top.display(ua,u2,h_shift_v,v_shift_v)

#test_seq = [0,0,0,1,0,1,0,2,2,0,2,0]
#crystal_top_test(3.1,test_seq)

def generate_ls_crystal_top(std_vecs,test_seq_):
	# Obtain standard vectors and loop indices for orientation
	if std_vecs[4]==0 or std_vecs[5]==0:
		return [None,False]
	[ua,u1,u2,u3,nua,nu2,nu3] = std_vecs[:7]
	print [nu2,nua]	
	h_shift= np.array([np.sign(nu2)*nua,abs(nu2)+1],dtype=int)
	h_shift_v = h_shift[0]*ua+(h_shift[1]-1)*u2+np.sign(nu2*nu3)*u3
	self_intersect = True
	test_seq__=test_seq_
	v_shift_v = np.zeros(2)
	v_shift = np.zeros(2,dtype=int)
	for i in range(len(test_seq__)):
		v_shift+=ls_incr[test_seq__[i]]
		v_shift_v+=ls_incr[test_seq__[i]][0]*ua+ls_incr[test_seq__[i]][1]*u2 
	# Check if v_shift requires augmentation
	# Look for self-intersections
	test_seq = test_seq__
	# Define the crystal domain and topology that will accommodate the wrapped line slip defect
	# Add the 'left' and 'right' sides of the line slip to crystal_top
	# (also, compute the index shift associated with the line slip configuration)
	cr_top = crystal_top()
	node = np.zeros(2,dtype=int)
	node_seq = []
	pnode = np.zeros(2,dtype=int)
	nnode = np.copy(h_shift)
	nnode[1]-=1
	nnode_seq = [tuple(nnode)]
	nnode_set = 1
	for i in range(len(test_seq)):
		node_seq.append(tuple(node))
		# Add 'node' to cr_top:
		cr_top.add_node_bdry(tuple(node))
		# Add neighbors of 'node' from across the boundary to cr_top
		if test_seq[i]==0:
			if test_seq[i-1]<2:
				nnode[0]+=1
				nnode_seq.append(tuple(nnode))
				nnode_set+=1
			else:
				# keep nnode and nnode_seq unchanged
				pass
		elif test_seq[i]==1:
			if test_seq[i-1]<2 and test_seq[i-2]<2:
				# Set value of 'occ' to 0 at current site:
				cr_top.occ[tuple(nnode)]=0
				# modify nnode_seq[i-1] 
				nnode[1]-=1
				nnode_seq[nnode_set-1]=tuple(nnode)
				nnode[0]+=1
				nnode_seq.append(tuple(nnode))
				nnode_set+=1
			else:
				if test_seq[i-1]==2:
					# This case should not occur (it is contractible)
					print "Warning: line slip configuration is contractible"
					nnode[1]-=1
				if test_seq[i-2]==2:
					# A sheared vacancy is formed along the line slip.
					cr_top.occ[tuple(nnode)]=0
					nnode[0]+=1
					nnode[1]-=1
					nnode_seq[nnode_set-1]=tuple(nnode)
		elif test_seq[i]==2:
			if test_seq[i-1]<2:
				nnode[0]+=1
				# Add the vertical bank
				nnode_seq.append(tuple(nnode))
				nnode[1]+=1
				# Add horizontal bank
				nnode_seq.append(tuple(nnode))
				nnode_set+=2
			else:
				nnode[1]+=1
				nnode_seq.append(tuple(nnode))
				nnode_set+=1
		# increment 'node' according to test_seq:
		node+=ls_incr[test_seq[i]]
		cr_top.add_edge(tuple(pnode),tuple(node))
		pnode = np.copy(node)
	#cr_top.add_edge((0,0),tuple(node_seq[-1]))
	#cr_top.display(ua,u2,h_shift_v,v_shift_v)
	self_intersect = self_intersection(cr_top,h_shift,v_shift)
	if self_intersect:
		success_flag=False
		return [cr_top,success_flag]
	for i in range(len(node_seq)):
		cr_top.add_edge(tuple(node_seq[i-1]),tuple(node_seq[i]))
	# Add edges along nnode_seq
	for i in range(len(nnode_seq)-1):
		cr_top.add_node_bdry(tuple(nnode_seq[i]))
	for i in range(1,len(nnode_seq)-1):
		cr_top.add_edge(tuple(nnode_seq[i]),tuple(nnode_seq[i-1]))
	cr_top.add_edge(tuple(nnode_seq[0]),tuple(nnode_seq[-2]))
	# Display current state of cr_top:
	# Expand cr_top to full interior domain, accounting for periodicity
	expand_crystal_topology(cr_top,h_shift,v_shift)
	success_flag=True
	return [cr_top,success_flag]

def display3D(idxs):
	xs = np.zeros([len(idxs),3])
	#for i in range(len(idxs))
#		xs[i,0] = idxs[i]['x']
#		xs[i,1] = idxs[i]['y']
#		xs[i,2] = idxs[i]['z']
	plt.scatter(idxs[:]['x'],idxs[:]['z'])
	plt.show() 

# NOTE: fix this!
def write_init_conditions(curve_index,Nsamples,ls_state,pars):
	stretch = pars[0]
	# Additional simulation parameters
	sim_pars = pars[1]
	depl_str = sim_pars[0]
	dpr = sim_pars[1]
	pdiffuse= sim_pars[2]
	rho_eq = sim_pars[3]
	Tsteps = sim_pars[4]
	dt = sim_pars[5]
	# Directory names
	dirname = pars[2]
	# (including data folders if necessary)
	mkdir_s("data")
	mkdir_s("data/parameters")
	mkdir_s("data/"+dirname)
	mkdir_s("data/parameters/"+dirname)
	curve = curve_list[curve_index]
	# Define a list of [v2d, na, nb] for each sample point on the given curve
	samples = exact_v2d_samples(curve,Nsamples)
	for i in range(len(samples)):
		std_vecs = standard_vectors(samples[i][1],samples[i][2])
		[cr_top,success_flag] = generate_ls_crystal_top(std_vecs,ls_state)
		if success_flag:
			[ua,u1,u2,u3,nua,nu2,nu3] = std_vecs[:7]
			# Write configuration to a coordinate file (incorporating the stretch factor)
			# Also write parameter file
			idxs = np.zeros(len(cr_top.bulk),dtype=idx)
			R_true = samples[i][0]*stretch
			R_cyl = R_true-1
			# Compute the length and twist of the cylinder from ls_seq and ua,u2
			v_shift_v = np.zeros(2)
			for ii in range(len(ls_state)):
				lat_disp = ls_incr[ls_state[ii]]
				v_shift_v+=lat_disp[0]*ua+lat_disp[1]*u2
			L_true = stretch*abs(v_shift_v[1])
			twist = np.sign(v_shift_v[1])*stretch*v_shift_v[0]/R_true
			cos_twist = np.cos(twist)
			sine_twist = np.sin(twist)
			ii = 0
			# Display cr_top: cr_top.display(ua,u2,h_shift_v,v_shift_v)
			for index in cr_top.bulk:
				x2d = index[0]*ua+index[1]*u2
				x2d*=stretch
				idxs[ii]['id'] = ii
				idxs[ii]['x'] = R_true*np.cos(x2d[0]/R_true)
				idxs[ii]['y'] = R_true*np.sin(x2d[0]/R_true)
				idxs[ii]['z'] = x2d[1]-0.5*L_true*np.sign(v_shift_v[1])
				if idxs[ii]['z']>0.5*L_true:
					# shift down, with counter rotation
					idxs[ii]['z']-=L_true
					old_x=idxs[ii]['x']
					old_y=idxs[ii]['y']
					# CHECK THIS! (and the opposite case as well)
					idxs[ii]['x'] = cos_twist*old_x+sine_twist*old_y
					idxs[ii]['y'] = cos_twist*old_y-sine_twist*old_x
				elif idxs[ii]['z']<-0.5*L_true:
					# shift up, with rotation (CHECK THIS!)
					idxs[ii]['z']+=L_true
					old_x=idxs[ii]['x']
					old_y=idxs[ii]['y']
					idxs[ii]['x'] = cos_twist*old_x-sine_twist*old_y
					idxs[ii]['y'] = cos_twist*old_y+sine_twist*old_x
				ii+=1
			# Display the 3D coords:
			display3D(idxs)
			# Make directory folder for current trial
			mkdir_s("data/"+dirname+"/trial"+str(i))
			# Save coordinate file
			np.savetxt("data/"+dirname+"/trial"+str(i)+"/f.0",idxs,"%d %f %f %f")
			# Save parameters (first compute the height and twist from ls_seq and std_vecs)
			pars = np.zeros(1,dtype=par_dtype_sc)
			pars['Np']=-1 
			pars['R'] = R_cyl 
			pars['L'] = 0.5*L_true 
			pars['depl_str'] = depl_str
			pars['dpr'] = dpr
			pars['pdiffuse'] = pdiffuse
			pars['rho_eq'] = rho_eq
			pars['Tsteps'] = Tsteps
			pars['dt'] = dt
			pars['twist'] = twist
			np.savetxt("data/parameters/"+dirname+"/trial"+str(i)+".par",pars,"%d\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f")
											   # 1  2   3   4   5   6   7   8   9   10 

# 8.28.2018
# Generate a directory name from a ls_seq:
def make_dirname(ls_seq):
	dirname = ""
	for i in range(len(ls_seq)):
		dirname+=str(ls_seq[i])
	return dirname

# 8.28.2018
# Generate a single fractional vacancy configuration
def generate_fv_init_conditions(curve_index,mode,Nsamples,sim_pars,extra=10):
	# Define a line slip corresponding to mode:
	ls_seq = []
	pre = extra/2
	for i in range(pre):
		ls_seq.append(0)
	if mode==1:
		ls_seq.append(1)
	elif mode==2:
		ls_seq.append(2)
	else:
		print "Error: mode must be 1 (fast) or 2 (slow)"
		quit()
	for i in range(extra-pre):
		ls_seq.append(0)
	# Generate initial conditions for this line slip configuration
	dirname = "curve_"+str(curve_index)+"_sfv_t"+str(mode)+"_ext"+str(extra)
	pars = [1.02,sim_pars,dirname]
# write_init_conditions(curve_index,Nsamples,ls_state,pars)
	write_init_conditions(curve_index,Nsamples,ls_seq,pars)

# 8.28.2018
# Test generate_fv_init_conditions (check curve indices first)
#generate_fv_init_conditions()

# Generate points on a lattice with fractional vacancies for v2d at radius R:
def fractional_vacancies(ex_entry,R,L,mode,dirname,trinum):
	v2d = np.copy(ex_entry[0])
	[na,nb] = ex_entry[1:3]
	C = 2*np.pi*R
	[ua,u1,u2,u3]=standard_vectors(v2d)
	[x_list,x_array,occ]=generate_pts(v2d,na,nb,R,-L,L)
	(Nr,Nc)=x_array.shape[:2]
	if min(Nr,Nc)<2:
		return
	# x[:,0,:] contains one side of the line slip defect. 
	v1 = x_array[1,0,:]-x_array[0,0,:] # generator of translational symmetry unbroken by line slip
	v2 = x_array[0,1,:]-x_array[0,0,:] # transverse, aligned with l.s.g.
	
	if np.dot(v1,v2)>0:
		v3 = v2-v1 # transverse, obtuse angle with l.s.g.
	else:
		v3 = v2+v1
		aux_v = v2
		v2 = v3
		v3 = aux_v
	# Check if the frame is properly aligned with the line slip transverse:
	vls = v2d[4:]
	if np.dot(v3,vls)>0:
		vls*=-1
	# Place the f.v. at the center of the frame:
	# Create a 'fast' or slow f.v. by translating one side of the line slip defect or the other
	mlist = []
	for j in range(Nr):
		if occ[j,0]:
			mlist.append(x_array[j,0,:])
	if mode==0:
		# Fast case:
		if np.dot(v1,vls)>0:
			# Define the offset in terms of v1,v2,v3:
			offset_v = vls+v3
			sign = -1
		else:
			offset_v = vls+v2
			sign = 1
	else:
		# Slow case:
		if np.dot(v1,vls)>0:
			offset_v = vls+v2
			sign = 1
		else:
			offset_v = vls+v3
			sign = -1
	#print [mode,np.linalg.norm(x3d(offset_v,R))]
	# Find the first row at which the z-coordinate of the 
	# first entry changes sign
	x_ref = x_array[0,0,:]
	sign0 = np.sign(x_ref[1])
	for j in range(Nr):
		if np.sign(x_array[j,0,1])!=sign0:
			x_ref = x_array[j,0,:]
			break
	for j in range(Nr):
		if sign*np.dot(x_array[j,0,:]-x_ref,v1)>0:
			x_array[j,0,:]+=offset_v
	# transcribe positions to a new array:
	new_xs = []
	for j in range(Nr):
		for jj in range(Nc):
			if occ[j,jj]:
				nx = np.copy(x_array[j,jj,:])
				new_xs.append(nx)
	new_xs = np.array(new_xs)
	mlist = np.array(mlist)
	coords = np.zeros(len(new_xs),dtype=idx)
	for jj in range(len(new_xs)):
		coords[jj][0] = jj
		x = x3d(new_xs[jj,:],R)
		coords[jj][1] = x[0]
		coords[jj][2] = x[1]
		coords[jj][3] = x[2]
	np.savetxt("sweep/"+dirname+"/trial"+str(trinum)+"mode"+str(mode)+"/f.0",coords,fmt="%d %f %f %f")
	return new_xs
	
def test_fvgen(Ntrials):
	print "Loading data for interpolation"
	t0 = time.time()
	[Rsamples,tsamples,v2Ds,curve_features,curve_featuresc]=load_interp_data()
	t1 = time.time()
	time_lid = t1-t0
	t0 = time.time()
	print "Generating approximate curves for guesses"
	curve_list = make_guess_curves(Rsamples,tsamples,v2Ds,curve_features,curve_featuresc)
	t1 = time.time()
	time_mgc = t1-t0
	print "Runtime: load_interp_data, make_guess_curves"
	print [time_lid,time_mgc]
	for j in range(Ntrials):
		R = 2.5+2*np.random.rand()
		C = 2*np.pi*R
		L = 20
		# Get lattice vectors:
		t0 = time.time()
		print "Getting precise estimates for closed crystals at radius "+str(R)
		ex_list = get_closed_crystals(R,curve_list)
		t1 = time.time()
		time_gcc = t1-t0
		# Try a random entry in fv:
		print "Generating fractional vacancies for a random orientation"
		exi = np.random.randint(len(ex_list))
		t0 = time.time()
		fractional_vacancies(ex_list[exi],R,L,0)
		display_fv_pts(x_fast,C)
		fractional_vacancies(ex_list[exi],R,L,1)
		display_fv_pts(x_slow,C)
		t1 = time.time()
		time_fv = t1-t0
		print "Runtimes: get_closed_crystals(), fractional_vacancies()"
		print [time_gcc,time_fv]

def mkdir_s(dirname):
	try:
		os.mkdir(dirname)
	except:
		pass

def rm_s(name):
	try:
		os.execl("rm -r "+name)
	except:
		pass

# Define the data type for parameters:
par_dtype = np.dtype([('Np','i4'),('R','f4'),('L','f4'),('depl_str','f4'),('dpr','f4'),('pdiffuse','f4'),('rho_eq','f4'),('pbc_flag','i4'),('Tsteps','f4'),('dt','f4')])
# Parameters for 'single crystal' simulations (same as above, except with a 'twist' parameter at the end)
par_dtype_sc = np.dtype([('Np','i4'),('R','f4'),('L','f4'),('depl_str','f4'),('dpr','f4'),('pdiffuse','f4'),('rho_eq','f4'),('Tsteps','f4'),('dt','f4'),('twist','f4')])

def init_fv_configs(dirname,pars,theta_lims=None):
	R = pars[0][1]
	L = pars[0][2]
	# Make main data directory:
	mkdir_s("sweep/"+dirname)
	# Make a parameter directory:
	mkdir_s("sweep/parameters")
	mkdir_s("sweep/parameters/"+dirname)
	#print "Loading data for interpolation"
	#[Rsamples,tsamples,v2Ds,curve_features,curve_featuresc]=load_interp_data()
	#print "Generating approximate curves for guesses"
	#curve_list = make_guess_curves(Rsamples,tsamples,v2Ds,curve_features,curve_featuresc)
	# Get lattice vectors:
	print "Getting precise estimates for closed crystals at radius "+str(R)
	ex_list = get_closed_crystals(R,curve_list)
	print "Generating initial conditions for each type of fractional vacancy"
	if theta_lims!=None:
		# Create initial conditions for line slip generators within angle limits
		print "Interval:"+str(theta_lims)
	else:
		# Create initial conditions for every line slip angle
		proceed_flag=True
	trinum = 0
	for j in range(len(ex_list)):
		# Compute the angle made by va
		v2d = ex_list[j][0]
		theta = np.arccos(v2d[1]/2.) # Assumes radius 1
		if theta_lims!=None:
			if theta>min(theta_lims) and theta<max(theta_lims):
				proceed_flag=True
			else:
				proceed_flag=False
		if proceed_flag:
			# Generate a parameter file, create a trial directory, and create initial conditions:
			pfname = "sweep/parameters/"+dirname+"/trial_"+str(trinum)+".par"
			np.savetxt(pfname,pars,"%d\n%f\n%f\n%f\n%f\n%f\n%f\n%d\n%f\n%f")
			mkdir_s("sweep/"+dirname+"/trial"+str(trinum)+"mode0")
			# Generate initial coordinates
			fractional_vacancies(ex_list[j],R,L,0,dirname,trinum)
			mkdir_s("sweep/"+dirname+"/trial"+str(trinum)+"mode1")
			fractional_vacancies(ex_list[j],R,L,1,dirname,trinum)
			trinum+=1

def compute_gap_deficit(v2d,R):
	std_vecs = standard_vectors(v2d)
	[ua,u1,u2,u3]=std_vecs
	gap2 = 0.5*(u2-u3)
	gap1 = 0.5*(u3+u1)
	gap2_3d = np.array([R*np.cos(gap2[0]/R)-R,R*np.sin(gap2[0]/R),gap2[1]])
	gap1_3d = np.array([R*np.cos(gap1[0]/R)-R,R*np.sin(gap1[0]/R),gap1[1]])
	return np.dot(gap2_3d,gap2_3d)-np.dot(gap1_3d,gap1_3d)

# Functions for mobility of fractional vacancies:	
def find_equiv_radii():
	R_equiv = {}
	N_not_found = 0
	Req_not_found = []
	for i in range(len(curve_list)):
		[na,nb,ocase,icase] = curve_list[i][:4]
		case = [ocase,icase]
		# Solve for the radius or angle at which the 'gaps' associated with the fast and slow kinks are equivalent
		if curve_list[i][4]==0:
			R1 = curve_list[i][6]
			R2 = curve_list[i][7]
			v2d_of_R = lambda R: scipy.optimize.root(ls_cconstraint,curve_list[i][5](R),(na,nb,R,case)).x
			f = lambda R: compute_gap_deficit(v2d_of_R(R),R)
			R_equiv_r = scipy.optimize.root(f,0.5*(R1+R2))
			if R_equiv_r.success and min(R1,R2)<R_equiv_r.x<max(R1,R2):
				R_equiv[i] = R_equiv_r.x
			else: 
				Req_not_found.append(i)
				N_not_found+=1
		elif curve_list[i][4]==1:
			t1 = curve_list[i][6]
			t2 = curve_list[i][7]
			case = curve_list[i][2:4]
			v2dR_of_theta = lambda theta: scipy.optimize.root(ls_cconstraint_theta,curve_list[i][5](theta),(na,nb,theta,icase)).x
			v2d_of_theta = lambda theta: get_v2d_of_v2dR(v2dR_of_theta(theta),theta,case)
			f = lambda theta: compute_gap_deficit(v2d_of_theta(theta),v2dR_of_theta(theta)[4])
			theta_equiv = scipy.optimize.root(f,0.5*(t1+t2))
			if theta_equiv.success and min(t1,t2)<theta_equiv.x<max(t1,t2):
				R_equiv[i] = [v2dR_of_theta(theta_equiv.x)[4],theta_equiv.x]
			else: 
				Req_not_found.append(i)
				N_not_found+=1
	print "Unable to find equivalent radii for "+str(N_not_found)+" curves (check 'compute_gap_deficit')"
	print Req_not_found
	return R_equiv

#print "Finding radii at which fast and slow fractional vacancies match"
#R_equiv = find_equiv_radii()
#print "Done"

# Test the symmetrized product of lattice displacements as a classifier for orientation
def sym_product_test(R_interval,Nsamples=100):
	Rsamples = np.linspace(R_interval[0],R_interval[1],Nsamples)
	min_sp_seps = 1e99*np.ones(Nsamples)
	for j in range(Nsamples):
		print "Processing "+nth(j)+" sample"
		sp_samples = []
		ex_list = get_closed_crystals(Rsamples[j],curve_list)
		if len(ex_list)>0:
			for jj in range(len(ex_list)):
				v2djj= ex_list[jj][0]
				va = v2djj[0:2]
				vb = v2djj[2:4]
				if np.dot(va,vb)<0:
					vc = vb-va
				else:
					vc=vb+va
				vd=-va
				ve=-vb
				vf=-vc
				lattice_disps= np.array([va,vb,vc,vd,ve,vf])
				sp=symmetric_product(lattice_disps)
				sp_samples.append(sp.flatten())
				for jjj in range(jj):
					disp = sp_samples[jj]-sp_samples[jjj]
					sep = np.dot(disp,disp)
					if sep<min_sp_seps[j]:
						min_sp_seps[j]=sep
			min_sp_seps[j]=np.sqrt(min_sp_seps[j])
		else:
			min_sp_seps[j]=0
	plt.plot(Rsamples,min_sp_seps)
	plt.scatter(Rsamples,min_sp_seps,c='r',s=1)
	plt.show()


