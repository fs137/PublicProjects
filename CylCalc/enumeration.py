# AUTHOR: HENRY WILKIN
from common import *
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize
import scipy.interpolate
from algebra import *
from common import *


# Enumerate commensurate crystal geometries:
mkdir_s("LSDATA")
epsilon = 1e-6
twopi = 2*np.pi

def interval_overlap(ival1,ival2):
	return (ival2[0]<ival1[0]<ival2[1]) or (ival2[0]<ival1[1]<ival2[1])

def box_overlap(box1,box2):
	if len(box1)==len(box2):
		overlap = True
		for j in range(len(box1)):
			overlap=overlap and interval_overlap(box1[j],box2[j])
		return overlap
	else:
		print "Error (box_overlap): box dimensions do not agree"
		return

def zero_function(x):
	return np.zeros(6)

# Solve for R, theta given n0, n5
def flat_cconstraint(invRtheta,n0,n5):
	theta = invRtheta[1]
	invR = invRtheta[0]
	rval = np.zeros(2)
	if abs(np.cos(theta)*invR)<1 and abs(np.cos(theta-np.pi/3)*invR)<1:
		rval[0] = n0*np.arcsin(np.cos(theta)*invR)+n5*np.arcsin(np.cos(theta-np.pi/3)*invR)-np.pi
	else:
		rval[0] = (invR-0.05)+0.01*np.random.rand()
	rval[1] = n0*np.sin(theta)+n5*np.sin(theta-np.pi/3)
	return rval

def curved_cconstraint(vR,n0,n5):
	invR = 1./vR[6]
	hinvR = 0.5*invR
	# length constraints: x2
	# cylinder constraints: x2
	# difference length constraint: x1
	# closure: x2
	# total dof: 7 = 2+2+2+1
	rval = np.zeros(7)
	rval[0]= np.dot(vR[:3],vR[:3])-4
	rval[1]= np.dot(vR[3:6],vR[3:6])-4
	diff = vR[:3]-vR[3:6]
	rval[2] = np.dot(diff,diff)-4
	rval[3] = 2*vR[0]*vR[6]+vR[0]*vR[0]+vR[1]*vR[1]
	rval[4] = 2*vR[3]*vR[6]+vR[3]*vR[3]+vR[4]*vR[4]
	projxy0 = np.sqrt(np.dot(vR[:2],vR[:2]))
	projxy5 = np.sqrt(np.dot(vR[3:5],vR[3:5]))
	if abs(projxy0*hinvR)<=1 and abs(projxy5*hinvR)<=1:
		rval[5] = n0*np.arcsin(projxy0*hinvR)+n5*np.arcsin(projxy5*hinvR)-np.pi
	elif abs(invR*vR[1])<=1 and abs(invR*vR[4])<=1:
		dtheta0 = np.arcsin(vR[1]*invR)
		dtheta5 = np.arcsin(vR[4]*invR)
		if vR[0]<-vR[6]:
			dtheta0 = 0.5*np.pi-dtheta0
		if vR[3]<-vR[6]:
			dtheta5 = 0.5*np.pi-dtheta5
		rval[5] = n0*dtheta0+n5*dtheta5-2*np.pi
		#if invR*vR[0]>=-1 and invR*vR[3]>=-1:
		#	rval[5] = n0*np.arcsin(vR[1]*invR)+n5*np.arcsin(vR[4]*invR)-2*np.pi
		#else:
		#	rval[5] = n0*np.arccos(1+vR[0]*invR)+n5*np.arccos(1+vR[3]*invR)
	elif abs(invR*vR[0]+1)<=1 and (invR*vR[3]+1)<=1:
		dtheta0 = np.arccos(1-abs(invR*vR[0]))
		dtheta5 = np.arccos(1-abs(vR[3]*invR))
		rval[5] = n0*dtheta0+n5*dtheta5-2*np.pi
	else:
		rval[5] = 0.1+invR 
	if invR<0:
		rval[5] = 0.1*vR[6]-3 # Chosen to return to positive R values
	rval[6] = n0*vR[2]+n5*vR[5]
	return rval

def get_v2d_of_v2dR(v2dR,theta,case):
	v2d = np.zeros(6)
	v2d[:2] = np.copy(v2dR[:2])
	v2d[4:6] = np.copy(v2dR[2:4])
	v0 = np.zeros(2)
	v0[0] = v2dR[4]*2*np.arccos(np.cos(theta)/v2dR[4])
	v0[1] = 2*np.sin(theta)
	#v0 = np.array([v2dR[4]*2*np.arccos(np.cos(theta)/v2dR[4]),2*np.sin(theta)])
	#print v2dR.shape
	#print v2d.shape
	#print v0.shape
	if case==[2,0]:
		# va = v1, vb = v0
		v2d[2:4] = np.copy(v0)
	elif case==[2,1]:
	 	# va = v1, vb = v5
		v2d[2:4] = v0-v2d[:2]
	return v2d

max_attempts = 5e3
# v0, v5, R in 3D: start with a simplified (flat) 2D model, then expand to 3D
def enumerate_cc(n_max,Rmax=None):
	# Determine if n_max is sufficient to determine all commensurate crystals up to R_max:
	if Rmax!=None:
		n_max = max(n_max,int(2*np.pi/(2*np.arcsin(1./Rmax)))+1)
	sol_list = []
	indices = []
	failed_indices = []
	succ_array = -np.ones([n_max,n_max],dtype=int)
	sol_ht = {}
	sol_array = -2*np.ones([n_max,n_max,7]) # Array containing v0,v5, and R for each commensurate solution
	for n0 in range(n_max):
		for n5 in range(n_max):
			if(n0+n5==0):
				continue
			to_be_found = True
			n_attempts = 0
			while to_be_found and n_attempts<max_attempts:
				n_attempts+=1
				invRtheta0 = np.random.rand(2)
				invRtheta0[0]*=0.1
				if n5==0:
					invRtheta0[1]=0.0
				elif n0==0:
					invRtheta0[1]=np.pi/3
				else:
					invRtheta0[1]=0.1+invRtheta0[1]
				invRtheta_root = scipy.optimize.root(flat_cconstraint,invRtheta0,(n0,n5),tol=1e-7)
				to_be_found = not invRtheta_root.success
			
			if not to_be_found:
				n_attempts=0
				invRtheta1 = invRtheta_root.x
				invRthetaval = np.dot(invRtheta_root.fun,invRtheta_root.fun)
				if invRthetaval>1e-6:
					#print "Incorrect flat root "+str([n0,n5])+" err:"+str(invRthetaval)
					continue
				theta = invRtheta1[1]
				R = 1./invRtheta1[0]
				# Use invRtheta1 as an initial guess for the exact constraint:
				to_be_found = True
				while n_attempts<max_attempts and to_be_found:
					n_attempts+=1
					vR0 = np.zeros(7)
					vR0[0] = 0
					vR0[1] = 2*np.cos(theta)
					vR0[2] = 2*np.sin(theta)
					vR0[3] = 0
					vR0[4] = 2*np.cos(theta-np.pi/3)
					vR0[5] = 2*np.sin(theta-np.pi/3)
					vR0[6] = R
					vR0[:6]+=np.random.rand(6)*0.1
					vR1root = scipy.optimize.root(curved_cconstraint,vR0,(n0,n5),tol=1e-6)
					to_be_found = not vR1root.success
				if not to_be_found:
					vR1 = vR1root.x
					if np.dot(vR1root.fun,vR1root.fun)<1e-6:
						# Add vR1 to a list of solutions
						if Rmax==None:
							pass
						else:
							if vR1[6]>Rmax:
								succ_array[n0,n5]=-3
								continue
						sol_list.append(vR1)
						sol_array[n0,n5,:] = vR1
						indices.append([n0,n5])
						succ_array[n0,n5]=len(sol_list)-1
					else:
						succ_array[n0,n5]=-2
						failed_indices.append([n0,n5])
						#print "Failed to converge at n0,n5="+str([n0,n5])
				else:
					#print "Failed to converge at n0,n5="+str([n0,n5])
					succ_array[n0,n5]=-2
					failed_indices.append([n0,n5])
			else:
				#print "Failed to converge at n0,n5="+str([n0,n5])+" (flat)"
				failed_indices.append([n0,n5])
				succ_array[n0,n5]=-2
	# Loop over failed indices until succ_array no longer changes
	print "Refining solution list"
	changed = True
	iter = 0
	Nfailed = len(failed_indices)
	while(changed):
		#print iter
		iter+=1
		changed=False
		i = 0
		while i < Nfailed:
			#print "Verifying "+str(failed_indices[i])+" of "+str(Nfailed)
			# Check if failed_indices[i] is adjacent to a successful index in succ_array:
			[n0,n5] = failed_indices[i]
			starting_guesses = []	
			for ii in range(max(n0-1,0),min(n0+1,n_max)):
				if succ_array[ii,n5]>-1:
					starting_guesses.append(succ_array[ii,n5])
			for ii in range(max(n5-1,0),min(n5+1,n_max)):
				if succ_array[n0,ii]>-1:
					starting_guesses.append(succ_array[n0,ii])
			if len(starting_guesses)>0:
				#print "Nearby solutions: "+str(starting_guesses)
				# Run max_attempts random guesses near each starting guess:
				n_attempts = 0
				not_yet_found = True
				while not_yet_found and n_attempts<1e2:
					n_attempts+=1
					sgi = np.random.randint(len(starting_guesses))
					vR0 = sol_list[starting_guesses[sgi]]+0.01*np.random.rand(7)
					vR1root = scipy.optimize.root(curved_cconstraint,vR0,(n0,n5),tol=1e-6)
					if vR1root.success:
						changed=True
						not_yet_found = False
						vR1 = vR1root.x 
						sol_list.append(vR1)
						sol_array[n0,n5,:] = np.copy(vR1)
						indices.append([n0,n5])
						print "Added "+str([n0,n5])
						failed_indices[i] = failed_indices[len(failed_indices)-1]
						Nfailed-=1
						succ_array[n0,n5]=len(sol_list)-1
				if not_yet_found:
					#print "still failed to converge at "+str(failed_indices[i])
					i+=1
			else:
				i+=1
				
	# (End of loop 'while changed')				
	# Visualize solutions:
	Rthetas = []
	for i in range(len(sol_list)):
		# check each solution:
		check = curved_cconstraint(sol_list[i],indices[i][0],indices[i][1])
		if np.dot(check,check)<0.00001:
			Rthetas.append([sol_list[i][6],np.arcsin(0.5*sol_list[i][2])])
		else:
			print "Error at n0,n5="+str(indices[i])+". cc value = "+str(np.dot(check,check))
	Rthetas = np.array(Rthetas)
	plt.scatter(Rthetas[:,0],Rthetas[:,1])
	for i in range(len(Rthetas)):
		plt.text(Rthetas[i,0],Rthetas[i,1],str(indices[i]),fontsize=8)
	plt.savefig("commensurate_nmax"+str(n_max)+".png")
	plt.clf()
	return [sol_list,indices,sol_array]



# Solve for line slip geometries
# For a given radius R, n0, n5, first solve for the orientation in which n0v0z+n5v5z=0
#	Choose a 'ribbon' by selecting a generator for the line-slip defect. 
#	Call the line slip generator va, and the standard dual vector with opposite z component vb
#		(This also selects an 'index' nb for the crystal, which is used to express
#		the commensurability/closure constraint.)
# 	Next, allow the ribbon to tilt until the 'top' and 'bottom' edges meet.
#		Replace the n0v0z+n5v5z=0 constraint with na*
#	Call the displacement vector between neighbors in the top and bottom edges vls
# Once a set of starting points has been found, integrate them along the associated tangents
#	(in va,vb,vls space)
# If va=v0, (Case 1)
#	(i) vb=v1=v0-v5=-vls. Then na*va+nb*vb-v0+v5 = n0*v0+n5*v5 => na+nb-1=n0, nb-1=-n5, nb=-n5+1 if n5>0,
#		(n5=0, va = v0, is an 'attractor': the radius does not change.) 
#		The next commensurate state should then be [n0+1, n5]:  
#		(vls: -v0+v5-> v3=v5, na*va+nb*vb+vls-> (n0+n5)*v0+(-n5+1)*(v0-v5)+v5=(n0+1)*v0+n5*v5)
# 	(ii) vb=v5. Then na*v0+(nb+1)*v5=n0*v0+n5*v5=>na=n0, nb=n5-1. 
#		Since vls'=v4=v5-v0, we have [n0,n5]-> [n0,n5-1]+[-1,1]=[n0-1,n5]
# if va=v5, (Case 2)
# 	(i) vb=v0=vls. Then nb'=n0-1 if n0>0, and the radius decreases
#		The next index should be [n0,n5-1] (vls=v0 goes to vls=v1=v0-v5)
#	(ii) vb=v1=vls. Then nb' = n0-1, and the radius increases
#		Next index: [n0,n5+1]
# if va=v1, (Case 3)
#	(i) vb=v0=vls. na*(v0-v5)+nb*v0+v0=n0*v0+n5*v5, na=-n5, nb'=n0+n5, nb=n0+n5-1
#	( vls' = v5, and -n5*(v0-v5)+(n0+n5-1)*v0+v5 = (n0-1)*v0+(n5+1)*v5 , [n0,n5]-> [n0-1,n5+1])
#	(ii) vb=v5=vls, na*(v0-v5)+nb*v5+v5=n0*v0+n5*v5, na=n0, nb'=n5+n0, nb=n0+n5-1, 
#	( vls' = v0, and n0*(v0-v5)+(n0+n5-1)*v5+v0 = (n0+1)*v0+(n5-1)*v5, [n0,n5]-> [n0+1,n5-1])
def ls_cconstraint(v2d,na,nb,R,case):
	C = 2*np.pi*R
	invR = 1./R
	rval = np.zeros(6)
	try:
		rval[0] = na*v2d[0]+nb*v2d[2]+v2d[4]-C 
		rval[1] = na*v2d[1]+nb*v2d[3]+v2d[5]
		rval[2] = (R*np.cos(v2d[0]*invR)-R)**2+(R*np.sin(v2d[0]*invR))**2+v2d[1]**2-4
		rval[3] = (R*np.cos(v2d[2]*invR)-R)**2+(R*np.sin(v2d[2]*invR))**2+v2d[3]**2-4
		if (case[0]==1 and case[1]==1) or (case[0]==2 and case[1]==1):
			# Revise this case! The distance is only approximately 2*sqrt(3). 
			# The norm of va+vb should also be 4:
			rval[4] = (R*np.cos((v2d[0]+v2d[2])*invR)-R)**2+(R*np.sin((v2d[0]+v2d[2])*invR))**2+(v2d[3]+v2d[1])**2-4
		else:
			rval[4] = (R*np.cos((v2d[0]-v2d[2])*invR)-R)**2+(R*np.sin((v2d[0]-v2d[2])*invR))**2+(v2d[3]-v2d[1])**2-4
		rval[5] = (R*np.cos(v2d[4]*invR)-R)**2+(R*np.sin(v2d[4]*invR))**2+v2d[5]**2-4
		return rval
	except:
		# check types:
		for i in range(6):
			if str(type(v2d[i]))!="<type 'numpy.float64'>":
				print "Error: v2d has wrong type "+str(type(v2d[i]))
				quit()
		if type(na)!=int or type(nb)!=int:
			print "Error: na, nb have types "+str(type(na))+" "+str(type(nb))
			quit()
		if type(R)!=float:
			print "Error: R has type "+str(type(R))
			quit()

def ls_cconstraint_all_cases(v2d,na,nb,R):
	cases = [[0,0],[0,1],[1,0],[1,1],[2,0],[2,1]]
	rval = 1
	for j in range(len(cases)):
		rval*=ls_cconstraint(v2d,na,nb,R,cases[j])
	return rval

def get_ls_curve(v2D,na,nb,R,Rnext,Nsamples,case):
	Rsamples = np.linspace(R,Rnext,Nsamples)
	dR = Rsamples[1]-Rsamples[0]
	v2Dvals = np.zeros([Nsamples,6])
	v2Dvals[0,:]=v2D
	success = True
	endpoint = Nsamples
	# Check that v2D satisfies ls_cconstraint at the initial radius:
	test = ls_cconstraint(v2D,na,nb,R,case)
	if np.dot(test,test)>1e-2:
		print "Error: incorrect starting loop (starting ls_cconstraint value: "+str(test)+" case: "+str(case)
		quit()
	for ii in range(1,len(Rsamples)):
		# Solve for va, vb, and vls at radius Rsamples[ii] 
		# (with previous sol'n as initial guess)
		init_guess = np.copy(v2Dvals[ii-1,:])
		v2D_sr = scipy.optimize.root(ls_cconstraint,init_guess,(na,nb,Rsamples[ii],case))
		if v2D_sr.success:
			v2Dvals[ii,:]=v2D_sr.x
			# Check for continuity:
			disp = v2Dvals[ii,:]-v2Dvals[ii-1,:]
			dispsq = np.dot(disp,disp)
			if dispsq>10*abs(dR):
				print "Likely discontinuous jump at R="+str(Rsamples[ii])+". Resolution="+str(dR)+", jump="+str(dispsq)
				success = False
				endpoint = ii
				break
		else:
			print "Unable to find root at R="+str(Rsamples[ii])+" with resolution dR="+str(dR)
			success = False
			endpoint = ii
			break
	# And interpolate between the endpoints.
	Rcsamples = np.copy(Rsamples[:endpoint])
	v2Dcvals = np.copy(v2Dvals[:endpoint,:])
	if Rnext<R:
		Rcsamples = Rcsamples[::-1]
		v2Dcvals = v2Dcvals[::-1,:]
	if endpoint>3:
		curve = scipy.interpolate.interp1d(Rcsamples,v2Dcvals,kind='cubic',axis=0)
	else:
		curve = zero_function
	return [curve,Rcsamples[endpoint-1],success]

# Test get_ls_curve:
def test_get_ls_curve(nmax_cc=10):
	[cc_list,cc_indices,cc_array] = enumerate_cc(nmax_cc)
	start_index = np.random.randint(len(cc_list))
	n0 = cc_indices[start_index][0]
	n5 = cc_indices[start_index][1]
	vR0 = cc_list[start_index]
	R = vR0[6]
	case = [0,0]
	v2D = np.zeros(6)
	projxy = np.sqrt(np.dot(vR0[:2],vR0[:2]))
	v2D[0] = 2*R*np.arcsin(projxy/(2*R))
	v2D[1] = vR0[2]
	projxy = np.sqrt(np.dot(vR0[3:5],vR0[3:5]))
	v2D[2] = 2*R*np.arcsin(projxy/(2*R))
	v2D[3] = vR0[5]
	v2D[2:4] = -np.copy(v2D[2:4])+np.copy(v2D[:2])
	v2D[4:] = -np.copy(v2D[2:4])
	print v2D
	na = n0+n5
	nb = -n5+1
	[n0p,n5p] = [n0+1,n5]
	print [n0,n5,na,nb]
	if n0p<10:
		if cc_array[n0p,n5p,6]>0:
			Nsamples = 100
			R = vR0[6]
			Rnext = cc_array[n0p,n5p,6]
			success = False
			n_attempts = 0
			print "Attempting to define ls curve"
			while not success and n_attempts<6:
				n_attempts+=1
				[curve,terminal_R,success] = get_ls_curve(v2D,na,nb,R,Rnext,Nsamples,case)
				Nsamples*=2
			# Plot the curve:
			if success:
				print "Plotting curve"
				Rmin = min(R,terminal_R)
				Rmax = max(R,terminal_R)
				Rsamples = np.linspace(Rmin,Rmax,Nsamples)
				thetas = np.zeros(Nsamples)
				for i in range(Nsamples):
					v2D = curve(Rsamples[i])
					thetas[i] = np.arcsin(v2D[1]*0.5)
				plt.plot(Rsamples,thetas)
				plt.show()
			else:
				print "Failed to converge"
				Rmin = min(R,Rnext)
				Rmax = max(R,Rnext)
				Rsamples = np.linspace(Rmin,Rmax,Nsamples)
				thetas = np.zeros(Nsamples)
				for i in range(Nsamples):
					v2D = curve(Rsamples[i])
					thetas[i] = np.arcsin(v2D[1]*0.5)
				plt.plot(Rsamples,thetas)
				plt.show()
	quit()

#test_get_ls_curve()

def get_v0_of_vabls(v2D,case123,caseab):
	v0 = np.zeros(2)
	if case123==0:
		v0 = np.copy(v2D[:2])
	elif case123==1:
		if caseab==0:
			v0 = v2D[2:4]
			# Check that v0 is in the correct quadrant:
		elif caseab==1:
			v0 = v2D[2:4]+v2D[:2]
		else:
			print "Error: caseab must be 0 or 1"
	elif case123==2:
		if caseab==0:
			v0 = v2D[2:4]
		elif caseab==1:
			v0 = v2D[:2]+v2D[2:4]
		else:
			print "Error: caseab must be 0 or 1"
	else:
		print "Error: case123 must be 0, 1, or 2"
	if (v0[0]<-epsilon) or (v0[1]<-epsilon):
		print "Error: incorrect quadrant "+str([v0,case123,caseab])
	return v0


# Obtain the initial two dimensional lattice vectors of the crystal with line slip
#	v2d or v2D:
#	va: [0:2] line slip defect generator (i.e. unbroken translation symmetry)
#	vb: [2:4] transverse vector used in commensurability condition
#	vls: [4:] displacement vector between adjcent neighbors across the line slip defect
#		(vls is chosen to form an accute angle with va initially.)
#	v2dR: (when parametrizing by theta instead of R)
#	va: [0:2]
#	vls: [2:4] 
#	R: [4]
def get_v2D_of_vR(vR):
	tR = 2*vR[6]
	itR = 1./tR
	proj0xy = np.sqrt(np.dot(vR[:2],vR[:2]))
	proj5xy = np.sqrt(np.dot(vR[3:5],vR[3:5]))
	dtheta0 = np.arcsin(proj0xy*itR)
	dS0 = tR*dtheta0
	dtheta5 = np.arcsin(proj5xy*itR)
	dS5 = tR*dtheta5
	rval = np.zeros([3,2,6])
	v2D = np.zeros(6)
	# case[0]==0 and case[1]==0:
	# va = v0, vb = v1
	v2D[0] = dS0
	v2D[1] = np.copy(vR[2])
	projxy = np.sqrt(np.dot(vR[3:5],vR[3:5]))
	v2D[2] = dS0-dS5
	v2D[3] = vR[2]-vR[5]
	v2D[4:] = -v2D[2:4]
	rval[0,0,:] = np.copy(v2D)
	# case[0]==0 and case[1]==1:
	# Case 1b: va = v0, vb = v5, [n0p,n5p]=[n0-1,n5], na=n0, nb=n5-1
	v2D[2] = dS5
	v2D[3] = np.copy(vR[5])
	v2D[4:] = np.copy(v2D[2:4])
	rval[0,1,:]=np.copy(v2D)
	# case[0]==1 and case[1]==0:
	# Case 2a: va = v5, vb = v0 
	v2D[0]=dS5
	v2D[1]= np.copy(vR[5])
	v2D[2] = dS0
	v2D[3] = np.copy(vR[2])
	v2D[4:] = np.copy(v2D[2:4])
	rval[1,0,:]=np.copy(v2D)
	# case[0]==1 and case[1]==1:
	# Case 2b: va = v5, vb = v1
	v2D[2] = dS0-dS5
	v2D[3] = vR[2]-vR[5]
	v2D[4:] = np.copy(v2D[2:4])
	rval[1,1,:] = np.copy(v2D)
	# Case 3a: va = v1, vb = v0 = vls
	v2D[0] = dS0-dS5
	v2D[1] = vR[2]-vR[5]
	v2D[2] = dS0
	v2D[3] = vR[2]
	v2D[4:] = np.copy(v2D[2:4])
	rval[2,0,:] = np.copy(v2D)
	# Case 3b: va = v1, vb = v5 = vls
	v2D[2] = dS5
	v2D[3] = vR[5]
	v2D[4:] = np.copy(v2D[2:4])
	rval[2,1,:] = np.copy(v2D)
	return rval

def get_v2DR_of_vR(vR):
	tR = 2*vR[6]
	itR = 1./tR
	proj0xy = np.sqrt(np.dot(vR[:2],vR[:2]))
	proj5xy = np.sqrt(np.dot(vR[3:5],vR[3:5]))
	dtheta0 = np.arcsin(proj0xy*itR)
	dS0 = tR*dtheta0
	dtheta5 = np.arcsin(proj5xy*itR)
	dS5 = tR*dtheta5
	v2DRs = np.zeros([2,5])
	# Case 3a: va = v1, vb = v0 = vls
	v2DRs[0,0] = dS0-dS5
	v2DRs[0,1] = vR[2]-vR[5]
	v2DRs[0,2] = dS0
	v2DRs[0,3] = vR[2]
	v2DRs[0,4] = vR[6]
	# Case 3b: va = v1, vb = v5 = vls
	v2DRs[1,:2] = np.copy(v2DRs[0,:2])
	v2DRs[1,2] = dS5
	v2DRs[1,3] = vR[5]
	v2DRs[1,4] = vR[6]
	return v2DRs

def get_nbnas_of_n0n5(n0,n5):
	rval = np.zeros([3,2,2],dtype=int)
	rval[0,0,:] = np.array([n0+n5,-n5+1])
	rval[0,1,:] = np.array([n0,n5-1])
	rval[1,0,:] = np.array([n5,n0-1])
	rval[1,1,:] = np.array([n0+n5,n0-1])
	rval[2,0,:] = np.array([-n5,n0+n5-1])
	rval[2,1,:] = np.array([n0,n0+n5-1])
	return rval

def get_next_n0n5(n0,n5):
	rval = np.zeros([3,2,2],dtype=int)
	if n0>0 and n5>0:
		rval[0,0,:]=np.array([n0+1,n5])
		rval[0,1,:]=np.array([n0-1,n5])
		rval[1,0,:]=np.array([n0,n5-1])
		rval[1,1,:]=np.array([n0,n5+1])
		rval[2,0,:]=np.array([n0-1,n5+1])
		rval[2,1,:]=np.array([n0+1,n5-1])
	elif n0>0:
		rval[0,0,:]=np.array([n0,n5])
		rval[0,1,:]=np.array([n0,n5])
		rval[1,0,:]=np.array([n0,n5])
		rval[1,1,:]=np.array([n0,n5+1])
		rval[2,0,:]=np.array([n0-1,n5+1])
		rval[2,1,:]=np.array([n0,n5])
	else:
		rval[0,0,:]=np.array([n0+1,n5])
		rval[0,1,:]=np.array([n0-1,n5])
		rval[1,0,:]=np.array([n0,n5])
		rval[1,1,:]=np.array([n0,n5])
		rval[2,0,:]=np.array([n0,n5])
		rval[2,1,:]=np.array([n0+1,n5-1])
	return rval

def ls_cconstraint_theta(v2dR,na,nb,theta,caseab):
	R = v2dR[4]
	# sweeping over v0(theta): assume case 2
	v02dx = 2*R*np.arcsin(np.cos(theta)/v2dR[4])
	v02dz = 2*np.sin(theta)
	v0 = np.zeros(3)
	v0[0]=2*np.cos(theta)
	# 2*R*vx+vx^2+vy^2=0
	v0[2] = 2*np.sin(theta)
	projxysq = 4-v0[2]*v0[2]
	v0[0]=-projxysq/(2*R)
	v0[1]=np.sqrt(projxysq-v0[0]*v0[0])
	# va=v1, vb = v0 or v5
	rval = np.zeros(5)
	# length constraint of va:=v1:
	va = np.zeros(3)
	va[2] = v2dR[1]
	va[0] = R*(np.cos(v2dR[0]/R)-1)
	va[1] = R*np.sin(v2dR[0]/R)
	vls = np.zeros(3)
	vls[2] = v2dR[3]
	vls[0] = R*(np.cos(v2dR[2]/R)-1)
	vls[1] = R*np.sin(v2dR[2]/R)
	if caseab==0:
		vb = np.copy(v0)
		vb2d = np.array([v02dx,v02dz])
	elif caseab==1:
		vb = np.zeros(3)
		vb2d = -v2dR[:2]
		vb2d[0]+=v02dx
		vb2d[1]+=v02dz
		vb[2] = vb2d[1]
		vb[0] = R*(np.cos(vb2d[0]/R)-1)
		vb[1] = R*np.sin(vb2d[0]/R)
	else:
		print "Error: caseab must be 0 or 1"
		quit()
	rval[0] = np.dot(va,va)-4
	# length constraint of vls
	rval[1] = np.dot(vls,vls)-4
	# length constraint of va from v0:
	rval[2] = np.dot(va-v0,va-v0)-4
	rval[3] = na*v2dR[0]+nb*vb2d[0]+v2dR[2]-2*np.pi*v2dR[4]
	rval[4] = na*v2dR[1]+nb*vb2d[1]+v2dR[3]
	return rval


# Modify this to incorporate 'get_R2valmax' and ensure that the 
# range of the approximant includes that of the exact curve
def get_ls_curve_theta(v2dR,na,nb,theta,theta_next,Nsamples,case):
	# Check that the case is appropriate:
	if case[0]!=2:
		print "Error: get_ls_curve_theta should only be used in case 3 (va=v1,vb=v0,v5)"
		quit()
	theta_samples = np.linspace(theta,theta_next,Nsamples)
	dtheta = theta_samples[1]-theta_samples[0]
	rtdtheta = np.sqrt(dtheta)
	v2dRvals = np.zeros([Nsamples,5])
	v2dRvals[0,:] = np.copy(v2dR)
	theta_term = Nsamples
	success = True
	test = ls_cconstraint_theta(v2dR,na,nb,theta,case[1])
	if np.dot(test,test)>1e-3:
		print "Error in get_ls_curve_theta: check preparation of v2dR (it does not start at a root. val = "+str(test)+")"
		quit()
	imax = 0
	Rmax = v2dR[4]
	v2dRmax = np.copy(v2dR)
	for i in range(1,Nsamples):
		# Use each root as an initial guess for the next
		init_guess = np.copy(v2dRvals[i-1,:])
		v2dR_sr = scipy.optimize.root(ls_cconstraint_theta,init_guess,(na,nb,theta_samples[i],case[1]))
		if v2dR_sr.success:
			v2dRvals[i,:]=np.copy(v2dR_sr.x)
			# Check if imax and Rmax need to be updated:
			if v2dRvals[i,4]>v2dRmax[4]:
				v2dRmax = np.copy(v2dRvals[i,:])
				imax = i
			# check for continuity
			diff = v2dRvals[i,:]-v2dRvals[i-1,:]
			diffsq = np.dot(diff,diff)
			if diffsq>10*abs(dtheta):
				print "Likely discontinuity at theta="+str(theta_samples[i])+" of "+str(np.sqrt(diffsq))
				success = False
				theta_term = i
				break
			elif diffsq==0:
				print "Likely cusp (or worse) in curve [na,nb]="+str([na,nb])
		else:
			print "Unable to find root at theta="+str(theta_samples[i])
			success = False
			theta_term = i
			break
		
	# Define a curve to approximate v2dR(theta)
	if theta_term>3:
		thetaCsamples = np.copy(theta_samples[:theta_term])
		v2dRCvals = np.copy(v2dRvals[:theta_term,:])
		if theta_next<theta:
			thetaCsamples = thetaCsamples[::-1]
			v2dRCvals = v2dRCvals[::-1,:]
		curve = scipy.interpolate.interp1d(thetaCsamples,v2dRCvals,kind='cubic',axis=0)
		# Check if v2dRvals includes a local maximum:
		if Rmax>v2dRvals[0,4] and Rmax>v2dRvals[theta_term-1,4]:
			# include a precise estimate of the maximum near Rmax:
			thetaRm = theta_samples[imax]
			theta1 = theta_samples[max(imax-10,0)]
			theta2 = theta_samples[min(imax+10,theta_term-1)]
			theta_min = min(theta1,theta2)
			theta_max = max(theta1,theta2)
			[v2DRmax,tvm,success] = get_R2valmax(v2DRm_init,thetaRm,na,nb,case[1],theta_min,theta_max,curve)
			# If if tvm is less than theta_max, create a new interpolant that includes v2DRmax and tvm:
			if theta_min<tvm<theta_max:
				print "Attempting to update branched curve to extend over full range. (Case = "+str(case)+")"
				v2dRvals_new = np.zeros([theta_term+1,5])
				tsamples = np.zeros([theta_term+1,5])
				i = 0
				while theta_samples[i]<tvm:
					v2dRvals_new[i,:] = np.copy(v2dRvals[i,:])
					tsamples[i]=theta_samples[i]
					i+=1
				tsamples[i]=tvm
				v2dRvals_new[i,:] = np.copy(v2DRmax)
				while i<theta_term+1:
					v2dRvals_new[i,:] = np.copy(v2dRvals[i-1,:])
				# update 'curve':
				thetaCsamples = tsamples[:theta_term+1]
				v2dRCvals = v2dRvals_new[:theta_term+1]
				if theta_next<theta:
					thetaCsamples = thetaCsamples[::-1]
					v2dRCvals = v2dRCvals[::-1,:]
				curve = scipy.interpolate.interp1d(thetaCsamples,v2dRCvals,kind='cubic',axis=0)
				theta_term+=1
	else:
		curve = zero_function
		thetaCsamples = np.copy(theta_samples[:theta_term])
		success = False
	return [curve,thetaCsamples[theta_term-1],success]

def get_R2valmax(v2DRm_init,thetaRm,na,nb,icase,theta_min,theta_max,ls_curve,tol=1e-9):
	# Start by sampling v2DR near v2DRm_init:
	v2DRm_init_r = scipy.optimize.root(ls_cconstraint_theta,v2DRm_init,(na,nb,thetaRm,icase))
	if v2DRm_init_r.success:
		v2DRm_init = v2DRm_init_r.x
	else:
		print "Error: bad initial guess for [get_R2valmax]"
		quit()
	# Converge on maximum
	err = 1
	R0 = v2DRm_init[4]
	R1 = R0
	t1 = theta_min
	t2 = theta_max
	Nsamples = 7
	Ntries = 0
	local_curve = ls_curve
	success = False
	v2DRmax = np.copy(v2DRm_init)
	tvm = thetaRm
	while Ntries<20:
		Ntries+=1
		samples = np.linspace(t1,t2,Nsamples)
		v2DRs = np.zeros([Nsamples,5])
		iRmax = 0
		for i in range(Nsamples):
			# Solve for v2DR:
			v2DR_r = scipy.optimize.root(ls_cconstraint_theta,local_curve(samples[i]),(na,nb,samples[i],icase))
			if v2DR_r.success:
				v2DRs[i,:] = np.copy(v2DR_r.x)
			else:
				print "Error: refine line-slip interpolation"
				quit()
			if v2DRs[i,4]>R1:
				R1 = v2DRs[i,4]
				iRmax = i
				v2DRmax = np.copy(v2DRs[i,:])
				tvm = samples[i]
		err = abs(R1-R0)
		if err<tol:
			success = True
			break
		else:
			# Update the local curve, and samples:
			local_curve = scipy.interpolate.interp1d(samples,v2DRs,axis=0,kind='cubic')
			t1 = samples[max(iRmax-1,0)]
			t2 = samples[min(iRmax+1,Nsamples-1)]
			#samples = np.linspace(t1,t2,Nsamples)
			R0 = R1
	return [v2DRmax,tvm,success]

def enumerate_ccls(nmax,Nsamples,Rmax=None,dpi=1024):
	curve_list = []
	# Solve for commensurate crystals
	[cc_list,indices,cc_array] = enumerate_cc(nmax,Rmax)
	# If va=v0, (Case 1 [icase=0])
	#	(i) vb=v1=v0-v5=-vls. Then na*va+nb*vb-v0+v5 = n0*v0+n5*v5 => na+nb-1=n0, nb-1=-n5, nb=-n5+1 if n5>0,
	#		(n5=0, va = v0, is an 'attractor': the radius does not change.) 
	#		The next commensurate state should then be [n0+1, n5]:  
	#		(vls: -v0+v5-> v3=v5, na*va+nb*vb+vls-> (n0+n5)*v0+(-n5+1)*(v0-v5)+v5=(n0+1)*v0+n5*v5)
	# 	(ii) vb=v5. Then na*v0+(nb+1)*v5=n0*v0+n5*v5=>na=n0, nb=n5-1. 
	#		Since vls'=v4=v5-v0, we have [n0,n5]-> [n0,n5-1]+[-1,1]=[n0-1,n5]
	# if va=v5, (Case 2 [icase=1])
	# 	(i) vb=v0=vls. Then nb'=n0-1 if n0>0, and the radius decreases
	#		The next index should be [n0,n5-1] (vls=v0 goes to vls=v1=v0-v5)
	#	(ii) vb=v1=vls. Then nb' = n0-1, and the radius increases
	#		Next index: [n0,n5+1]
	# if va=v1, (Case 3 [icase=2])
	#	(i) vb=v0=vls. na*(v0-v5)+nb*v0+v0=n0*v0+n5*v5, na=-n5, nb'=n0+n5, nb=n0+n5-1
	#	( vls' = v5, and -n5*(v0-v5)+(n0+n5-1)*v0+v5 = (n0-1)*v0+(n5+1)*v5 , [n0,n5]-> [n0-1,n5+1])
	#	(ii) vb=v5=vls, na*(v0-v5)+nb*v5+v5=n0*v0+n5*v5, na=n0, nb'=n5+n0, nb=n0+n5-1, 
	#	( vls' = v0, and n0*(v0-v5)+(n0+n5-1)*v5+v0 = (n0+1)*v0+(n5-1)*v5, [n0,n5]-> [n0+1,n5-1])
	# Keep track of which curves have already been defined
	defined = {}
	for i in range(len(cc_list)):
		vR = cc_list[i] 
		theta_init = np.arcsin(0.5*vR[2])
		[n0,n5] = indices[i]
		R = vR[6]
		plt.scatter([R],[theta_init])
		plt.text(R,theta_init,str([n0,n5]),fontsize=3)
		v2Ds = get_v2D_of_vR(vR)
		v2DRs = get_v2DR_of_vR(vR)
		nbnas = get_nbnas_of_n0n5(n0,n5)
		n0pn5ps = get_next_n0n5(n0,n5)
		for ocase in range(2):
			for icase in range(2):
				#if (ocase==0 and n5==0) or (ocase==1 and n0==0): continue
				[na,nb]=nbnas[ocase,icase]
				[n0p,n5p] = n0pn5ps[ocase,icase]
				if ([n0p,n5p] in indices) and ([n0p,n5p]!=[n0,n5]) and not ((n0p,n5p,n0,n5) in defined):
					if cc_array[n0p,n5p,6]>0:
						v2D = np.copy(v2Ds[ocase,icase,:])
						R = np.copy(vR[6])
						Rnext = cc_array[n0p,n5p,6]
						if Rnext==R:
							print "Cases "+str([n0,n5])+" and "+str([n0p,n5p])+" have equival radii"
						# Solve for va,vb, and vls for a sequence of points between R and Rnext,
						Ncsamples = Nsamples
						success = False
						Nattempts = 0
						while not success and Nattempts<5:
							[curve,Rterm,success] = get_ls_curve(v2D,na,nb,R,Rnext,Ncsamples,[ocase,icase])
							Ncsamples*=2
							Nattempts+=1
						# Store these approximate curves for future reference (they can provide initial guesses.)
						if success:
							defined[(n0,n5,n0p,n5p)] = len(curve_list)
							curve_list.append([na,nb,ocase,icase,curve,R,Rnext,0,[(n0,n5,n0p,n5p)]])
							print "Success at branch "+str([n0,n5,n0p,n5p,ocase,icase])
				elif (n0p,n5p,n0,n5) in defined:
					# Recognize that the curve that starts at n0p,n5p is also relevant to n0, n5
					# find n0p and n5p in curve_list:
					curve_index = defined[(n0p,n5p,n0,n5)]
					# Include the connection between [n0,n5] and [n0p,n5p], with the associated case
					curve_list[curve_index][-1].append([i,(na,nb,ocase,icase)]) 
		# Case 3: parametrize with theta instead of R
		for icase in range(2):
			#if ((icase==0 and n0==0) or (icase==1 and n5==0)): continue
			[na,nb]=nbnas[2,icase]
			[n0p,n5p] = n0pn5ps[2,icase]
			if ([n0p,n5p] in indices) and ([n0p,n5p]!=[n0,n5]) and not ((n0p,n5p,n0,n5) in defined):
					v2DR = np.copy(v2DRs[icase,:])
					theta = np.arcsin(0.5*vR[2])
					theta_next = np.arcsin(0.5*cc_array[n0p,n5p,2])
					if theta_next==theta:
						print "Cases "+str([n0,n5])+" and "+str([n0p,n5p])+" have equal orientations"
					Ncsamples = Nsamples
					success = False
					Nattempts = 0
					while not success and Nattempts<5:
						[curve,theta_term,success] = get_ls_curve_theta(v2DR,na,nb,theta,theta_next,Ncsamples,[2,icase])
						Ncsamples*=2
						Nattempts+=1
					if success:
						defined[(n0,n5,n0p,n5p)] = len(curve_list)
						curve_list.append([na,nb,2,icase,curve,theta,theta_next,1,[(n0,n5,n0p,n5p)]])
						print "Success at branch "+str([n0,n5,2,icase,theta,theta_next])
			elif (n0p,n5p,n0,n5) in defined:
				curve_index = defined[(n0p,n5p,n0,n5)]
				curve_list[curve_index][-1].append([i,(na,nb,2,icase)])
			
			# Case 1a: va = v0, vb = v1
			# Case 1b: va = v0, vb = v5, [n0p,n5p]=[n0-1,n5], na=n0, nb=n5-1
			# Case 2a: va = v5, vb = v0
			# Case 2b: vb = v1
			# Case 3a: va = v1, vb = v0 = vls
			# Case 3b: va = v1, vb = v5 = vls
	# Plot each curve segment:
	print "Plotting segments ("+str(len(curve_list))+" total curves)"
	Rsamples_full = []
	thetas_full = []
	v2Ds_full = []
	features_full = [] # integer features of each curve
	features_cont = [] # continuous features of curves (mainly relevant to case 3 curves)
	# elements of features_cont have the form [Rmin,Rmax,R2valmin,R2valmax]
	for i in range(len(curve_list)):
		# Get n0, n5 associated with curve_list
		nnp_tuple = curve_list[i][8][0]
		n0=nnp_tuple[0]
		n5=nnp_tuple[1]
		n0p=nnp_tuple[2]
		n5p=nnp_tuple[3]
		#[n0,n5,n0p,n5p] = curve_list[i][8]
		# NOTE: better to store the curve index than the curve indices (easier for subsequent calculations.)
		if len(curve_list[i])>9:
			[nap,nbp,ocasep,icasep] = curve_list[i][8][2]
			curve_index_p=curve_list[i][8][1]
		else:
			[nap,nbp,ocasep,icasep] = [-1,-1,-1,-1]
			curve_index_p=i # i.e. the mirror image of the curve is itself. Check remaining parts for consistency
		if curve_list[i][7]==0:
			[na,nb,ocase,icase] = curve_list[i][:4]
			#print str(i)+": "+str(curve_list[i][:4])+str(curve_list[i][5:])
			R1 = curve_list[i][5]
			R2 = curve_list[i][6]
			Rmin = min(R1,R2)
			Rmax = max(R1,R2)
			if Rmax==Rmin:
				print "Error: curve at "+str([n0,n5,ocase,icase])+" is actually a point"
				continue
			Rsamples = np.linspace(Rmin+epsilon*(Rmax-Rmin),Rmax-epsilon*(Rmax-Rmin),2*Nsamples)
			thetas = np.zeros(len(Rsamples))
			imin = len(Rsamples_full)
			for ii in range(len(Rsamples)):
				v2D = curve_list[i][4](Rsamples[ii]) # NOTE: va is not always v0
				v02D = get_v0_of_vabls(v2D,curve_list[i][2],curve_list[i][3])
				# Compute the orientation associated with v0:
				thetas[ii] = np.arcsin(0.5*v02D[1])
				Rsamples_full.append([Rsamples[ii],n0,n5,ocase,icase])
				thetas_full.append(thetas[ii])
				v2Ds_full.append(v2D)
			imax = len(Rsamples_full)
			features_full.append([na,nb,ocase,icase,imin,imax,0,n0,n5,n0p,n5p,curve_index_p,nap,nbp,ocasep,icasep])
			features_cont.append([thetas[0],thetas[1],R1,R2,0,0,-10])
			plt.plot(Rsamples,thetas)
			if (curve_list[i][2]==0 and curve_list[i][3]==0) or curve_list[i][2:4]==[1,1]:
				plt.text(Rsamples[Nsamples],thetas[Nsamples],str(curve_list[i][2:4]),fontsize=3)
		elif curve_list[i][7]==1:
			[na,nb,ocase,icase] = curve_list[i][:4]
			print str(i)+": "+str(curve_list[i][:4])+str(curve_list[i][5:])
			theta1 = curve_list[i][5]
			theta2 = curve_list[i][6]
			if theta1==theta2:
				print "Error: curve at "+str([n0,n5,ocase,icase])+" is actually a point"
			tmin = min(theta1,theta2)
			tmax = max(theta1,theta2)
			tsamples0 = np.linspace(tmin+epsilon*(tmax-tmin),tmax-epsilon*(tmax-tmin),2*Nsamples)
			Rsamples0 = np.zeros(len(tsamples0))
			imin = len(Rsamples_full)
			v2DR0 = curve_list[i][4](tsamples0[0])
			v2DR1 = curve_list[i][4](tsamples0[len(tsamples0)-1])
			R0 = v2DR0[4]
			R1 = v2DR1[4]
			Rmax = max(R0,R1)
			iRmax = 0
			v2DRmax = v2DR0
			v2DRs0 = np.zeros([len(tsamples0),5])
			for ii in range(len(tsamples0)):
				v2DR = curve_list[i][4](tsamples0[ii])
				v2DRs0[ii,:] = np.copy(v2DR)
				Rsamples0[ii] = v2DR[4]
				if Rsamples0[ii]>Rmax:
					Rmax = Rsamples0[ii]
					iRmax = ii
					v2DRmax = np.copy(v2DR)
				Rsamples_full.append([Rsamples0[ii],n0,n5,ocase,icase])
				thetas_full.append(tsamples0[ii])
				# Get v2D from v2DR:
				v2D = np.zeros(6)
				v2D[:2] = np.copy(v2DR[:2])
				v2D[3] = 2*np.sin(tsamples0[ii])
				projxy = np.cos(tsamples0[ii])
				v2D[2] = 2*v2DR[4]*np.arcsin(projxy/v2DR[4])
				if curve_list[i][3]==1: # then vb = v5 = -v1+v0
					v2D[2:4] = -v2D[:2]+v2D[2:4]
				v2D[4:] = np.copy(v2DR[2:4])
				v2Ds_full.append(v2D)
			R2valmin = max(R0,R1)
			theta_min = tsamples0[max(iRmax-20,0)]
			theta_max = tsamples0[min(iRmax+20,len(tsamples0)-1)]
			thetaRm = tsamples0[iRmax]
			# (v2DR: initial guess, na, nb, icase, theta_min, theta_max: (bounds for maximum))
			[v2DR2valm,theta_vm,R2valm_success] = get_R2valmax(v2DRmax,thetaRm,na,nb,curve_list[i][3],theta_min,theta_max,curve_list[i][4]) 
			if theta_vm==0:
				print "Likely error on curve (n0,n5,n0p,n5p)="+str([n0,n5,n0p,n5p])
			R2valmax = v2DR2valm[4]
			
			imax = len(Rsamples_full)
			plt.plot(Rsamples0,tsamples0)
			if R2valmax>max(R1,R0)+epsilon:
				plt.scatter([R2valmax],[thetaRm])
				plt.text(R2valmax,thetaRm,"c"+str(i)+":"+str([n0,n5,ocase,icase]),fontsize=3)
				features_full.append([na,nb,ocase,icase,imin,imax,1,n0,n5,n0p,n5p,curve_index_p,nap,nbp,ocasep,icasep])
				features_cont.append([tsamples0[0],tsamples0[len(tsamples0)-1],R0,R1,R2valmin,R2valmax,theta_vm])
			elif curve_list[i][3]==0: 
				plt.text(Rsamples0[Nsamples/3],tsamples0[Nsamples/3],"(r)c"+str(i)+":"+str([ocase,icase]),fontsize=3)
				features_full.append([na,nb,ocase,icase,imin,imax,0,n0,n5,n0p,n5p,curve_index_p,nap,nbp,ocasep,icasep])
				features_cont.append([theta1,theta2,R0,R1,0,0,-10])
			else:
				features_full.append([na,nb,ocase,icase,imin,imax,0,n0,n5,n0p,n5p,curve_index_p,nap,nbp,ocasep,icasep])
				features_cont.append([theta1,theta2,R0,R1,0,0,-10])
	Rsamples_full = np.array(Rsamples_full)
	thetas_full = np.array(thetas_full)
	v2Ds_full = np.array(v2Ds_full)
	features_full = np.array(features_full,dtype=int)
	features_cont = np.array(features_cont)
	# Save curve data for future reference:
	np.savetxt("LSDATA/ls_curvef.dat",features_full,fmt="%d")
	np.savetxt("LSDATA/ls_curvefc.dat",features_cont)
	np.savetxt("LSDATA/ls_Rsamples.dat",Rsamples_full)
	np.savetxt("LSDATA/ls_v2Ds.dat",v2Ds_full)
	np.savetxt("LSDATA/ls_thetas.dat",thetas_full)
	plt.xlabel("Cylinder radius ($r_{sphere}=1$)")
	plt.ylabel("Crystal orientation")
	plt.savefig("ls_curves_nmax"+str(nmax)+".png",dpi=dpi)
	plt.clf()

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




enumerate_ccls(5,50)
