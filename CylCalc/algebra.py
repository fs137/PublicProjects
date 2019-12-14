import numpy as np
# Elementary rotations in 3D:
def qtimes(a,b):
	rval = np.array([a[0]*b[0]-np.dot(a[1:],b[1:]),a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1]])
	rval[1:]+=a[0]*b[1:]+b[0]*a[1:]
	return rval

def qconj(c):
	rval = np.copy(c)
	rval[1:]*=-1
	return rval

def qrotate(c,x):
	qx = np.zeros(4)
	qx[1:]=np.copy(x)
	cc = qconj(c)
	qrval = qtimes(c,qtimes(qx,cc))
	rval = qrval[1:]
	return rval

def factorial(N):
	if N>1:
		return N*factorial(N-1)
	else:
		return 1

def choose(N,M):
	if N>M:
		num = choose(N-1,M)
		return (N*num)/(N-M)
	elif N==M:
		return 1

