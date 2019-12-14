import os
import numpy as np
import pickle
# Miscellaneous operations
def user_input_continue_yn(message):
	ui = raw_input(message+"\n")
	if ui=='y':
		return True
	elif ui=='n':
		return False
	else:
		rval= user_input_continue_yn(message)
		return rval
			

def custom_dtconv(expr,dtflag):
	if dtflag=='i':
		return int(dtflag)
	elif dtflag=='f':
		return float(dtflag)
	elif dtflag=='b':
		return bool(dtflag)
	elif dtflag=='s' or dtflag=='u':
		return dtflag
	else:
		print "Warning: unrecognized dec type flag "+str(dtflag)+" encountered in custom_dtconv. (Recommend checking syntax in reading/writing node/edge decorations)"
		return dtflag

def transcribe_list(list_):
	new_list = []
	for j in range(len(list_)):
		new_list.append(list_[j])
	return new_list


def transcribe_dict(dict_):
	new_dict = {}
	for j in dict_:
		new_dict[j]=dict_[j]
	return new_dict
def mkdir_s(dirname):
	try:
		os.mkdir(dirname)
	except:
		pass

nth_exceptions = {}
nth_exceptions[1] = "1st"
nth_exceptions[2] = "2nd"
nth_exceptions[3] = "3rd"
def nth(j):
	if j not in nth_exceptions:
		return str(j)+"th"
	else:
		return nth_exceptions[j]

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

def add_dl(dl_instance,label,value):
	if label in dl_instance:
		dl_instance[label].append(value)
	else:
		dl_instance[label]=[value]

def incr_dict_tally(dt_instance,label,amt=1):
	if label in dt_instance:
		dt_instance[label]+=amt
	else:
		dt_instance[label]=amt

alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
# Load a list of colorful nouns and adjectives

cn_file= open("colorful_nouns.dat","r")
colorful_nouns = {}
for line in cn_file:
	line = line.split()
	if len(line)>1:
		add_dl(colorful_nouns,line[0],line[1])
cn_file.close()

ca_file=open("colorful_adjectives.dat","r")
colorful_adjec = {}
for line in ca_file:
	line = line.split()
	if len(line)>1:
		add_dl(colorful_adjec,line[0],line[1])
ca_file.close()

def random_name():
	adj_i = np.random.randint(len(colorful_adjec))
	n1_i = np.random.randint(len(colorful_nouns))
	n2_i = np.random.randint(len(colorful_nouns))
	while n2_i==n1_i:
		n2_i = np.random.randint(len(colorful_nouns))
	name = colorful_adjec[adj_i]+" "+colorful_nouns[n1_i]+" and "+colorful_nouns[n2_i]
	return name

def fill_map(domain,partial_map,fill_value=-1):
	for key in domain:
		if key in partial_map:
			pass
		else:
			partial_map=fill_value
	return partial_map

class unordered_list:
	def __init__(self,list_=None):
		if np.any(list_==None):
			self.e=[]
			self.state={}
		else:
			self.e = list_
			self.c = {}
			for j in list_:
				incr_dict_tally(self.c,list_[j])
	def equiv(self,ulist):
		return self.c==ulist.c
	def add(self,element):
		self.e.append(element)
		incr_dict_tally(self.c,element)

def type_s(x):
	tx = type(x)
	if tx!=str:
		return np.dtype(tx)
	else:
		return 'S'+str(len(x))

def dtype(A,types=None):
	Atype = []
	if np.any(types!=None):
		pass
	else:
		types={}
		for key in A:
			try:
				types[key]=type_s(A[key])
			except:
				print "Warning (dtype, common.py): Uncertain data type encountered in input "+str(A)+". (Check that input is actually a dictionary)"
	for key in A:
		Atype.append((str(key),types[key]))
	Atype=np.dtype(Atype)
	return Atype

def advance_multi_index(mi,mi_max):
	if mi[0]<mi_max[0]:
		mi[0]+=1
		complete = False
	elif len(mi)>1:
		mi[0]=0
		[mi[1:],complete]=advance_multi_index(mi[1:],mi_max[1:])
	else:
		mi[0]=0
		complete = True
	return [mi,complete]

# Establish a correspondence between two collections of subsets based on overlap

def overlap_corr(col1,col2):
	overlaps = {}
	for j in range(len(col1)):
		sc1j = set(col1[j])
		for jj in range(len(col2)):
			sc2jj= set(col2[jj])
			overlaps[(j,jj)]=len(sc1j.intersection(sc2jj))
	corr = {}
	invcorr = {}
	for j in range(len(col1)):
		max_overlap = 0
		corr[j]=-1
		for jj in range(len(col2)):
			if overlaps[(j,jj)]>max_overlap:
				corr[j]=jj
				max_overlap=overlaps[(j,jj)]
		if corr[j]>-1:
			add_dl(invcorr,corr[j],j)
	for jj in range(len(col2)):
		if jj in invcorr:
			pass
		else:
			# (i.e. new domain spontaneously formed)
			invcorr[jj]=-1 
	return [corr,invcorr]

def merge_id_filtered_dictionaries(dlist):
	idd = {}
	visited = {}
	next_id = 0
	for j in range(len(dlist)):
		for id_ in dlist[j]:
			if id_ not in visited:
				idd[id_]=dlist[j][id_]
			else:
				while visited[next_id]:
					next_id+=1
				idd[next_id]=dlist[j][id_]
	return idd

graded_array = np.dtype([('score','f8'),('array',type(np.array([])))])
graded_int = np.dtype([('score','f8'),('int','i4')])


