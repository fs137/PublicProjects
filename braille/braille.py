import numpy as np
from matplotlib import pyplot as plt

braille = {}
braille['a']=[[1,0],[0,0],[0,0]]
braille['b']=[[1,0],[1,0],[0,0]]
braille['c']=[[1,1],[0,0],[0,0]]
braille['d']=[[1,1],[0,1],[0,0]]
braille['e']=[[1,0],[0,1],[0,0]]
braille['f']=[[1,1],[1,0],[0,0]]
braille['g']=[[1,1],[1,1],[0,0]]
braille['h']=[[1,0],[1,1],[0,0]]
braille['i']=[[0,1],[1,0],[0,0]]
braille['j']=[[0,1],[1,1],[0,0]]
braille['k']=[[1,0],[0,0],[1,0]]
braille['l']=[[1,0],[1,0],[1,0]]
braille['m']=[[1,1],[0,0],[1,0]]
braille['n']=[[1,1],[0,1],[1,0]]
braille['o']=[[1,0],[0,1],[1,0]]
braille['p']=[[1,1],[1,0],[1,0]]
braille['q']=[[1,1],[1,1],[1,0]]
braille['r']=[[1,0],[1,1],[1,0]]
braille['s']=[[0,1],[1,0],[1,0]]
braille['t']=[[0,1],[1,1],[1,0]]
braille['u']=[[1,0],[0,0],[1,1]]
braille['v']=[[1,0],[1,0],[1,1]]
braille['x']=[[1,1],[0,0],[1,1]]
braille['y']=[[1,1],[0,1],[1,1]]
braille['z']=[[1,0],[0,1],[1,1]]
braille['w']=[[0,1],[1,1],[0,1]]

braille['cap']=[[0,0],[0,0],[0,1]]
braille['capterm']=[[0,0],[0,0],[1,0]]
braille["'"]=braille['capterm']

bks = braille.keys()

for key in braille:
	braille[key]=np.array(braille[key],dtype=int)

sc = {}
sc[0]=" "
sc[1]="o"

def display(braille_symbol,mode='text'):
	l1 = braille_symbol[0,:]
	l2 = braille_symbol[1,:]
	l3 = braille_symbol[2,:]
	if mode=='text':
		print sc[l1[0]]+sc[l1[1]]
		print sc[l2[0]]+sc[l2[1]]
		print sc[l3[0]]+sc[l3[1]]
	else:
		plt.matshow(braille_symbol)
	plt.show()

class braille_learner_vis:
	def __init__(self,N_letters):
		self.letters = {}
		for j in range(N_letters):
			i = np.random.randint(N_letters)
			key_i = bks[i]
			while key_i in self.letters:
				i = np.random.randint(len(bks))
				key_i = bks[i]
			self.letters[key_i]={}
			self.letters[key_i]['N_attempts']=0
			self.letters[key_i]['N_correct']=0
			self.letters[key_i]['weight']=1
		# Confidence in ability to learn a new letter without forgetting past material
		self.confidence = 0
		# Confidence threshold above which adding a new letter is suggested
		self.thr = 0.6
	def add_letter(self):
		nya = {}
		for key in braille:
			if key not in self.letters:
				nya[key]=True
		if len(nya)>0:
			nya_keys = nya.keys()
			i = np.random.randint(len(nya_keys))
			key_i = nya_keys[i]
			self.letters[key_i]={}
			self.letters[key_i]['N_attempts']=0
			self.letters[key_i]['N_correct']=0
			self.letters[key_i]['weight']=1
			self.confidence*=float(len(self.letters)-1)/len(self.letters)
		else:
			print "No more letters to be added!"
	def test(self):
		if self.confidence>self.thr:
			self.add_letter()
		keys = self.letters.keys()
		wts = []
		for i in range(len(keys)):
			wts.append(self.letters[keys[i]]['weight'])
		cwts = np.cumsum(wts)
		max_cwt = cwts[-1]
		omega = np.random.rand()*max_cwt
		i=0
		while omega>cwts[i]:
			i+=1
		#i = np.random.randint(len(keys))
		key_i = keys[i]
		braille_sym = braille[key_i]
		display(braille_sym)
		guess = raw_input()
		if self.letters[key_i]['N_correct']>0:
			self.confidence-=float(self.letters[key_i]['N_correct'])/(self.letters[key_i]['N_attempts']*len(keys))
		if guess==key_i:
			print "Correct!"
			self.letters[key_i]['N_correct']+=1
		else:
			print "Correct answer: "+key_i
		self.letters[key_i]['N_attempts']+=1
		self.confidence+=float(self.letters[key_i]['N_correct'])/(self.letters[key_i]['N_attempts']*len(keys))
		# Account for the total number of correct guesses in weighting each letter so that newer letters are practiced after
		#	the number of letters becomes substantial.
		self.letters[key_i]['weight']=np.exp(-self.letters[key_i]['N_correct'])
		
