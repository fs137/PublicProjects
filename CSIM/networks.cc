// NOTE: fix 'path_exists' function
// (Latest version of the networks.cc file)
#include "cstdlib"
#include "networks.hh"

// Elementary structures

VertexSet::VertexSet() {
	mem = dmem;
	Nv = 0;
	v = new int [mem];
}
VertexSet::VertexSet(int * v_,int Nv_): v(v_),Nv(Nv_) {
	mem = (int(Nv/dmem)+1)*dmem;
}
VertexSet::VertexSet(const VertexSet& vs) {
	v = new int [vs.mem];
	Nv = vs.Nv;
	mem = vs.mem;
	for(int i=0;i<Nv;i++) v[i]=vs.v[i];
}
VertexSet::~VertexSet() {
	delete [] v;
}

void VertexSet::operator=(const VertexSet& vs) {
	delete [] v;
	v = new int [vs.mem];
	Nv = vs.Nv;
	mem = vs.mem;
	for(int i=0;i<Nv;i++) v[i]=vs.v[i];
}
void add_v(int vv, VertexSet& vs) {
	// check if there is enough memory:
	if(vs.Nv==vs.mem) {
		vs.mem+=dmem;
		int * auxv = new int [vs.mem];
		for(int i=0;i<vs.Nv;i++) auxv[i]=vs.v[i];
		delete [] vs.v;
		vs.v = auxv;
	}
	vs.v[vs.Nv]=vv;
	vs.Nv++;
}
void add_v_safe(int vv, VertexSet& vs) {
	if(vs.Nv==vs.mem) {
		vs.mem+=mem_inc;
		int * auxv = new int [vs.mem];
		for(int i=0;i<vs.Nv;i++) auxv[i]=vs.v[i];
		delete [] vs.v;
		vs.v = auxv;
	}
	// check if vv is already in the list:
	bool found = false;
	for(int i=0;i<vs.Nv&&!found;i++) {
		if(vs.v[i]==vv) found = true;
	}
	if(!found) {
		vs.v[vs.Nv++]=vv;
	}
}
void remove_v(int vv, VertexSet& vs) {
	// look for vv in the list:
	bool found = false;
	int ii;
	for(ii=0;ii<vs.Nv&&!found;ii++) {
		if(vs.v[ii]==vv) {
			found=true;
			break;
		}
	}
	if(ii==vs.Nv) {
		printf("Error in remove_v(int, VertexSet&): vertex does not appear in VertexSet\n");
		return;
	}
	// copy the final element of vs.v to vs.v[ii]:
	vs.Nv--;
	vs.v[ii]=vs.v[vs.Nv];
}

int_decor::int_decor(): v(new int [mem_inc]), mem(mem_inc),Nv(0),Nremoved(0),decor(new int [mem_inc]) {}
int_decor::int_decor(int * v_, int * decor_, int Nv_):Nv(Nv_) {
	Nremoved = 0;
	mem = (int(Nv_/mem_inc)+1)*mem_inc;
	v = new int [mem];
	decor = new int [mem];
	for(int i=0;i<Nv;i++) {
		v[i]=v_[i];
		decor[i]=decor_[i];
	}
}
int_decor::int_decor(const int_decor& idec) {
	Nv = idec.Nv;
	mem = idec.mem;
	Nremoved = idec.Nremoved;
	v = new int [mem];
	decor = new int [mem];
	for(int i=0;i<Nv;i++) {
		v[i]=idec.v[i];
		decor[i]=idec.decor[i];
	}
}
void int_decor::operator=(const int_decor& idec) {
	delete [] v;
	delete [] decor;
	Nv = idec.Nv;
	mem = idec.mem;
	Nremoved = idec.Nremoved;
	v = new int [mem];
	decor = new int [mem];
	for(int i=0;i<Nv;i++) {
		v[i]=idec.v[i];
		decor[i]=idec.decor[i];
	}
}
void int_decor::reset() {
	delete [] v;
	delete [] decor;
	mem = mem_inc;
	Nv = 0;
	Nremoved =0;
	v = new int [mem];
	decor = new int [mem];
}

void add_v(int vv,int vvid, int_decor& idec, int& found_at) {
	// check if vv is already in the list:
	found_at = -1;
	for(int i=0;i<idec.Nv&&found_at==-1;i++) {
		if(idec.v[i]==vv) found_at = i;
	}
	if(found_at==-1) {
		// check if there is enough memory:
		if(idec.Nv==idec.mem) {
			int * auxv = new int [idec.mem];
			int * auxdecor = new int [idec.mem];
			for(int i=0;i<idec.mem;i++) { 
				auxv[i]=idec.v[i];
				auxdecor[i]=idec.decor[i];
			}
			delete [] idec.v;
			delete [] idec.decor;
			idec.mem+=mem_inc;
			idec.v = new int [idec.mem];
			idec.decor = new int [idec.mem];
			for(int i=0;i<idec.Nv;i++) {
				idec.v[i]=auxv[i];
				idec.decor[i]=auxdecor[i];
			}
			delete [] auxv;
			delete [] auxdecor;
		}
		idec.v[idec.Nv]=vv;
		idec.decor[idec.Nv++]=vvid;
	}
}

void remove_v(int vv, int_decor& idec) {
	// look for vv in the list:
	bool found = false;
	int ii;
	for(ii=0;ii<idec.Nv&&!found;ii++) {
		if(idec.v[ii]==vv) found=true;
	}
	idec.Nv--;
	for(int i=ii;i<idec.Nv;i++) {
		idec.v[i]=idec.v[i+1];
		idec.decor[i]=idec.decor[i+1];
	}
	idec.Nremoved++;
	if(idec.Nremoved>mem_inc) {
		int * auxv = new int [idec.Nv];
		int * auxdecor = new int [idec.Nv];
		for(int i=0;i<idec.Nv;i++) {
			auxv[i]=idec.v[i];
			auxdecor[i]=idec.decor[i];
		}
		delete [] idec.v;
		delete [] idec.decor;
		idec.mem-=mem_inc;
		idec.v = new int [idec.mem];
		for(int i=0;i<idec.Nv;i++) {
			idec.v[i]=auxv[i];
			idec.decor[i]=auxdecor[i];
		}
		delete [] auxv;
		delete [] auxdecor;
		idec.Nremoved = 0;
	}
}

// Graph Structures
Graph::Graph() {
	Nv = 0;
	Nremoved = 0;
	N_edges = 0;
	mem = mem_inc;
	id = new int [mem];
	nbors = new int * [mem];
	nbors_index = new int*[mem]; 
	edge_dec = new int*[mem];
	Nnbors = new int [mem];
	nbors_mem = new int [mem];
	for(int i=0;i<mem;i++) {
		Nnbors[i]=0;
		nbors[i] = new int [dmem];
		nbors_index[i]=new int[dmem];
		nbors_mem[i]=dmem;
		edge_dec[i] = new int[dmem];
	}
	total_mem = 3*mem*dmem+3*mem+4;
}
Graph::Graph(int Nv_, int * id_,int * Nnbors_, int ** nbors_) {
	Nv = Nv_;
	N_edges = 0;
	mem = mem_inc*(int(Nv_/mem_inc)+1);
	Nnbors = new int [mem];
	nbors = new int * [mem];
	nbors_index = new int*[mem];
	edge_dec = new int *[mem];
	nbors_mem = new int [mem];
	id = new int [mem];
	total_mem = (mem-Nv)*dmem+3*mem+4;
	for(int i=Nv;i<mem;i++) {
		nbors[i]=new int [dmem];
		nbors_index[i] = new int[dmem];
		edge_dec[i]=new int[dmem];
		nbors_mem[i]=dmem;
	}
	for(int i=0;i<Nv;i++) {
		if(Nnbors_!=NULL) {
			Nnbors[i]=Nnbors_[i];
			N_edges+=Nnbors[i];
			nbors_mem[i]=dmem*(int(Nnbors[i]/dmem)+1);
			nbors[i]=new int[nbors_mem[i]];
			nbors_index[i]=new int[nbors_mem[i]];
			edge_dec[i] = new int[nbors_mem[i]];
		}
		else {
			Nnbors[i]=0;
			nbors[i]=new int[dmem];
			nbors_index[i]=new int[dmem];
			edge_dec[i] = new int[dmem];
			nbors_mem[i]=dmem;
		}
		total_mem+=3*nbors_mem[i];
		if(nbors_!=NULL) {
			for(int j=0;j<Nnbors[i];j++) {
				nbors[i][j]=nbors_[i][j];
				int np = nbors_[i][j];
				int index_i_np = 0;
				for(index_i_np=0;index_i_np<Nnbors[np];index_i_np++) {
					if(nbors_[np][index_i_np]==i) {
						nbors_index[i][j]=index_i_np;
						break;
					}
				}
			}
		}
	}
	if(Nnbors_!=NULL) {
		if(N_edges%2!=0) printf("Error: total sum of number of neighbors should be even\n");
		N_edges/=2;
	}
	if(id_!=NULL) {
		for(int i=0;i<Nv;i++) {
			id[i]=id_[i];
		}
	} else {
		for(int i=0;i<Nv;i++) {
			id[i]=i;
		}
	}
	Nremoved = 0;
}
Graph::Graph(const Graph& g) {
	total_mem = g.total_mem;
	Nv = g.Nv;
	N_edges = g.N_edges;
	mem = g.mem;
	Nnbors = new int [mem];
	nbors = new int * [mem];
	nbors_index = new int*[mem];
	edge_dec = new int*[mem];
	nbors_mem = new int [mem];
	id = new int [mem];
	for(int i=0;i<mem;i++) {
		nbors_mem[i]=g.nbors_mem[i];
		nbors[i]= new int [nbors_mem[i]];
		nbors_index[i] = new int[nbors_mem[i]];
		edge_dec[i]=new int[nbors_mem[i]];
	}
	for(int i=0;i<Nv;i++) {
		Nnbors[i]=g.Nnbors[i];
		id[i] = g.id[i];
		for(int j=0;j<Nnbors[i];j++) {
			nbors[i][j]=g.nbors[i][j];
			nbors_index[i][j] = g.nbors_index[i][j];
			edge_dec[i][j]=g.edge_dec[i][j];
		}
	}
	Nremoved = g.Nremoved;
	for(int i=0;i<Nremoved;i++) rem_list[i]=g.rem_list[i];
}
void Graph::operator=(const Graph& g) {
	if(Nv>0) {
		for(int i=0;i<mem;i++) {
			delete [] nbors[i];
			delete [] nbors_index[i];
			delete [] edge_dec[i];
		}
		delete [] nbors;
		delete [] nbors_index;
		delete [] edge_dec;
		delete [] nbors_mem;
		delete [] Nnbors;
		delete [] id;
	}
	total_mem = g.total_mem;
	Nv = g.Nv;
	N_edges = g.N_edges;
	mem = g.mem;
	Nnbors = new int [mem];
	nbors = new int * [mem];
	nbors_index = new int*[mem];
	edge_dec = new int*[mem];
	nbors_mem = new int [mem];
	for(int i=Nv;i<mem;i++) {
		nbors[i]=new int [dmem];
		nbors_index[i] = new int[dmem];
		nbors_mem[i]=dmem;
		edge_dec[i]=new int[dmem];
	}
	id = new int [mem];
	for(int i=0;i<Nv;i++) {
		Nnbors[i]=g.Nnbors[i];
		id[i]=g.id[i];
		nbors_mem[i]=g.nbors_mem[i];
		nbors[i]=new int [nbors_mem[i]];
		nbors_index[i] = new int[nbors_mem[i]];
		edge_dec[i] = new int[nbors_mem[i]];
		for(int j=0;j<Nnbors[i];j++) {
			nbors[i][j]=g.nbors[i][j];
			nbors_index[i][j]=g.nbors_index[i][j];
			edge_dec[i][j]=g.edge_dec[i][j];
		}
	}
	Nremoved = g.Nremoved;
	for(int i=0;i<Nremoved;i++) rem_list[i]=g.rem_list[i];
}
Graph::~Graph() {
	for(int i=0;i<mem;i++) {
		if(nbors_mem[i]>0) {
			delete [] nbors[i];
			delete [] edge_dec[i];
			delete [] nbors_index[i];
		}
		else printf("Error: nbors_mem should always be positive\n");
	}
	delete [] nbors;
	delete [] edge_dec;
	delete [] nbors_index;
	delete [] Nnbors;
	delete [] id;
}

// Check if an id appears in a graph:
bool new_id(const int& id_,const Graph& g) {
	bool nf = true; // nf: 'not found', or 'new feature'
	for(int i=0;i<g.Nv&&nf;i++) {
		if(g.id[i]==id_) nf=false;
	}	
	return nf;
}

// Graph operations: adding or removing nodes and edges:
void remove_e(const int& in1,const int& in2, Graph& g, int mode) {
	int v1,v2;
	int iv2_in_v1,iv1_in_v2;
	if(mode==0) { // in1 and in2 are array indices for vertices.
		v1 = in1;
		v2 = in2;
		if(g.Nnbors[v1]<g.Nnbors[v2]) {
			iv2_in_v1=0;
			while(g.nbors[v1][iv2_in_v1]!=v2) {
				iv2_in_v1+=1;
			}
			iv1_in_v2=g.nbors_index[v1][iv2_in_v1];
		} else {
			iv1_in_v2=0;
			while(g.nbors[v2][iv1_in_v2]!=v1) {
				iv1_in_v2+=1;
			}
			iv2_in_v1=g.nbors_index[v2][iv1_in_v2];
		}
	} 
	else if(mode==12321) {
		// in1 is an array index for vertex 1, and in2 is a neighbor list index
		v1 = in1;
		v2 = g.nbors[v1][in2];
		iv2_in_v1=in2;
		iv1_in_v2=g.nbors_index[v1][in2];
	}
	else { // in1 and in2 are the ids of vertices: need to find associated array indices.
		bool v1_tbf=true;
		bool v2_tbf=true;
		int ni = 0;
		while(v1_tbf||v2_tbf) {
			if((g.id[ni]!=in1)&&(g.id[ni]!=in2)){}
			else if(g.id[ni]==in1) {
				v1=ni;
				v1_tbf=false;
			} else {
				v2=ni;
				v2_tbf=false;
			}
		}
		if(g.Nnbors[v1]<g.Nnbors[v2]) {
			iv2_in_v1=0;
			while(g.nbors[v1][iv2_in_v1]!=v2) iv2_in_v1+=1;
			iv1_in_v2=g.nbors_index[v1][iv2_in_v1];
		} else {
			iv1_in_v2=0;
			while(g.nbors[v2][iv1_in_v2]!=v1) iv1_in_v2+=1;
			iv2_in_v1=g.nbors_index[v2][iv1_in_v2];
		}
	}
	// find v2 in g.nbors[v1] and remove it by copying the last element of g.nbors[v1] over it.
	int j=0;
	while(j<g.Nnbors[v1]) {
		if(g.nbors[v1][j]==v2) break;
		j++;
	}
	g.Nnbors[v1]--;
	// Update the neighbor lists
	g.nbors[v1][iv2_in_v1]=g.nbors[v1][g.Nnbors[v1]];
	g.Nnbors[v2]-=1;
	g.nbors[v2][iv1_in_v2]=g.nbors[v2][g.Nnbors[v2]];
	// Update edge decorations
	g.edge_dec[v1][iv2_in_v1]=g.edge_dec[v1][g.Nnbors[v1]];
	g.edge_dec[v2][iv1_in_v2]=g.edge_dec[v2][g.Nnbors[v2]];
	// update nbors_index
	int v1p = g.nbors[v1][iv2_in_v1];
	int iv1_in_v1p = g.nbors_index[v1][g.Nnbors[v1]];
	g.nbors_index[v1][iv2_in_v1]=iv1_in_v1p;
	g.nbors_index[v1p][iv1_in_v1p] = iv2_in_v1;
	int v2p = g.nbors[v2][iv1_in_v2];
	int iv2_in_v2p = g.nbors_index[v2][g.Nnbors[v2]];
	g.nbors_index[v2][iv1_in_v2]=iv2_in_v2p;
	g.nbors_index[v2p][iv2_in_v2p] = iv1_in_v2;
	g.N_edges-=1;
}

void remove_v(const int& in, Graph& g,int mode) {	
	int v;
	if(mode==0) {
		v=in;
		if(v<g.Nv) {}
		else {
			printf("Attempting to remove nonexistent vertex\n");
			return;
		}
	}
	else { // in = id of particle to be removed.
		v=0;
		while(v<g.Nv) {
			if(g.id[v]==in) break;
			v++;
		}
		if(v<g.Nv) {}
		else {
			printf("Attempting to remove nonexistent vertex\n");
			return;
		}
	}
	g.N_edges-=g.Nnbors[v];
	// Remove the vertex from the neighbor lists of other vertices.
	for(int i=0;i<g.Nnbors[v];i++) {
		remove_e(v,i,g,12321);
	}
	// Copy the data for the last vertex onto the removed vertex:
	g.Nv--; 
	if(v<g.Nv) {
		g.Nnbors[v]=g.Nnbors[g.Nv];
		// swap pointers of v and the last vertex:
		int * nbors_temp = g.nbors[g.Nv];
		g.nbors[g.Nv]=g.nbors[v];
		g.nbors[v]=nbors_temp;
		int * nbors_index_temp = g.nbors_index[g.Nv];
		g.nbors_index[g.Nv]=g.nbors_index[v];
		g.nbors_index[v]=nbors_index_temp;
		int * edge_dec_temp = g.edge_dec[g.Nv];
		g.edge_dec[g.Nv]=g.edge_dec[v];
		g.edge_dec[v]=edge_dec_temp;
		int temp_nbors_mem = g.nbors_mem[v];
		g.nbors_mem[v]=g.nbors_mem[g.Nv];
		g.nbors_mem[g.Nv]=temp_nbors_mem;
		g.Nnbors[g.Nv]=0;
		// change (-,Nv) edges to (-,v) edges:	
		for(int i=0;i<g.Nnbors[v];i++) {
			int vi = g.nbors[v][i];
			int index_v_vi = g.nbors_index[v][i];
			g.nbors[vi][index_v_vi] = v;
		}
		//	change the id at site v to that at site g.Nv:
		g.id[v]=g.id[g.Nv];
	} else {
		g.Nnbors[g.Nv]=0;
	}
}

void add_e(Graph& g, const int& in1, const int& in2,int mode) {
	// in1, int2 are either vertex indices or vertex identities within g.
	// mode 0: vertex indices (add edge between vertex n1=in1 and vertex n2=in2.)
	// mode 1: vertex ids (add edge between vertices n1 and n2, where g.id[n1]=in1 and g.id[n2]=in2.)
	int n1,n2;
	if(mode==0) {
		n1 = in1;
		n2 = in2;
		if(n1<g.Nv&&n2<g.Nv) {}
		else {
			printf("Error: attempting to connect nonexistent vertices. n1=%d, n2=%d, g.Nv=%d\n",n1,n2,g.Nv);
			return;
		}
		if(n1!=n2){}
		else {
			printf("Warning: Graph class was not designed with explicit support for self-edges.\n");
		}
	}
	else {
		// find the ids of the particles: 
		for(n1=0;n1<g.Nv;n1++){
			if(g.id[n1]==in1||g.id[n1]==in2) break;
		}
		if(g.id[n1]==in1) {
			for(n2=n1+1;n2<g.Nv;n2++) {
				if(g.id[n2]==in2) break;
			}
		} else if(g.id[n1]==in2) {
			for(n2=n1+1;n2<g.Nv&&g.id[n2]!=in1;n2++){
				if(g.id[n2]==in1) break;
			}
		} else {
			return;
		}
		if(n1==g.Nv||n2==g.Nv) {
			// these ids do not occur in g.
			return;
		}
	}
	// check if there already exists an edge between n1 and n2:
	for(int i=0;i<g.Nnbors[n1];i++) {
		if(g.nbors[n1][i]==n2) {
			return;
		}
	}
	g.N_edges+=1;
	// check if more memory needs to be allocated:
	if(g.Nnbors[n1]==g.nbors_mem[n1]) {
		g.total_mem+=2*g.nbors_mem[n1];
		g.nbors_mem[n1]*=2;
		int * temp_nbors = new int [g.nbors_mem[n1]];
		int * temp_edge_dec = new int [g.nbors_mem[n1]];
		int * temp_nbors_index = new int[g.nbors_mem[n1]];
		for(int i=0;i<g.Nnbors[n1];i++) {
			temp_nbors[i]=g.nbors[n1][i];
			temp_edge_dec[i]=g.edge_dec[n1][i];
			temp_nbors_index[i]=g.nbors_index[n1][i];
		}
		delete [] g.nbors[n1];
		delete [] g.edge_dec[n1];
		delete [] g.nbors_index[n1];
		g.nbors[n1]=temp_nbors;
		g.edge_dec[n1] = temp_edge_dec;
		g.nbors_index[n1] = temp_nbors_index;
	}
	if(g.Nnbors[n2]==g.nbors_mem[n2]) {
		g.total_mem+=2*g.nbors_mem[n2];
		g.nbors_mem[n2]*=2;
		int * temp_nbors = new int [g.nbors_mem[n2]];
		int * temp_edge_dec = new int[g.nbors_mem[n2]];
		int * temp_nbors_index = new int[g.nbors_mem[n2]];
		for(int i=0;i<g.Nnbors[n2];i++) {
			temp_nbors[i]=g.nbors[n2][i];
			temp_edge_dec[i] = g.edge_dec[n2][i];
			temp_nbors_index[i] = g.nbors_index[n2][i];
		}
		delete [] g.nbors[n2];
		delete [] g.edge_dec[n2];
		delete [] g.nbors_index[n2];
		g.nbors[n2]=temp_nbors;
		g.edge_dec[n2] = temp_edge_dec;
		g.nbors_index[n2]=temp_nbors_index;
	}
	g.nbors[n1][g.Nnbors[n1]]=n2;
	g.nbors[n2][g.Nnbors[n2]]=n1;
	g.nbors_index[n1][g.Nnbors[n1]]=g.Nnbors[n2];
	g.nbors_index[n2][g.Nnbors[n2]]=g.Nnbors[n1];
	g.Nnbors[n2]+=1;
	g.Nnbors[n1]+=1;
}

void add_v(Graph& g, const int& id, int * nvs, int Nnvs) {
	// trivially add a new vertex to the graph (optionally: connect it to a list of vertices).
	if(g.Nv == g.mem) {
		// increase the memory allocated to g
		g.mem+=mem_inc;
		int * temp_Nnbors= new int[g.mem]; 
		int ** temp_nbors = new int*[g.mem];
		int ** temp_nbors_index = new int*[g.mem]; 
		int * temp_nbors_mem = new int [g.mem];
		int * temp_id = new int [g.mem];
		int ** temp_edge_dec = new int*[g.mem];
		for(int i=g.Nv;i<g.mem;i++) {
			temp_nbors[i]=new int[dmem];
			temp_nbors_mem[i] = dmem;
			temp_edge_dec[i]=new int[dmem];
			temp_nbors_index[i]=new int[dmem];
		}
		for(int i=0;i<g.Nv;i++) {
			temp_nbors_mem[i]=g.nbors_mem[i];
			temp_edge_dec[i] = new int[temp_nbors_mem[i]];
			temp_nbors[i]=new int [temp_nbors_mem[i]];
			temp_nbors_index[i]=new int[temp_nbors_mem[i]];
			temp_Nnbors[i]=g.Nnbors[i];
			for(int j=0;j<g.Nnbors[i];j++) {
				temp_nbors[i][j]=g.nbors[i][j];
				temp_nbors_index[i][j]=g.nbors_index[i][j];
				temp_edge_dec[i][j] = g.edge_dec[i][j];
			}
			delete [] g.nbors[i];
			delete [] g.nbors_index[i];
			delete [] g.edge_dec[i];
			temp_id[i]=g.id[i];
		}
		delete [] g.edge_dec;
		delete [] g.nbors;
		delete [] g.nbors_index;
		delete [] g.nbors_mem;
		delete [] g.Nnbors;
		delete [] g.id;
	
		g.nbors = temp_nbors;
		g.nbors_index = temp_nbors_index;
		g.edge_dec = temp_edge_dec;
		g.id = temp_id;
		g.Nnbors = temp_Nnbors;
		g.nbors_mem = temp_nbors_mem;
		g.total_mem+=2*(mem_inc-1)*dmem+3*mem_inc;
	} 
	if(g.nbors_mem[g.Nv]<=Nnvs) {
		delete [] g.nbors[g.Nv];
		delete [] g.nbors_index[g.Nv];
		delete [] g.edge_dec[g.Nv];
		g.nbors_mem[g.Nv]=dmem*(int(Nnvs/dmem)+1);
		g.total_mem+=2*g.nbors_mem[g.Nv];
		g.nbors[g.Nv] = new int [g.nbors_mem[g.Nv]];
		g.nbors_index[g.Nv]= new int[g.nbors_mem[g.Nv]];
		g.edge_dec[g.Nv] = new int[g.nbors_mem[g.Nv]];
	}
	g.Nnbors[g.Nv]=0;
	g.id[g.Nv]=id;
	if(nvs!=NULL) {
		for(int i=0;i<Nnvs;i++) {
			add_e(g,g.Nv,nvs[i],0);
		}
	}
	g.Nv++;
}

void get_cc(const Graph& g,VertexSet *& g_cc,int& ncc, int& mem,int *& g2g_cc) {
	g2g_cc = new int [g.Nv];
	ncc = 0;
	mem = 1;
	g_cc = new VertexSet[mem];
	bool *visited=new bool[g.Nv];
	for(int i=0;i<g.Nv;i++) visited[i]=false;
	int start_vertex = 0;
	while(start_vertex<g.Nv) {
		// look for an unvisited site in g:
		while(visited[start_vertex]&&(start_vertex<g.Nv)) start_vertex++;
		if(start_vertex>=g.Nv) {
			break;
		}
		VertexSet bdry_cc;
		// check if more memory is needed:
		if(ncc==mem) {
			mem+=1;
			VertexSet * aux_g_cc = new VertexSet [mem];
			for(int i=0;i<ncc;i++) {
				aux_g_cc[i]=g_cc[i];
			}
			delete [] g_cc;
			g_cc = aux_g_cc;
		}
		// add start_vertex to bdry_cc[ncc], 'decorated' by its index in g_cc[ncc]
		add_v(start_vertex,bdry_cc);
		if(start_vertex>=g.Nv) printf("Error: start_vertex=%d, g.Nv=%d",start_vertex,g.Nv);
		g2g_cc[start_vertex]=0;
		add_v(start_vertex,g_cc[ncc]);
		visited[start_vertex]=true;
		bool exp_possible =true;
		while(exp_possible) {
			exp_possible=false;
			// for each site on bdry_cc[ncc], add neighboring sites that do not also appear in g_cc[ncc]:
			VertexSet new_bdry_cc;
			for(int i=0;i<bdry_cc.Nv;i++) {
				int vi = bdry_cc.v[i];
				for(int j=0;j<g.Nnbors[vi];j++) {
					// check if g.nbors[vi][j] has been visited:
					int vj = g.nbors[vi][j];
					if(visited[vj]) {
						continue;
					}
					else {
						// add g.nbors[vi][j] to new_bdry_cc:
						add_v(vj,new_bdry_cc);
						g2g_cc[vj]=g_cc[ncc].Nv;
						add_v(vj,g_cc[ncc]);
						visited[vj]=true;
						exp_possible=true;
					}
				}
			}
			bdry_cc = new_bdry_cc;
		}
		ncc++;
	}
	delete [] visited;
}

bool path_exists(const Graph &g,const int &in1, const int & in2, int mode) {
	if(in1==in2) return true;
	// determine vertices:
	int v1,v2;
	if(mode==0) {
		// check if in1 and in2 are within graph:
		if(in1<g.Nv||in2<g.Nv) {
			v1 = in1;
			v2 = in2;
		}
		else {
			printf("Error in path_exists: specified vertices do not exist in graph.\n");
			return false;
		}
	}
	else {
		// find vertices corresponding to ids in1 and in2:
		v1 = 0;
		while(v1<g.Nv) {
			if(g.id[v1]==in1||g.id[v1]==in2) {
				break;
			}
			v1++;
		}
		if(v1==g.Nv) return false;
		if(g.id[v1]==in1) {
			v2 = v1+1;
			while(v2<g.Nv) {
				if(g.id[v2]==in2) break;
				v2++;
			}
		}
		else if(g.id[v1]==in2) {
			v2 = v1+1;
			while(v2<g.Nv) {
				if(g.id[v2]==in1) break;
			}
		}
		else {
			printf("Error in path_exists\n");
			return false;
		}
		if(v2==g.Nv) return false;
	}
	// Check if v2 is contained in the component of v1.
	bool* visited=new bool[g.Nv];
	for(int i=0;i<g.Nv;i++) visited[i]=false;
	bool exp_possible = true;
	VertexSet bdry_cc;
	add_v(v1,bdry_cc);
	visited[v1]=true;
	while(exp_possible) {
		VertexSet new_bdry_cc;
		exp_possible = false;
		for(int i=0;i<bdry_cc.Nv;i++) {
			int vi = bdry_cc.v[i];
			for(int j=0;j<g.Nnbors[vi];j++) {
				if(visited[g.nbors[vi][j]]) continue;
				else {
					exp_possible = true;
					add_v(g.nbors[vi][j],new_bdry_cc);
					visited[g.nbors[vi][j]]=true;	
				}
			}
		}
		if(visited[v2]) {
			delete [] visited;
			return true;
		}
		bdry_cc = new_bdry_cc;
	}
	if(visited[v2]) {
		return true;
	}
	else return false;
}

