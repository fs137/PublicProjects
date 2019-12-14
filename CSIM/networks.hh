#ifndef NETWORKS_HH
#define NETWORKS_HH
#include "stddef.h"
#include "stdio.h"
#include "common_tasks.hh"
const int MaxNremoved = 128;
const int mem_inc = 256;
const int dmem = 16; // 

// Elementary structures

class VertexSet {
	public:
	int * v;
	int mem;
	int Nv;
	VertexSet();
	VertexSet(int * v_,int Nv_);
	VertexSet(const VertexSet& vs); 
	~VertexSet(); 
	void operator=(const VertexSet& vs); 
};
void add_v(int vv, VertexSet& vs);
void remove_v(int vv, VertexSet& vs);

class int_decor {
	public:
	int * v;
	int mem;
	int Nv;
	int Nremoved;
	int * decor;
	int_decor();
	int_decor(int * v_, int * decor_, int Nv_);
	int_decor(const int_decor& idec);
	void operator=(const int_decor& idec);
	void reset(); 
};

void add_v(int vv,int vvid, int_decor& idec, int& found_at);
void remove_v(int vv, int_decor& idec); 
// Graph Structures
class Graph {
	private:
	// a list of nodes that have been removed:
	public:
	int total_mem;
	int Nremoved;
	int rem_list [MaxNremoved];
	int mem; // NOTE: use this to avoid frequent shifts of memory.
	int Nv; // number of vertices;
	int N_edges; // number of edges
	int * Nnbors; // number of neighbors
	int ** nbors; // neighbor list
	int ** nbors_index; // indices of each vertex in the neighbor lists of its neighbors.
	int ** edge_dec; // An integer property for each edge
	int * nbors_mem; // memory for neighbor list
	int * id; // integer labels for each vertex.
	Graph();
	Graph(int Nv_, int * id_=NULL,int * Nnbors_=NULL, int ** nbors_=NULL);
	Graph(const Graph& g); 
	void operator=(const Graph& g); 
	~Graph();
};

// Check if an id appears in a graph:
bool new_id(const int& id_,const Graph& g);

// Graph operations: adding or removing a node:
// Problems: 
// 	(i) g.Nnbors[v] eventually becomes negative for some vertices.
// 	(ii) 
void remove_v(const int& in, Graph& g,int mode=0); 
void remove_e(const int& in1,const int& in2, Graph& g, int mode);

void add_e(Graph& g, const int& in1, const int& in2,int mode=0);
void add_v(Graph& g, const int& id, int * nvs=NULL, int Nnvs=0);

void get_cc(const Graph& g,VertexSet *& g_cc, int& ncc, int& mem,int*& g2g_cc);

bool path_exists(const Graph& g, const int& v1, const int& v2, int mode=0);

#endif
