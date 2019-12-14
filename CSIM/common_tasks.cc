#include "common_tasks.hh"
void sort0(double * x, int N) {
	double auxmin = x[0];
	int kmin;
	int ind= 0;
	while(ind<N) {
		auxmin = x[ind];
		kmin = ind;
		for(int k=ind+1;k<N;k++) {
			if(x[k]<auxmin) {
				auxmin = x[k];
				kmin = k;
			}
		}
		x[kmin] = x[ind];
		x[ind] = auxmin;
		ind++;
	}
}

void merge_adjacent(double * x, int N1, int N2) {
	double * xtemp = new double [N1+N2];
	int ind1 = 0;
	int ind2 = 0;
	int ind = 0;
	while(ind<N1+N2) {
		while(ind1<N1&&ind2<N2) {
			if(x[ind1]<=x[N1+ind2]) {
				while(x[ind1]<=x[N1+ind2]&&ind1<N1) {
					xtemp[ind++] = x[ind1++];
				}
			}
			else {
				while(x[ind1]>x[N1+ind2]&&ind2<N2) {
					xtemp[ind++] = x[N1+(ind2++)];
				}
			}
		}
		while(ind1<N1) {
			xtemp[ind++] = x[ind1++];
		}
		while(ind2<N2) {
			xtemp[ind++] = x[N1+ind2++];
		}
	}
	for(int i=0;i<N1+N2;i++) x[i]=xtemp[i];
	delete [] xtemp;
}

void sort1(double * x, int N) {
	int blocklen =8;
	int Nblocks =N/blocklen;
	for(int i=0;i<Nblocks-1;i++) {
		sort0(x+blocklen*i,blocklen);
	}
	sort0(x+blocklen*(Nblocks-1),blocklen+N%blocklen);
	int exceptlen = blocklen+N%blocklen;
	while(Nblocks>1) {
		int newblocklen = 2*blocklen;
		int newNblocks = Nblocks/2;
		for(int i=0;i<newNblocks-1;i++) {
			merge_adjacent(x+i*newblocklen,blocklen,blocklen);
		}
		merge_adjacent(x+(newNblocks-1)*newblocklen,blocklen,exceptlen);
		blocklen = newblocklen;
		exceptlen *= 2;
		Nblocks = newNblocks;
	}
}

void sort1(int * x, int N) {
	double * y = new double [N];
	for(int i=0;i<N;i++) y[i]=double(x[i]);
	sort1(y,N);
	for(int i=0;i<N;i++) x[i]=int(y[i]+0.01); // offset to prevent rounding errors
}
