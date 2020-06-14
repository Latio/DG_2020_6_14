#pragma once
#include"MeshUnion_dim.h"

class Cell
{
public:
	Cell();
	~Cell();

	static int* Nq;
	static int* Nv;

	static double* Dr;
	static double* Ds;
	static double* Dt;
	static double* FToV;
	static double* faceType;
	static double* Fmask;
	static double* invM;
	static double* LAV;
	static double* M;
	static int* N;
	static double* Nface;
	static int* Nfp;
	static double* Nfv;
	static int* Np;
	//static double* Nq;
	//static double* Nv;
	static double* r;
	static double* rq;
	static double* s;
	static double* sq;
	static double* t;
	static int* TNfp;
	static double* tq;
	static double* type;
	static double* V;
	static double* Vq;
	static double* vr;
	static double* vs;
	static double* vt;
	static double* wq;
};
