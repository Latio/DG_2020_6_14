#pragma once

#include"BoundaryEdge.h"
#include"InnerEdge.h"
#include"Cell.h"
extern "C" {
	void c_GetMeshIntegralValue(double *nodeVal_, double *wq_, double *J_, double *Vq_, int *Np_, int *K_, int *Nq_, double *integralValue_);
}


class MeshUnion
{
public:
	MeshUnion();
	~MeshUnion();
	void GetMeshAverageValue(double *nodeVal, double *averageValue);


	InnerEdge inneredge;
	static InnerEdge *inneredge_p;

	BoundaryEdge boundarydge;
	static BoundaryEdge *boundarydge_p;

	Cell cell;
	static Cell *cell_p;
	//notice the order of definition,which will cause the of constructor. 

	static int*K;
	static int*Nv;
	static int Nfield;

	static double*charLength;
	static double*EToE;
	static double*EToF;
	static double*EToM;
	static signed char*EToR;
	static double*EToV;
	static double*ind;
	static double*J;
	//static double*K;
	static double*LAV;
	//static double*Nv;
	static double*rx;
	static double*ry;
	static double*rz;
	static signed char *status;
	static double*sx;
	static double*sy;
	static double*sz;
	static double*tx;
	static double*ty;
	static double*type;
	static double*tz;
	static double*vx;
	static double*vy;
	static double*vz;
	static double*x;
	static double*xc;
	static double*y;
	static double*yc;
	static double*z;
	static double*zc;
};


