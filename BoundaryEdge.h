#pragma once
#include"Bcell.h"
extern "C" {
	void c_boundary_EvaluateSurfValue(const double *FToM_, const double *FToE_, const double *FToN1_, const double *FToN2_, double *fphys_, double *fm_, double *fp_, const int *Nfp_, const  int *Ne_, const  int *Np_, const  int *K_, const  int Nfield_);
	void c_boundary_EvaluateStrongFormEdgeRHS(double *invM_, double *M_, double *FToE_, double *FToN1_, double *Js_, double *J_, double *fluxM__, double *fluxS__, int *Np_, int *K_, int *Nfp_, int *Ne_, int Nfield_, double *frhs_temp_);
}

class BoundaryEdge
{
public:
	BoundaryEdge();
	~BoundaryEdge();
	void EvaluateSurfValue(double *fphys, double *fm, double *fp, int *Np, int *K, int Nfield);
	void EvaluateStrongFromEdgeRHS(double *invM_, double *J, double *fluxM, double *fluxS, int *Np, int *K, int Nfield, double *const frhs_temp);

	Bcell bcell;

	static double* FToE;
	static double* FToF;
	static double* FToM;
	static double* FToN1;
	static double* FToN2;
	static double* FToV;
	static signed char* ftype;
	static double* Js;
	static double* LAV;
	static double* M;
	static int* Ne;
	static int* Nfp;
	static double* nx;
	static double* ny;
	static double* nz;
	static double* r;
	static double* xb;
	static double* yb;
};

