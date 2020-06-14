#pragma once
#include"Icell.h"
extern "C"
{
	void c_inner_EvaluateSurfValue(double *FToE_, double *FToN1_, double *FToN2_, double *fphys_, double *fm_, double *fp_, int *Nfp_, int *Ne_, int *Np_, int *K_, int Nfield_);
	void c_inner_EvaluateStrongFromEdgeRHS(double *invM_, double *M_, double *FToE_, double *FToN1_, double *FToN2_, double *Js_, double *J_, double *fluxM_, double *fluxP_, double *fluxS_, double *frhs_, int *Np_, int *K_, int *Nfp_, int *Ne_, int Nfield_);

};

class InnerEdge
{
public:
	InnerEdge();
	~InnerEdge();

	void EvaluateSurfValue(double *fphys, double *fm, double *fp, int *Np, int *K, int Nfield);
	void EvaluateStrongFromEdgeRHS(double *fluxM, double *fluxP, double *fluxS, double *frhs, double *invM, double *J, int *Np, int *K,int Nfield);

	Icell icell;

	static double* FToE;
	static double* FToF;
	static double* FToM;
	static double* FToN1;
	static double* FToN2;
	static double* FToV;
	static double* Js;
	static double* LAV;
	static double* M;
	static int* Ne;
	static int* Nfp;
	static double* nx;
	static double* ny;
	static double* nz;
	static double* r;

};

