#pragma once
#include<iostream>
#include"MeshUnion.h"
#include"cblas.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

extern "C" {
	void c_EvaluateSurfValue(double *fm_, double  *fp_, double *FToE_, double *FToN1_, double *FToN2_, double *Kappa_, int Nfp_, int Ne_, int Np_, int K_);

}

class AbstractDiffSolver
{
public:
	AbstractDiffSolver();
	~AbstractDiffSolver();

	void EvaluateSurfValue(double *fm_, double  *fp_, double *FToE_, double *FToN1_, double *FToN2_, double *Kappa_, int Nfp_, int Ne_, int Np_, int K_);
	void Evaluate_rdhuv(double *rs, double *d, double *rdhuv, double *temp);

	double *nv;
};


