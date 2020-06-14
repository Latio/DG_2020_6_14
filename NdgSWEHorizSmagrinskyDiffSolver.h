#pragma once

#include"NdgHorizDiffSolver.h"



class NdgSWEHorizSmagrinskyDiffSolver :public NdgHorizDiffSolver
{

public:
	NdgSWEHorizSmagrinskyDiffSolver(double c_);
	~NdgSWEHorizSmagrinskyDiffSolver();

	void EvaluateDiffRHS(double *fphys, double *frhs, double *InnerEdgefm2d, double * InnerEdgefp2d, double * BoundaryEdgefm2d, double * BoundaryEdgefp2d);

	void UpdateViscosity(double *hu, double *hv, double *h);
	//void Evaluate_rdhuv(double *r, double *d, double *rdhuv, double *temp);


	double C;

};

