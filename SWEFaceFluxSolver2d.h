#pragma once
#include"MeshUnion.h"
extern MeshUnion mesh;
extern const MeshUnion *meshunion;


extern "C" {
	void c_EvaluateSurfFlux(double hmin_, double gra_, double *nx_, double *ny_, double *fm_, double *fluxM_, int *fm_Nfp, int *fm_Ne, int fm_Nfield);
}

class SWEFaceFluxSolver2d
{
public:
	SWEFaceFluxSolver2d();
	~SWEFaceFluxSolver2d();

	void surfluxSolver_evaluate(double hmin, double gra, double *nx, double *ny, double *fm, double *fluxM, int *Nfp, int *Ne);
};

