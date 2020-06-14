#pragma once
#include"MeshUnion.h"
extern MeshUnion mesh;
extern const MeshUnion *meshunion;

extern "C" {
	void c_Evaluate(double hmin_, double gra_, double *nx_, double *ny_, double *fm_, double *fp_, double *fluxS_, int *TNfp_, int *K_);
}//hmin, gra, nx, ny, fm, fp, fluxS_

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class SWEHLLNumFluxSolver2d
{
public:
	SWEHLLNumFluxSolver2d();
	~SWEHLLNumFluxSolver2d();

	void numfluxSolver_evaluate(double hmin, double gra, double *nx, double *ny, double *fm, double *fp, double *fluxS, int *Nfp, int *Ne);
};

