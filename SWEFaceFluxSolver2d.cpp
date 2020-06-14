#include "SWEFaceFluxSolver2d.h"



SWEFaceFluxSolver2d::SWEFaceFluxSolver2d()
{
}


SWEFaceFluxSolver2d::~SWEFaceFluxSolver2d()
{
}

void SWEFaceFluxSolver2d::surfluxSolver_evaluate(double hmin, double gra, double *nx, double *ny, double *fm, double *fluxM,int *Nfp,int *Ne)
{

	int Nfield = meshunion->Nfield;

	c_EvaluateSurfFlux(hmin, gra, nx, ny, fm, fluxM, Nfp, Ne, Nfield);
};