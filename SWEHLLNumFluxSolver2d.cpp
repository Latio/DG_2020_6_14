#include "SWEHLLNumFluxSolver2d.h"



SWEHLLNumFluxSolver2d::SWEHLLNumFluxSolver2d()
{
}


SWEHLLNumFluxSolver2d::~SWEHLLNumFluxSolver2d()
{
}

void SWEHLLNumFluxSolver2d::numfluxSolver_evaluate(double hmin, double gra, double *nx, double *ny, double *fm, double *fp, double *fluxS,int *Nfp,int *Ne)
{
	//double *nx = nx_;
	//double *ny = ny_;
	//double *fm = fm_;
	//double *fp = fp_;
	//double hmin = hmin_;
	//double gra = gra_;
	//
	//int *TNfp_=meshunion->inneredge_p->Nfp;
	//int *K_ = meshunion->inneredge_p->Ne;


	c_Evaluate(hmin, gra, nx, ny, fm, fp, fluxS,Nfp,Ne);
};