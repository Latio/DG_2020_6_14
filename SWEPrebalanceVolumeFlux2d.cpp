#include "SWEPrebalanceVolumeFlux2d.h"



SWEPrebalanceVolumeFlux2d::SWEPrebalanceVolumeFlux2d()
{
}


SWEPrebalanceVolumeFlux2d::~SWEPrebalanceVolumeFlux2d()
{
}

void SWEPrebalanceVolumeFlux2d::evaluate(double hmin, double gra, double *fphys, double *E, double *G)
{
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;
	signed char *status = meshunion->status;
	c_EvaluateFlux2d(hmin, gra, status, fphys, Np, K, E, G);
};