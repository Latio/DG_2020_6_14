#include "TopographySourceTerm.h"



TopographySourceTerm::TopographySourceTerm()
{
}


TopographySourceTerm::~TopographySourceTerm()
{
}

void TopographySourceTerm::EvaluateTopographySourceTerm(double gra, double *fphys, double *zGrad, double *frhs_temp)
{
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;
	int Nfield = meshunion->Nfield;
	signed char *status = meshunion->status;


	c_EvaluateSourceTopography2d(gra, status, fphys, zGrad, frhs_temp, Np, K, Nfield);
	//double gra_, signed char *status_, double *fphys_, double *zGrad_, double *frhs_temp, int *Np_, int *K_, int Nfield_
}