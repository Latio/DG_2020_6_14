#pragma once
#include"MeshUnion.h"
extern MeshUnion mesh;
extern const MeshUnion *meshunion;

extern "C" {
	void c_EvaluateSourceTopography2d(double gra_, signed char *status_, double *fphys_, double *zGrad_, double *frhs_temp, int *Np_, int *K_, int Nfield_);

}

class TopographySourceTerm
{
public:
	TopographySourceTerm();
	~TopographySourceTerm();
	void EvaluateTopographySourceTerm(double gra, signed char *status, double *fphys, double *zGrad, double *frhs_temp);
};

