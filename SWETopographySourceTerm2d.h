#pragma once
#include"MeshUnion.h"
#include"cblas.h"

extern const MeshUnion *meshunion;
extern MeshUnion mesh;

extern "C" {
	void c_EvaluateSourceTopography2d(double gra_, signed char *status_, double *fphys_, double *zGrad_, double *frhs_temp, int *Np_, int *K_, int Nfield_);
}

class SWETopographySourceTerm2d
{
public:
	SWETopographySourceTerm2d();
	~SWETopographySourceTerm2d();
	void EvaluateTopographySourceTerm(double gra, double *fphys, double *zGrad, double *frhs);
};

