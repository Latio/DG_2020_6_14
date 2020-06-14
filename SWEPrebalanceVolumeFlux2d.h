#pragma once
#include"MeshUnion.h"
extern MeshUnion mesh;
extern const MeshUnion *meshunion;

extern "C" {
	void c_EvaluateFlux2d(double hmin_, double gra_, signed char *status_, double *fphys_, int *Np_, int *K_, double *E_, double *G_);
}
class SWEPrebalanceVolumeFlux2d
{
public:
	SWEPrebalanceVolumeFlux2d();
	~SWEPrebalanceVolumeFlux2d();
	void evaluate(double hmin, double gra, double *fphys, double *E, double *G);
};

