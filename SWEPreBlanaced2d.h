#pragma once
#include "SWEConventional2d.h"
#include"SWEPrebalanceVolumeFlux2d.h"

class SWEPreBlanaced2d :
	public SWEConventional2d
{
public:
	SWEPreBlanaced2d();
	~SWEPreBlanaced2d();
	void EvaluateFlux(double *fphys, double *E, double *G);
	SWEPrebalanceVolumeFlux2d sweprebalancevolumeflux2d;

};

