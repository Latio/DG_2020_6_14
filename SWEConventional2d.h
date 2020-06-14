#pragma once
#include "SWEAbstract2d.h"
extern "C" {
	void c_EvaluatePostFunc2d(double hmin_, double *fphys_, double *hc_, double *qxc_, double *qyc_, int *K_, int *Np_);
	void c_UpdateWDWetDryState(double hmin_, double *fphys_, signed char *status_, int *Np_, int *K_, int Nfield_);

}
class SWEConventional2d :
	public SWEAbstract2d
{
public:
	SWEConventional2d();
	~SWEConventional2d();
	void EvaluatePostFunc(double *fphys);
	void UpdateWetDryState(double *fphys);
	void pre_UpdateWetDryState(double *fphys);
};

