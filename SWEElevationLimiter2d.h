#pragma once

#include<iostream>
#include"MeshUnion.h"
#include"cblas.h"

extern "C" {

	void c_EvaluateVertAverage(double *cvar_, int *Nv_, double *Nvc_, double *VToM_, double *VToK_, double *VToW_, double *fvert_, double *fvmin_, double *fvmax_, int Nvcmax_);
	void c_VertLimit2d(double *fphys_, double *x_, double *y_, double *xc_, double *yc_, double *vx_, double *vy_, double *fvert_, double *fvmin_, double *fvmax_, double  *cvar_, double *EToV_, double *Fmask_, int *Np_, int *Nv_, int *K_, int *Nfp_);
}

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class SWEElevationLimiter2d
{
public:
	SWEElevationLimiter2d();
	~SWEElevationLimiter2d();

	void apply(double *fphys);
	void matLimit(double *fphys, int fieldId);
	void EvaluateVertAverage(double *fphys, int fieldId, double *fvert, double *fvmin, double *fvmax, double *cvar);
protected:
	int Nvcmax;
	double *Nvc, *VToK, *VToM, *VToW;

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet = 4,		//well cell(SWE)
		Dry = 5,		//dry cell(SWE)
		PartialWet = 6,
		PartialWetFlood = 7,
		PartialWetDamBreak = 8
	} enumsweregion;
};


