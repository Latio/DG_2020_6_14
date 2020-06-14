#pragma once
//#include "NdgPhysMat.h"
#include "SWEFaceFluxSolver2d.h"
#include"SWEHLLNumFluxSolver2d.h"
#include"SWEPrebalanceVolumeFlux2d.h"
#include"SWETopographySourceTerm2d.h"
#include"SWEElevationLimiter2d.h"
#include"FrictionTermSolver.h"
#include"RollerWaveRadiationSolver.h"
#include"new_delete.h"
//#include"NdgWaveCurrentVisSolver2d.h"

extern "C" {
	//void surfluxSolver_evaluate(double hmin_, double gra_, double *nx_, double *ny_, MeshUnion *mesh_, InnerEdge *edge_);//(hmin, gra, nx, ny, fm, mesh, edge)
	double c_UpdateTimeInterval2d(double hmin_, double gra_, int N_, double *dx_, signed char *status_, double * const fphys_, int *Np_, int *K_, int Nfield_);
	void c_ImposeBoundaryCondition(double gra_, double *nx_, double *ny_, double *fp_, double *fext_, signed char *ftype_, int *Nfp_, int* Ne_, int Nfield_);
	void c_HydrostaticReconstruction(double hmin_, double *fm_, double *fp_, const int *Nfp_, const int *Ne_, const int Nfield_);
	//void c_EvaluateSurfaceValue(double *fp_, double *fm_, double hmin_, double gra_, double *eidM_, double *eidP_, double *nx_, double *ny_, signed char *eidtype_, double *fphys_, double *fext_, int *Np_, int *K_, double *TNfp_);
}

class SWEAbstract2d
{

public:

	const double gra;
	const double hmin;
	const double cfl;
	const int Nfield;
	const int Nvar;
	double *dx;

	void EvaluateSurfFlux(double *nx, double *ny, double *fm, double *fluxM, int *Nfp, int *Ne);
	void EvaluateSurfNumFlux(double *nx, double *ny, double *fm, double *fp, double *fluxS, int *Nfp, int *Ne);
	void ImposeBoundaryCondition(double *nx, double *ny, double *fm, double *fp, double *fext);
	void EvaluateSourceTerm(double *fphys, double *frhs, double *zGrad, double time);
	double UpdateTimeInterval(double *fphys);

	//void EvaluateSurfaceValueOld(double *fp, double *fm, double *fphys, double *fext);
	//void EvaluateDerivativeX(double *qx, double *fm, double *fp, double hmin, double *field, int n);
	//void EvaluateDerivativeY(double *qx, double *fm, double *fp, double hmin, double *field, int n);
	void multiply(double *const matrix1, double *const matrix2, double *result, int M_, int N_, int L_);
	void dotmul(int num, double *const matrix1, double *matrix2, double *result);
	void dotdiv(int num, double *const matrix1, double *matrix2, double *result);

	SWEAbstract2d();
	~SWEAbstract2d();

	SWEFaceFluxSolver2d swefacefluxsolver2d;
	SWEHLLNumFluxSolver2d swehllnumfluxsolver2d;
	SWEPrebalanceVolumeFlux2d sweprebalancevolumeflux2d;
	SWETopographySourceTerm2d swetopographysourceterm2d;
	FrictionTermSolver frictiontermsolver2d;
	RollerWaveRadiationSolver rollerwaveradiationsolver;
	SWEElevationLimiter2d sweelevationlimiter2d;

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet = 4,		//well cell(SWE)
		Dry = 5,		//dry cell(SWE)
		PartialWet = 6,
		PartialWetFlood = 7,
		PartialWetDamBreak = 8
	} enumsweregion;

	//NdgWaveCurrentVisSolver2d ndgwavecurrentvissolver2d;
};

