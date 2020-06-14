#pragma once
#include"MeshUnion.h"
#include"cblas.h"
#include<math.h>

extern const MeshUnion *meshunion;
extern MeshUnion mesh;

class FrictionTermSolver
{
public:
	//FrictionTermSolver();
	FrictionTermSolver(double r_, double n_, double ks_, double den_, double TimeInterval_, double Time1_, double Time2_);
	~FrictionTermSolver();

	//void evaluateFrictionTermRHS(double gra, double hmin, double *fphys, double *frhs, double *Wave, double time);
	void evaluateFrictionTermRHS(double gra, double hmin, double * fphys, double * frhs);
	void EvaluateWaveFrictionCoefficient(double *fphys, double *Wave, double time, double hmin, double gra);
	int  EvaluateTimeStep(double t);

	double r;
	double n;//roughness coefficient
	double ks;//Equivalent roughness height
	double den; //Density of water
	double *WaveFrictionCoefficient;
	double Time1;
	double Time2;
	double TimeInterval;

	//double g;

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet = 4,		//well cell(SWE)
		Dry = 5,		//dry cell(SWE)
		PartialWet = 6,
		PartialWetFlood = 7,
		PartialWetDamBreak = 8
	} enumsweregion;
};