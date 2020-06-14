#pragma once
#include"MeshUnion.h"
#include"cblas.h"
#include<math.h>

extern const MeshUnion *meshunion;
extern MeshUnion mesh;

class RollerWaveRadiationSolver
{


public:
	double den;
	double TimeInterval;
	double Time1;
	double Time2;
	double *WaveField;

	RollerWaveRadiationSolver();
	~RollerWaveRadiationSolver();

	void valuateWaveRadiationRHS(double time, double *frhs, double *fphys, double *Wave, double hmin, double gra);
	int EvaluateTimeStep(double t);
	void Evaluate_SR(double *r, double *d, double *SR, double *temp);

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet = 4,		//well cell(SWE)
		Dry = 5,		//dry cell(SWE)
		PartialWet = 6,
		PartialWetFlood = 7,
		PartialWetDamBreak = 8
	} enumsweregion;



};

