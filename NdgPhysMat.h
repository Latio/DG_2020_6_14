#pragma once
#include"NdgQuadFreeStrongFormAdvSolver2d.h"
//#include"NdgWaveCurrentVisSolver2d.h"
#include"NdgSWEHorizSmagrinskyDiffSolver.h"
#include"AbstractOutputFile.h"
//#include<vector>
#include<fstream>
//#include<iomanip>
//#include<algorithm>
#include<time.h>

extern "C" {

}
class NdgPhysMat :public SWEPreBlanaced2d
{
public:
	NdgPhysMat();
	~NdgPhysMat();

	void matSolver();
	void matEvaluateSSPRK22();

	void UpdateExternalField(double tloc, double *fphys);

	void EvaluateRHS(double *fphys, double *frhs, double time);
	void UpdateOutputResult(double& time, double *fphys, int Nvar);
	void EvaluateLimiter(double *fphys);
	//void evaluateViscosityRHS(double *fphys, double time);

	//int Calculate_Volume(double *huvxy, double *r, double *s, double *dr, double *fphys_n, double *ds, double *lift, double *n, double *js, double *fm_n, double *varflux_n, double *j);

	//double* EvaluatePostFunc(double *fphys);
	//double UpdateTimeInterval(double *fphys);
	//void EvaluateSourceTerm(double *fphys);

protected:

	SWEAbstract2d sweabstract2d;
	SWEConventional2d sweconventional2d;
	NdgQuadFreeStrongFormAdvSolver2d ndgquadfreestrongformadvsolver2d;
	AbstractOutputFile abstractoutputfile;

	NdgSWEHorizSmagrinskyDiffSolver ndgswehorizsmagrinskydiffsolver;
	//NdgWaveCurrentVisSolver2d ndgwavecurrentvissolver2d;

	//double *fext;
	//int *outputfile;
	//double *limiter;
	//double *davectionSolver;
	//double *viscositySolver;
	//double *NonhydrostaticSolver;
	//double ftime;
	std::string casename;

	double *frhs;
	double *fext;
	//double *fphys0;
	double *fphys;
	double *zGrad;
	std::vector<double> tidal;
	std::vector<int> obeindex;

	double tidalinterval;

	//double ftime;
	int outputIntervalNum;
	int *Np;
	int *K;
	//int Nfield;
	int *Nv;
	//int Nvar;
	int *boundarydge_Nfp;
	int *boundarydge_Ne;

	double startTime, finalTime;



	double *InnerEdgefm2d, *InnerEdgefp2d, *BoundaryEdgefm2d, *BoundaryEdgefp2d;
	//double gra;
	//double hmin;
	//int Np;// dimension 
	//int K;// dimension 
	//int Nfield;// dimension 
	//int Nvar;// dimension 

};

