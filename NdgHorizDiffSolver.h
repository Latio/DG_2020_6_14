#pragma once
#include "AbstractDiffSolver.h"
class NdgHorizDiffSolver :
	public AbstractDiffSolver
{
public:
	NdgHorizDiffSolver();
	~NdgHorizDiffSolver();


	void UpdatePenaltyParameter(double *h);
	void CalculateAuxialaryVariable(double *fphys, double *Kappa, int VarIndex, double *InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp);
	void evaluate_pxy(double *r, double *s, double *fphys, double *pxy);
	void evaluate_inner_pxy(double *n, double *KappaM, double *KappaP, double *fm, double *fp, double *pxy_n, int Ne, int Nfp);
	void evaluate_bound_pxy(double *n, double *KappaM, double *fm, double *fp, double *pxy_n, int Ne, int Nfp);
	//void CalculatePartDerivTermX(double *frhs_tmp, double * pfield, double * Kappa, double * fphys, int Prantl, double * InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp);
	//void CalculatePartDerivTermY(double *frhs_tmp, double * pfield, double * Kappa, double * fphys, int Prantl, double * InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp);

	void CalculatePartDerivTermXY(double *frhs_tmp, double * pfield, double * Kappa, double * fphys, int Prantl, double * InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp, double *r, double *s, double *n_inner, double *n_bound);

	void assembleMassMatrix();


	double *InnerEdgeTau;
	double *BoundaryEdgeTau;
	double *Kappa;
	double *px;
	double *py;
	double *M;
	double *invM;

};

