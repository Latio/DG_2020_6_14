#include "AbstractDiffSolver.h"


AbstractDiffSolver::AbstractDiffSolver() {
	requestmemory(&nv, meshunion->cell_p->Np, meshunion->K);

};


AbstractDiffSolver::~AbstractDiffSolver() {

	freememory(&nv);
};


void AbstractDiffSolver::Evaluate_rdhuv(double *rs, double *d, double *rdhuv, double *temp)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	multiply(d, rdhuv, temp, Np, K, Np);
	dotmul(num, rs, temp, temp);
}

void AbstractDiffSolver::EvaluateSurfValue(double *fm_, double  *fp_, double *FToE_, double *FToN1_, double *FToN2_, double *Kappa_, int Nfp_, int Ne_, int Np_, int K_) {
	c_EvaluateSurfValue(fm_, fp_, FToE_, FToN1_, FToN2_, Kappa_, Nfp_, Ne_, Np_, K_);
};