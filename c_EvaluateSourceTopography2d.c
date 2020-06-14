//#include "mex.h"
#include "SWE2d.h"

//#define NRHS 4
//#define NLHS 1
#define NVAR 3

void c_EvaluateSourceTopography2d(double gra_, signed char *status_, double *fphys_, double *zGrad_, double *frhs_temp, int *Np_, int *K_, int Nfield_)
{

	double gra = gra_;
	signed char* regType = status_;
	// double* fphys = mxGetPr(prhs[2]);
	double* zgrad = zGrad_;

	PhysField fphys = convertMexToPhysFieldp(fphys_, Np_, K_, Nfield_);
	const int Np = fphys.Np;
	const int K = fphys.K;
	const int Ntmp = Np * K;

	double* bx = zgrad;
	double* by = zgrad + Ntmp;

	PhysField source = convertMexToPhysFieldp(frhs_temp, Np_, K_, NVAR);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		NdgRegionType type = (NdgRegionType)regType[k];
		if (type == NdgRegionDry) {
			continue;
		}
		if (type == NdgRegionPartialWetFlood) {
			continue;
		}

		for (int n = 0; n < Np; n++) {
			int sk = k * Np + n;
			const double eta_ = fphys.h[sk] + fphys.z[sk];
			source.hu[sk] = -gra * eta_ * bx[sk];
			source.hv[sk] = -gra * eta_ * by[sk];
		}
	}

	return;
}
