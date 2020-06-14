#include <math.h>
#include "SWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateSurfFluxTerm(const double hmin,  ///< water depth threshold
	const double gra,   ///< gravity acceleration
	const double h,     ///< water depth
	const double hu,    ///< water flux
	const double hv,    ///< water flux
	double* E,          ///< flux term
	double* G           ///< flux term
) {
	double u, v;
	evaluateFlowRateByDeptheThreshold(hmin, h, hu, hv, &u, &v);
	const double huv = h * u * v;
	const double h2 = h * h;
	E[0] = hu;
	G[0] = hv;
	E[1] = h * u * u + 0.5 * gra * h2;
	G[1] = huv;
	E[2] = huv;
	G[2] = h * v * v + 0.5 * gra * h2;
	return;
}


void c_EvaluateSurfFlux(double hmin_, double gra_, double *nx_, double *ny_, double *fm_, double *fluxM_, int *fm_Nfp, int *fm_Ne, int fm_Nfield)
{

	double hcrit = hmin_;
	double gra = gra_;
	double *nx = nx_;
	double *ny = ny_;
	int NVAR = 3;

	PhysField fm = convertMexToPhysFieldp(fm_, fm_Nfp, fm_Ne, fm_Nfield);

	const int TNfp = fm.Np;
	const int K = fm.K;

	PhysField surfFlux = convertMexToPhysField(fluxM_, TNfp, K, NVAR);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

	for (int k = 0; k < K; k++) {
		for (int n = 0; n < TNfp; n++) {
			const int sk = k * TNfp + n;

			const double nx__ = nx[sk];
			const double ny__ = ny[sk];

			double E[3], G[3];
			evaluateSurfFluxTerm(hcrit, gra, fm.h[sk], fm.hu[sk], fm.hv[sk], E, G);

			surfFlux.h[sk] = nx__ * E[0] + ny__ * G[0];
			surfFlux.hu[sk] = nx__ * E[1] + ny__ * G[1];
			surfFlux.hv[sk] = nx__ * E[2] + ny__ * G[2];
		}
	}
	return;
}
