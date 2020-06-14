#include "SWE2d.h"
#include <math.h>
//#include<stdint.h>
#include<stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateFlowRateByCellState(
	const NdgRegionType type,  ///< cell types
	const double h,            ///< water depth
	const double hu,           ///< flux
	const double hv,           ///< flux
	double* u,                 ///< velocity result
	double* v                  ///< velocity result
) {
	if (type == NdgRegionWet) {
		*u = hu / h;
		*v = hv / h;
	}
	else {
		*u = 0;
		*v = 0;
	}
	return;
}


double c_UpdateTimeInterval2d(double hmin_, double gra_, int N_, double *dx_, signed char *status_, double *const fphys_, int *Np_, int *K_, int Nfield_)
{
	double gra = gra_;
	int N = N_;
	double *dx = dx_;
	signed char* regionType =status_;

	PhysField fphys = convertMexToPhysFieldp(fphys_, Np_, K_, Nfield_);

	const int Np = fphys.Np;
	const int K = fphys.K;

	double dt = 1e6;

	for (int k = 0; k < K; k++) {
		NdgRegionType type = (NdgRegionType)regionType[k];
		if (type == NdgRegionDry) {
			continue;
		}
		double dx__ = dx[k];
		for (int n = 0; n < Np; n++) {
			const int sk = k * Np + n;
			const double h_ = fphys.h[sk];

			double u, v;
			evaluateFlowRateByCellState(type, h_, fphys.hu[sk], fphys.hv[sk], &u, &v);
			const double spe = sqrt(u * u + v * v);
			const double dtloc = dx__ / (spe + sqrt(gra * h_)) / (2 * N + 1);
			dt = min(dt, dtloc);
		}
	}
	//for (int i = 0; i < 360; i++)
	//{
	//	printf("dx[%d]:%f\n", i, dx[i]);
	//}

	return dt;
}

