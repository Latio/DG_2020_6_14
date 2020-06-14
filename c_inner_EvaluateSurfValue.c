#include "cblas.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void c_inner_EvaluateSurfValue(double *FToE_, double *FToN1_, double *FToN2_, double *fphys_, double *fm_, double *fp_, int *Nfp_, int *Ne_, int *Np_, int *K_, int Nfield_)
{
	double *FToE = FToE_;
	double *FToN1 = FToN1_;
	double *FToN2 = FToN2_;
	double *fphys = fphys_;
	const int Nfp = *Nfp_;
	const int Ne = *Ne_;
	const int Np = *Np_;
	const int K = *K_;
	const int Nfield = Nfield_;

#ifdef DEBUG
	printf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
#endif

	double *fM = fm_;
	double *fP = fp_;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

	for (int fld = 0; fld < Nfield; fld++) {
		double *fM_ = fM + Nfp * Ne * fld;
		double *fP_ = fP + Nfp * Ne * fld;
		double *fval = fphys + Np * K * fld;

		for (int k = 0; k < Ne; k++) {
			const int e1 = (int)FToE[2 * k] - 1;
			const int e2 = (int)FToE[2 * k + 1] - 1;

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + k * Nfp;

				const int n1 = (int)FToN1[sk] + e1 * Np - 1;
				const int n2 = (int)FToN2[sk] + e2 * Np - 1;

				fM_[sk] = fval[n1];
				fP_[sk] = fval[n2];
			}
		}
	}
}

