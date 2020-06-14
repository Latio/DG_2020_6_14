#include <math.h>
//#include "mex.h"
#include "SWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/** */
void evaluateNodalFlux(double hcrit, ///< water depth threshold
	double gra,   ///< gravity accerelation
	double h,     ///< water depth
	double qx,    ///< flux
	double qy,    ///< flux
	double z,     ///< bottom elevation
	double *Eh,   ///< flux term
	double *Eqx,  ///< flux term
	double *Eqy,  ///< flux term
	double *Gh,   ///< flux term
	double *Gqx,  ///< flux term
	double *Gqy   ///< flux term
) {
	if (h > hcrit) {
		double h2 = 0.5 * gra * (h * h - z * z);
		double huv = qx * qy / h;
		*Eh = qx;
		*Gh = qy;
		*Eqx = (qx * qx / h + h2);
		*Gqx = huv;
		*Eqy = huv;
		*Gqy = (qy * qy / h + h2);
	}
	else { // for dry nodes
		*Eh = 0;
		*Eqx = 0;
		*Eqy = 0;
		*Gh = 0;
		*Gqx = 0;
		*Gqy = 0;
	}
	return;
}

void c_EvaluateFlux2d(double hmin_, double gra_, signed char *status_, double *fphys_, int *Np_, int *K_, double *E_, double *G_)
{
	double hcrit = hmin_;
	double gra = gra_;
	signed char *regType = status_;
	double *fphys = fphys_;

	//const mwSize *dims = mxGetDimensions(prhs[3]);

	const int Np = *Np_;
	const int K = *K_;

	//const int ndimOut = 3;
	//const mwSize dimOut[3] = { Np, K, NVAR };
	//plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	//plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

	double *E = E_;
	double *G = G_;

	double *Eh = E;
	double *Ehu = E + K * Np;
	double *Ehv = E + 2 * K * Np;
	double *Gh = G;
	double *Ghu = G + K * Np;
	double *Ghv = G + 2 * K * Np;

	double *h = fphys;
	double *hu = fphys + K * Np;
	double *hv = fphys + 2 * K * Np;
	double *z = fphys + 3 * K * Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		NdgRegionType type = (NdgRegionType)regType[k];
		if (type == NdgRegionDry) {
			continue;
		}

		for (int n = 0; n < Np; n++) {
			int sk = k * Np + n;
			if (type == NdgRegionPartialWetFlood) {
				evaluateNodalFlux(hcrit, 0.0, h[sk], hu[sk], hv[sk], z[sk], Eh + sk,
					Ehu + sk, Ehv + sk, Gh + sk, Ghu + sk, Ghv + sk);
			}
			else {
				evaluateNodalFlux(hcrit, gra, h[sk], hu[sk], hv[sk], z[sk], Eh + sk,
					Ehu + sk, Ehv + sk, Gh + sk, Ghu + sk, Ghv + sk);
			}
		}
	}

	return;
}




///** */
//void evaluateNodalFlux(double hcrit, ///< water depth threshold
//	double gra,   ///< gravity accerelation
//	double h,     ///< water depth
//	double qx,    ///< flux
//	double qy,    ///< flux
//	double z,     ///< bottom elevation
//	double *Eh,   ///< flux term
//	double *Eqx,  ///< flux term
//	double *Eqy,  ///< flux term
//	double *Gh,   ///< flux term
//	double *Gqx,  ///< flux term
//	double *Gqy   ///< flux term
//) {
//	if (h > hcrit) {
//		double h2 = 0.5 * gra * (h * h - z * z);
//		double huv = qx * qy / h;
//		*Eh = qx;
//		*Gh = qy;
//		*Eqx = (qx * qx / h + h2);
//		*Gqx = huv;
//		*Eqy = huv;
//		*Gqy = (qy * qy / h + h2);
//	}
//	else { // for dry nodes
//		*Eh = 0;
//		*Eqx = 0;
//		*Eqy = 0;
//		*Gh = 0;
//		*Gqx = 0;
//		*Gqy = 0;
//	}
//	return;
//}
//
//
//void c_EvaluateFlux2d(double hmin_, double gra_, signed char *status_, double *fphys_, int *Np_, int *K_, double *E_, double *G_)
//{
//
//	double hcrit = hmin_;
//	double gra = gra_;
//	signed char *regType = (signed char *)status_;
//	double *fphys = fphys_;
//
//
//	const size_t Np = *Np_;
//	const size_t K = *K_;
//
//	double *E = E_;
//	double *G = G_;
//
//	double *Eh = E;
//	double *Ehu = E + K * Np;
//	double *Ehv = E + 2 * K * Np;
//	double *Gh = G;
//	double *Ghu = G + K * Np;
//	double *Ghv = G + 2 * K * Np;
//
//	double *h = fphys;
//	double *hu = fphys + K * Np;
//	double *hv = fphys + 2 * K * Np;
//	double *z = fphys + 3 * K * Np;
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int k = 0; k < K; k++) {
//		NdgRegionType type = (NdgRegionType)regType[k];
//		if (type == NdgRegionDry) {
//			continue;
//	}
//
//		for (int n = 0; n < Np; n++) {
//			size_t sk = k * Np + n;
//			if (type == NdgRegionPartialWetFlood) {
//				evaluateNodalFlux(hcrit, 0.0, h[sk], hu[sk], hv[sk], z[sk], Eh + sk,
//					Ehu + sk, Ehv + sk, Gh + sk, Ghu + sk, Ghv + sk);
//			}
//			else {
//				evaluateNodalFlux(hcrit, gra, h[sk], hu[sk], hv[sk], z[sk], Eh + sk,
//					Ehu + sk, Ehv + sk, Gh + sk, Ghu + sk, Ghv + sk);
//			}
//		}
//}
//
//	return;
//}
