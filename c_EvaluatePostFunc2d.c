#ifdef _OPENMP
#include <omp.h>
#endif

#define NVAR 3

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

void c_EvaluatePostFunc2d(double hmin_, double *fphys_, double *hc_, double *qxc_, double *qyc_, int *K_, int *Np_)
{

	/* get inputs */
	double hcrit = hmin_;
	double* fphys = fphys_;
	double* hc = hc_;
	double* qxc = qxc_;
	double* qyc = qyc_;

	/* get dimensions */
	const int Np = *Np_;
	const int K = *K_;

	double* h = fphys;
	double* qx = fphys + K * Np;
	double* qy = fphys + 2 * K * Np;

	double* h_pos = fphys;
	double* qx_pos = h_pos + K * Np;
	double* qy_pos = h_pos + 2 * K * Np;

	const double ksi = 0.0;
	// cell area and scalar averages
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		double hmean = hc[k];
		double qxmean = qxc[k];
		double qymean = qyc[k];

		if (hmean <= ksi) {
			for (int n = 0; n < Np; n++) {
				int sk = k * Np + n;
				h[sk] = 0;
				qx[sk] = 0;
				qy[sk] = 0;
			}
			continue;
		}

		double hmin = h[k * Np];
		for (int n = 0; n < Np; n++) {
			hmin = min(hmin, h[k * Np + n]);
		}

		double theta;
		if (hmin < hmean) {
			theta = min((hmean - ksi) / (hmean - hmin), 1.0);
		}
		else {
			theta = 0.0;
		}

		for (int n = 0; n < Np; n++) {
			int sk = k * Np + n;
			h_pos[sk] = theta * (h[sk] - hmean) + hmean;
			qx_pos[sk] = theta * (qx[sk] - qxmean) + qxmean;
			qy_pos[sk] = theta * (qy[sk] - qymean) + qymean;

			if (h_pos[sk] < hcrit) {  // dry nodes
				qx_pos[sk] = 0.0;
				qy_pos[sk] = 0.0;
			}
		}
	}
	return;
}
