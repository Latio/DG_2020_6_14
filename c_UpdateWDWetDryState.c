#include "SWE2d.h"
#include<stdbool.h>

NdgRegionType getCellWDType(const int Np, const double hmin, double *h,
	double *z) {
	NdgRegionType type = NdgRegionWet; // initialize to wet element
	int wetNodeNum = 0;
	bool partialWetFlood = true;
	for (int n = 0; n < Np; n++) {
		if (h[n] < hmin) { // dry nodes
			if ((h[n] + z[n] - hmin) > z[n]) {
				partialWetFlood = false;
			}
		}
		else if (h[n] > hmin) { // wet nodes
			wetNodeNum += 1;
		}
	}

	if (wetNodeNum == 0) { // dry element
		type = NdgRegionDry;
	}
	else if (wetNodeNum < Np) { // partial wet element
		if (partialWetFlood) {
			type = NdgRegionPartialWetFlood;
		}
		else {
			type = NdgRegionPartialWetDamBreak;
		}
	}
	return type;
}



void c_UpdateWDWetDryState(double hmin_, double *fphys_, signed char *status_, int *Np_, int *K_, int Nfield_)
{
	double hmin = hmin_;
	PhysField fphys = convertMexToPhysFieldp(fphys_, Np_, K_, Nfield_);

	//const int ndimOut = 1;
	//const mwSize dimLen[1] = { fphys.K };
	//plhs[0] = mxCreateNumericArray(ndimOut, dimLen, mxINT8_CLASS, mxREAL);
	signed char *cellType = status_;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

	for (int k = 0; k < fphys.K; k++) {
		cellType[k] = (signed char)getCellWDType(
			fphys.Np, hmin, fphys.h + k * fphys.Np, fphys.z + k * fphys.Np);
	}
	return;
}
