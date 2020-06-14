#ifndef __mxSWE_H__
#define __mxSWE_H__
//void myfree(double**arr);
//#include "mex.h"

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define EPS 1e-6

typedef enum {
	NdgRegionNormal = 1,
	NdgRegionRefine = 2,
	NdgRegionSponge = 3,
	NdgRegionWet = 4,
	NdgRegionDry = 5,
	NdgRegionPartialWet = 6,
	NdgRegionPartialWetFlood = 7,
	NdgRegionPartialWetDamBreak = 8
} NdgRegionType;

typedef struct {
	int Np;      ///< length of 1st dimension
	int K;       ///< length of 2nd dimension
	int Nfield;  ///< length of 3rd dimension
	double *h;
	double *hu;
	double *hv;
	double *z;
} PhysField;

/** convert mex variable to PhysVolField structure */
PhysField convertMexToPhysFieldp(const double *, int *Np_, int *K_, int Nfield_);
PhysField convertMexToPhysField(const double *, int Np_, int K_, int Nfield_);

/** Evaluate the flow rate depending on the depth threshold */
void evaluateFlowRateByDeptheThreshold(
	const double,  ///< depth threshold
	const double,      ///< depth
	const double,     ///< water flux
	const double,     ///< water flux
	double *,           ///< result velocity
	double *            ///< velocity
);

//void evaluateWetDryInterface(signed char *, const double *, double *);

#endif  //__mxSWE_H__




