//#include "mex.h"
#include "SWE2d.h"
#include <math.h>
#include<stdlib.h>
#include"cblas.h"
void myfree(double**arr);

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
	int Nfp;     ///< num of face nodes
	int Ne;      ///< num of edges
	int Nfield;  ///< num of field
	double *h;   ///< water depth at local node
	double *hu;  ///< flux at local node
	double *hv;  ///< flux at local node
	double *z;   ///< bottom elevation
} SurfNodeField;

//---------------------------------------------------------
/// \brief Implement the hydraostatic reconstruction
/// \details Details about the HR method can be found in Hou et. al. (2013)
void evaluateHydrostaticReconstructValue(
	const double hmin,                     ///< water depth threshold
	const int ind,                         ///< index
	SurfNodeField fM, SurfNodeField fP,    ///< input
	SurfNodeField fM1, SurfNodeField fP1)  ///< result
//---------------------------------------------------------
{
	double zstar = max(fM.z[ind], fP.z[ind]);
	double um, vm, up, vp;
	evaluateFlowRateByDeptheThreshold(hmin, fM.h[ind], fM.hu[ind], fM.hv[ind],
		&um, &vm);
	evaluateFlowRateByDeptheThreshold(hmin, fP.h[ind], fP.hu[ind], fP.hv[ind],
		&up, &vp);
	const double etaM = fM.h[ind] + fM.z[ind];
	const double etaP = fP.h[ind] + fP.z[ind];
	zstar = min(etaM, zstar);  // z* = min( \eta^-, z* )
	fM1.h[ind] = etaM - zstar;
	fP1.h[ind] = max(0, etaP - zstar) - max(0, fP.z[ind] - zstar);
	fM1.hu[ind] = fM.h[ind] * um;
	fP1.hu[ind] = fP.h[ind] * up;
	fM1.hv[ind] = fM.h[ind] * vm;
	fP1.hv[ind] = fP.h[ind] * vp;
	fM1.z[ind] = zstar;
	fP1.z[ind] = zstar;
	return;
}

//---------------------------------------------------------
SurfNodeField ConvertMexToSurfField_same(const double *mxField, const int *Nfp_, const int *Ne_, const int Nfield_)
//---------------------------------------------------------
{
	SurfNodeField field;
	//const mwSize *dims = mxGetDimensions(mxField);  // local phys field dimension
	const int Nfp = *Nfp_;                        // num of interp nodes
	const int Ne = *Ne_;                         // num of elements
	const int Nfield = Nfield_;
	//if (mxGetNumberOfDimensions(mxField) > 2) {
	//	Nfield = dims[2];
	//}
	//else {
	//	Nfield = 1;
	//}
	field.h = mxField;
	field.hu = field.h + Nfp * Ne;
	field.hv = field.hu + Nfp * Ne;
	field.z = field.hv + Nfp * Ne;
	field.Nfp = Nfp;
	field.Ne = Ne;
	field.Nfield = Nfield;

	return field;
}

//#define NRHS 3
//#define NLHS 2

void c_HydrostaticReconstruction(double hmin_, double *fm_, double *fp_, const int *Nfp_, const int *Ne_, const int Nfield_) {

	//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//  /* check input & output */
	//  if (nrhs != NRHS) {
	//    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
	//    mexPrintf("%d inputs required.\n", NRHS);
	//  }
	//
	//  if (nlhs != NLHS) {
	//    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
	//    mexPrintf("%d inputs required.\n", NLHS);
	//  }
	//
	//const int num = (*Nfp_)*(*Ne_)*Nfield_;

	//double *fM_t = (double*)malloc(sizeof(double)*num);
	//double *fP_t = (double*)malloc(sizeof(double)*num);

	//cblas_dcopy(num, fm_, 1, fM_t, 1);
	//cblas_dcopy(num, fp_, 1, fP_t, 1);
	double hmin = hmin_;

	//SurfNodeField fM_temp = ConvertMexToSurfField_same(fM_t, Nfp_, Ne_, Nfield_);
	//SurfNodeField fP_temp = ConvertMexToSurfField_same(fP_t, Nfp_, Ne_, Nfield_);
	SurfNodeField fP = ConvertMexToSurfField_same(fp_, Nfp_, Ne_, Nfield_);
	SurfNodeField fM = ConvertMexToSurfField_same(fm_, Nfp_, Ne_, Nfield_);
	const int Nfp = fM.Nfp;
	const int Ne = fM.Ne;

	//const size_t ndimOut = 3;
	//const mwSize dimOut[3] = { fM.Nfp, fM.Ne, fM.Nfield };
	//plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	//plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	//SurfNodeField fM1 = ConvertMexToSurfField(plhs[0]);
	//SurfNodeField fP1 = ConvertMexToSurfField(plhs[1]);

	for (int k = 0; k < Ne; k++) {
		for (int n = 0; n < Nfp; n++) {
			const int sk = k * Nfp + n;
			evaluateHydrostaticReconstructValue(hmin, sk, fM, fP, fM, fP);
		}
	}

	//myfree(&fM_t);
	//myfree(&fP_t);
}


//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//  /* check input & output */
//  if (nrhs != NRHS) {
//    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
//    mexPrintf("%d inputs required.\n", NRHS);
//  }
//
//  if (nlhs != NLHS) {
//    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
//    mexPrintf("%d inputs required.\n", NLHS);
//  }
//
//  double hmin = mxGetScalar(prhs[0]);
//  SurfNodeField fM = ConvertMexToSurfField(prhs[1]);
//  SurfNodeField fP = ConvertMexToSurfField(prhs[2]);
//  const int Nfp = fM.Nfp;
//  const int Ne = fM.Ne;
//
//  const size_t ndimOut = 3;
//  const mwSize dimOut[3] = {fM.Nfp, fM.Ne, fM.Nfield};
//  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//  SurfNodeField fM1 = ConvertMexToSurfField(plhs[0]);
//  SurfNodeField fP1 = ConvertMexToSurfField(plhs[1]);
//
//  for (int k = 0; k < Ne; k++) {
//    for (int n = 0; n < Nfp; n++) {
//      const int sk = k * Nfp + n;
//      evaluateHydrostaticReconstructValue(hmin, sk, fM, fP, fM1, fP1);
//    }
//  }
//}