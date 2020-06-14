#ifdef _OPENMP
#include <omp.h>
#endif
//
//#define DEBUG
//
//#include "mex.h"
//
//#define NRHS 5
//#define NLHS 2

void c_boundary_EvaluateSurfValue(const double *FToM_, const double *FToE_, const double *FToN1_, const double *FToN2_, double *fphys_, double *fm_, double *fp_, const int *Nfp_, const  int *Ne_, const  int *Np_, const  int *K_, const  int Nfield_)
{
	//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//	/* check input & output */
	//	if (nrhs != NRHS) {
	//		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
	//		mexPrintf("%d inputs required.\n", NRHS);
	//	}

	//	if (nlhs != NLHS) {
	//		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
	//		mexPrintf("%d inputs required.\n", NLHS);
	//	}

	double *FToM = FToM_;
	double *FToE = FToE_;
	double *FToN1 = FToN1_;
	double *FToN2 = FToN2_;
	//const int m1 = (int)FToM[0] - 1;
	//mxArray *mxPhysM = mxGetCell(prhs[4], m1);
	//double *mxPhysM = fphys_;
	double *fphysM = fphys_;  // local phys field

	const int Nfp = *Nfp_;
	const int Ne = *Ne_;                 // num of edges
	//const mwSize *dims = mxGetDimensions(mxPhysM);  // local phys field dimension
	const int Np = *Np_;                         // num of interp nodes
	const int K = *K_;                          // num of elements
	int Nfield = Nfield_;

	//if (mxGetNumberOfDimensions(mxPhysM) > 2) {
	//	Nfield = dims[2];
	//}
	//else {
	//	Nfield = 1;
	//}

	// #ifdef DEBUG
	//   mexPrintf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
	// #endif

	//const size_t ndimOut = 3;
	//const mwSize dimOut[3] = { Nfp, Ne, Nfield };

	//plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	//plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *fM = fm_;
	double *fP = fp_;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

	for (int fld = 0; fld < Nfield; fld++) {
		double *fM_ = fM + Nfp * Ne * fld;
		double *fP_ = fP + Nfp * Ne * fld;
		double *fval = fphysM + Np * K * fld;

		for (int k = 0; k < Ne; k++) {
			const int e1 = (int)FToE[2 * k] - 1;

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + k * Nfp;
				const int n1 = (int)FToN1[sk] - 1 + e1 * Np;
				fM_[sk] = fval[n1];
			}

			//const int m2 = (int)FToM[2 * k + 1] - 1;
			const int e2 = (int)FToE[2 * k + 1] - 1;
			//mxArray *mxPhysP = mxGetCell(prhs[4], m2);
			double *mxPhysP = fphys_;
			//dims = mxGetDimensions(mxPhysP);
			const int Np2 = *Np_;
			const int K2 = *K_;

			double *fval2 = mxPhysP + Np2 * K2 * fld;  // adjacent phys field

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + k * Nfp;
				const int n1 = (int)FToN2[sk] - 1 + e2 * Np2;
				fP_[sk] = fval2[n1];
			}
		}
	}
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
//  double *FToM = mxGetPr(prhs[0]);
//  double *FToE = mxGetPr(prhs[1]);
//  double *FToN1 = mxGetPr(prhs[2]);
//  double *FToN2 = mxGetPr(prhs[3]);
//  const int m1 = (int)FToM[0] - 1;
//  mxArray *mxPhysM = mxGetCell(prhs[4], m1);
//  double *fphysM = mxGetPr(mxPhysM);  // local phys field
//
//  const int Nfp = mxGetM(prhs[2]);
//  const int Ne = mxGetN(prhs[2]);                 // num of edges
//  const mwSize *dims = mxGetDimensions(mxPhysM);  // local phys field dimension
//  const int Np = dims[0];                         // num of interp nodes
//  const int K = dims[1];                          // num of elements
//  int Nfield;
//
//  if (mxGetNumberOfDimensions(mxPhysM) > 2) {
//    Nfield = dims[2];
//  } else {
//    Nfield = 1;
//  }
//
//  // #ifdef DEBUG
//  //   mexPrintf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
//  // #endif
//
//  const size_t ndimOut = 3;
//  const mwSize dimOut[3] = {Nfp, Ne, Nfield};
//
//  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//  double *fM = mxGetPr(plhs[0]);
//  double *fP = mxGetPr(plhs[1]);
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//
//  for (int fld = 0; fld < Nfield; fld++) {
//    double *fM_ = fM + Nfp * Ne * fld;
//    double *fP_ = fP + Nfp * Ne * fld;
//    double *fval = fphysM + Np * K * fld;
//
//    for (int k = 0; k < Ne; k++) {
//      const int e1 = (int)FToE[2 * k] - 1;
//
//      for (int n = 0; n < Nfp; n++) {
//        const int sk = n + k * Nfp;
//        const int n1 = (int)FToN1[sk] - 1 + e1 * Np;
//        fM_[sk] = fval[n1];
//      }
//
//      const int m2 = (int)FToM[2 * k + 1] - 1;
//      const int e2 = (int)FToE[2 * k + 1] - 1;
//      mxArray *mxPhysP = mxGetCell(prhs[4], m2);
//      dims = mxGetDimensions(mxPhysP);
//      const int Np2 = dims[0];
//      const int K2 = dims[1];
//
//      double *fval2 = mxGetPr(mxPhysP) + Np2 * K2 * fld;  // adjacent phys field
//
//      for (int n = 0; n < Nfp; n++) {
//        const int sk = n + k * Nfp;
//        const int n1 = (int)FToN2[sk] - 1 + e2 * Np2;
//        fP_[sk] = fval2[n1];
//      }
//    }
//  }
//}