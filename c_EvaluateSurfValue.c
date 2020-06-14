#ifdef _OPENMP
#include <omp.h>
#endif

//#include "mex.h"
//#include "blas.h"

//#define NRHS 4
//#define NLHS 2

// #define DEBUG

void c_EvaluateSurfValue(double *fm_, double  *fp_, double *FToE_, double *FToN1_, double *FToN2_, double *Kappa_, int Nfp_, int Ne_, int Np_, int K_)
{

	double *FToE = FToE_;
	double *FToN1 = FToN1_;
	double *FToN2 = FToN2_;
	double *fphys = Kappa_;

	//const mwSize *dims = mxGetDimensions(prhs[1]);
	const int Nfp = Nfp_;
	const int Ne = Ne_;  // num of edges
	//dims = mxGetDimensions(prhs[3]);
	const int Np = Np_;  // num of interp nodes
	const int K = K_;   // num of elements
	int Nfield = 1;

	//if (mxGetNumberOfDimensions(prhs[3]) > 2) {
	//	Nfield = dims[2];
	//}
	//else {
	//	Nfield = 1;
	//}

#ifdef DEBUG
	mexPrintf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
#endif

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


//
//
//
//
//
//
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//	/* check input & output */
//	if (nrhs != NRHS) {
//		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
//		mexPrintf("%d inputs required.\n", NRHS);
//	}
//
//	if (nlhs != NLHS) {
//		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
//		mexPrintf("%d inputs required.\n", NLHS);
//	}
//
//	double *FToE = mxGetPr(prhs[0]);
//	double *FToN1 = mxGetPr(prhs[1]);
//	double *FToN2 = mxGetPr(prhs[2]);
//	double *fphys = mxGetPr(prhs[3]);
//
//	const mwSize *dims = mxGetDimensions(prhs[1]);
//	const int Nfp = dims[0];
//	const int Ne = dims[1];  // num of edges
//	dims = mxGetDimensions(prhs[3]);
//	const int Np = dims[0];  // num of interp nodes
//	const int K = dims[1];   // num of elements
//	int Nfield = 1;
//
//	//if (mxGetNumberOfDimensions(prhs[3]) > 2) {
//	//	Nfield = dims[2];
//	//}
//	//else {
//	//	Nfield = 1;
//	//}
//
//#ifdef DEBUG
//	mexPrintf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
//#endif
//
//	//const size_t ndimOut = 3;
//	//const mwSize dimOut[3] = { Nfp, Ne, Nfield };
//
//	//plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//	//plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//	double *fM = fm_;
//	double *fP = fp_;
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//
//	for (int fld = 0; fld < Nfield; fld++) {
//		double *fM_ = fM + Nfp * Ne * fld;
//		double *fP_ = fP + Nfp * Ne * fld;
//		double *fval = fphys + Np * K * fld;
//
//		for (int k = 0; k < Ne; k++) {
//			const int e1 = (int)FToE[2 * k] - 1;
//			const int e2 = (int)FToE[2 * k + 1] - 1;
//
//			for (int n = 0; n < Nfp; n++) {
//				const int sk = n + k * Nfp;
//
//				const int n1 = (int)FToN1[sk] + e1 * Np - 1;
//				const int n2 = (int)FToN2[sk] + e2 * Np - 1;
//
//				fM_[sk] = fval[n1];
//				fP_[sk] = fval[n2];
//			}
//		}
//	}
//}


//#include "mex.h"
//#include "blas.h"
//
//#define NRHS 4
//#define NLHS 2
//
//// #define DEBUG
//
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//	/* check input & output */
//	if (nrhs != NRHS) {
//		mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
//		mexPrintf("%d inputs required.\n", NRHS);
//	}
//
//	if (nlhs != NLHS) {
//		mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
//		mexPrintf("%d inputs required.\n", NLHS);
//	}
//
//	double *FToE = mxGetPr(prhs[0]);
//	double *FToN1 = mxGetPr(prhs[1]);
//	double *FToN2 = mxGetPr(prhs[2]);
//	double *fphys = mxGetPr(prhs[3]);
//
//	const mwSize *dims = mxGetDimensions(prhs[1]);
//	const int Nfp = dims[0];
//	const int Ne = dims[1];  // num of edges
//	dims = mxGetDimensions(prhs[3]);
//	const int Np = dims[0];  // num of interp nodes
//	const int K = dims[1];   // num of elements
//	int Nfield;
//
//	if (mxGetNumberOfDimensions(prhs[3]) > 2) {
//		Nfield = dims[2];
//	}
//	else {
//		Nfield = 1;
//	}
//
//#ifdef DEBUG
//	mexPrintf("Nfp = %d, Ne = %d, Np = %d, K = %d\n", Nfp, Ne, Np, K);
//#endif
//
//	const size_t ndimOut = 3;
//	const mwSize dimOut[3] = { Nfp, Ne, Nfield };
//
//	plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//	plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//	double *fM = mxGetPr(plhs[0]);
//	double *fP = mxGetPr(plhs[1]);
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//
//	for (int fld = 0; fld < Nfield; fld++) {
//		double *fM_ = fM + Nfp * Ne * fld;
//		double *fP_ = fP + Nfp * Ne * fld;
//		double *fval = fphys + Np * K * fld;
//
//		for (int k = 0; k < Ne; k++) {
//			const int e1 = (int)FToE[2 * k] - 1;
//			const int e2 = (int)FToE[2 * k + 1] - 1;
//
//			for (int n = 0; n < Nfp; n++) {
//				const int sk = n + k * Nfp;
//
//				const int n1 = (int)FToN1[sk] + e1 * Np - 1;
//				const int n2 = (int)FToN2[sk] + e2 * Np - 1;
//
//				fM_[sk] = fval[n1];
//				fP_[sk] = fval[n2];
//			}
//		}
//	}
//}
