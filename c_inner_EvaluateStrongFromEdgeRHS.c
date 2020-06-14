#ifdef _OPENMP
#include <omp.h>
#endif

//#include "mex.h"
#include "cblas.h"
#include<stdlib.h>

// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

void myfree(double **arr) {
	if (*arr != NULL)
	{
		free(*arr);
		*arr = NULL;
	};
}

void c_inner_EvaluateStrongFromEdgeRHS(double *invM_, double *M_, double *FToE_, double *FToN1_, double *FToN2_, double *Js_, double *J_, double *fluxM__, double *fluxP__, double *fluxS__, double *frhs_, const int *Np_, int *K_, int *Nfp_, int *Ne_, int Nfield_)
{
	const double *invM = invM_;
	double *Mb = M_;
	double *FToE = FToE_;
	double *FToN1 = FToN1_;
	double *FToN2 = FToN2_;
	double *Js = Js_;
	double *J = J_;
	double *fluxM = fluxM__;
	double *fluxP = fluxP__;
	double *fluxS = fluxS__;

	const int Np = *Np_;  // num of interp nodes
	const int K = *K_;   // num of elements

	const int Nfp = *Nfp_;
	const int Ne = *Ne_;  // num of edges
	int Nfield = Nfield_;

	//if (Nfield_ > 2) {
	//	Nfield = Nfield_;
	//}
	//else {
	//	Nfield = 1;  // fluxM is a 2D matrix
	//}

	//const size_t ndimOut = 3;
	//const mwSize dimOut[3] = { Np, K, Nfield };

	//plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *frhs = frhs_;

	//char *chn = "N";
	//double one = 1.0, zero = 0.0;
	int oneI = 1;
	//int np = Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int fld = 0; fld < Nfield; fld++) {
		double *rhs = frhs + Np * K * fld;
		double *fluxM_ = fluxM + Nfp * Ne * fld;
		double *fluxP_ = fluxP + Nfp * Ne * fld;
		double *fluxS_ = fluxS + Nfp * Ne * fld;

		for (int k = 0; k < Ne; k++) {  // evaluate rhs on each edge
			const int e1 = (int)FToE[2 * k] - 1;
			const int e2 = (int)FToE[2 * k + 1] - 1;
			const int ind1 = e1 * Np - 1;
			const int ind2 = e2 * Np - 1;

			const int ind = k * Nfp;

			double *rhsM = (double*)malloc(sizeof(double)*Nfp);
			double *rhsP = (double*)malloc(sizeof(double)*Nfp);

			for (int n = 0; n < Nfp; n++) {
				rhsM[n] = 0;
				rhsP[n] = 0;
			}

			for (int n = 0; n < Nfp; n++) {
				const int sk = n + ind;
				double dfM = fluxM_[sk] - fluxS_[sk];
				double dfP = fluxP_[sk] - fluxS_[sk];
				double j = Js[sk];

				double *mb = Mb + n * Nfp;

				for (int m = 0; m < Nfp; m++) {
					rhsM[m] += mb[m] * j * dfM;
					rhsP[m] -= mb[m] * j * dfP;
				}
			}
			for (int n = 0; n < Nfp; n++) {
				const int sk = n + ind;
				const int m1 = (int)FToN1[sk] + ind1;
				const int m2 = (int)FToN2[sk] + ind2;
				rhs[m1] += rhsM[n];
				rhs[m2] += rhsP[n];
			}

			myfree(&rhsM);
			myfree(&rhsP);
		}

		double *temp = (double*)malloc(sizeof(double)*Np);
		for (int k = 0; k < K; k++) {
			double *rhs_ = rhs + k * Np;
			double *j = J + k * Np;

			const enum CBLAS_ORDER Order = CblasColMajor;
			const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
			const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
			const int M = Np;//A的行数，C的行数
			const int N = oneI;//B的列数，C的列数
			const int K = Np;//A的列数，B的行数
			const double alpha = 1.0;
			const float beta = 0.0;
			const int lda = M;//A的行        
			const int ldb = K;//B的行
			const int ldc = M;//C的行   //如果列优先，分别写ABC的行

			//cblas_dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp, &np);
			cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, invM, lda, rhs_, ldb, beta, temp, ldc);

			for (int n = 0; n < Np; n++) {
				rhs_[n] = temp[n] / j[n];
			}
		}
		myfree(&temp);
	}
	return;
}


//#define DEBUG 0
//
//#define NRHS 10
//#define NLHS 1
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
//	double *invM = mxGetPr(prhs[0]);
//	double *Mb = mxGetPr(prhs[1]);
//	double *FToE = mxGetPr(prhs[2]);
//	double *FToN1 = mxGetPr(prhs[3]);
//	double *FToN2 = mxGetPr(prhs[4]);
//	double *Js = mxGetPr(prhs[5]);
//	double *J = mxGetPr(prhs[6]);
//	double *fluxM = mxGetPr(prhs[7]);
//	double *fluxP = mxGetPr(prhs[8]);
//	double *fluxS = mxGetPr(prhs[9]);
//
//	// dims = mxGetDimensions(prhs[6]);
//	const int Np = mxGetM(prhs[6]);  // num of interp nodes
//	const int K = mxGetN(prhs[6]);   // num of elements
//	const mwSize *dims = mxGetDimensions(prhs[7]);
//	const int Nfp = dims[0];
//	const int Ne = dims[1];  // num of edges
//	int Nfield;
//
//	if (mxGetNumberOfDimensions(prhs[7]) > 2) {
//		Nfield = dims[2];
//	}
//	else {
//		Nfield = 1;  // fluxM is a 2D matrix
//	}
//
//	const size_t ndimOut = 3;
//	const mwSize dimOut[3] = { Np, K, Nfield };
//
//	plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
//	double *frhs = mxGetPr(plhs[0]);
//
//	char *chn = "N";
//	double one = 1.0, zero = 0.0;
//	ptrdiff_t oneI = 1;
//	ptrdiff_t np = Np;
//
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int fld = 0; fld < Nfield; fld++) {
//		double *rhs = frhs + Np * K * fld;
//		double *fluxM_ = fluxM + Nfp * Ne * fld;
//		double *fluxP_ = fluxP + Nfp * Ne * fld;
//		double *fluxS_ = fluxS + Nfp * Ne * fld;
//
//		for (int k = 0; k < Ne; k++) {  // evaluate rhs on each edge
//			const int e1 = (int)FToE[2 * k] - 1;
//			const int e2 = (int)FToE[2 * k + 1] - 1;
//			const int ind1 = e1 * Np - 1;
//			const int ind2 = e2 * Np - 1;
//
//			const int ind = k * Nfp;
//			double rhsM[Nfp], rhsP[Nfp];
//			for (int n = 0; n < Nfp; n++) {
//				rhsM[n] = 0;
//				rhsP[n] = 0;
//			}
//
//			for (int n = 0; n < Nfp; n++) {
//				const int sk = n + ind;
//				double dfM = fluxM_[sk] - fluxS_[sk];
//				double dfP = fluxP_[sk] - fluxS_[sk];
//				double j = Js[sk];
//
//				double *mb = Mb + n * Nfp;
//
//				for (int m = 0; m < Nfp; m++) {
//					rhsM[m] += mb[m] * j * dfM;
//					rhsP[m] -= mb[m] * j * dfP;
//				}
//			}
//
//			for (int n = 0; n < Nfp; n++) {
//				const int sk = n + ind;
//				const int m1 = (int)FToN1[sk] + ind1;
//				const int m2 = (int)FToN2[sk] + ind2;
//				rhs[m1] += rhsM[n];
//				rhs[m2] += rhsP[n];
//			}
//		}
//		double temp[Np];
//		for (int k = 0; k < K; k++) {
//			double *rhs_ = rhs + k * Np;
//			double *j = J + k * Np;
//
//			dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp,
//				&np);
//			//	const int M = 4;//A的行数，C的行数
//			//	const int N = 2;//B的列数，C的列数
//			//	const int K = 3;//A的列数，B的行数
//			//	const float alpha = 1;
//			//	const float beta = 0;
//			//	const int lda = M;//A的行        //如果列优先，分别写ABC的行
//			//	const int ldb = K;//B的行
//			//	const int ldc = N;//C的列
//			//	const double A[K*M] = { 1,4,7,8,2,5,8,7,3,6,9,6 };
//			//	const double B[K*N] = { 5,3,1,4,2,0 };
//			//	double C[M*N];
//			//
//			//	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
//				  // copy rhs
//			for (int n = 0; n < Np; n++) {
//				rhs_[n] = temp[n] / j[n];
//			}
//		}
//	}
//	return;
//}