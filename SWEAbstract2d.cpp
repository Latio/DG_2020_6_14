#include "SWEAbstract2d.h"

//double SWEAbstract2d::gra = 9.8;
//double SWEAbstract2d::hmin = 0.05;
//double SWEAbstract2d::cfl = 1;

void SWEAbstract2d::multiply(double *const matrix1, double *const matrix2, double *result, int M_, int N_, int L_)
{
	const enum CBLAS_ORDER Order = CblasColMajor;
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	const int M = M_;//A的行数，C的行数
	const int N = N_;//B的列数，C的列数
	const int L = L_;//A的列数，B的行数
	const double alpha = 1.0;
	const double beta = 0.0;
	const int lda = M;//A的行        
	const int ldb = L;//B的行
	const int ldc = M;//C的行   //如果列优先，分别写ABC的行

	cblas_dgemm(Order, TransA, TransB, M, N, L, alpha, matrix1, lda, matrix2, ldb, beta, result, ldc);
	//std::cout << "NdgQuadFreeStrongFormAdvSolver2d.cpp" << std::endl;
}

void SWEAbstract2d::dotmul(int num, double *const matrix1, double *matrix2, double *result)
{

	for (int i = 0; i < num; i++)
	{
		result[i] = matrix1[i] * matrix2[i];
	}
	return;
}

void SWEAbstract2d::dotdiv(int num, double *const matrix1, double *matrix2, double *result)
{

	for (int i = 0; i < num; i++)
	{
		result[i] = matrix1[i] / matrix2[i];
	}
	return;
}

SWEAbstract2d::SWEAbstract2d() :gra(9.8), hmin(0.5), cfl(1), Nfield(7), Nvar(3), frictiontermsolver2d(0.017, 0.0133, 1.0e-3, 1000, 5, 0, 0)
{
	int *K = meshunion->K;
	double *LAV = meshunion->LAV;

	requestmemory(&dx, K);
	for (int i = 0; i < *K; i++)
	{
		*(dx + i) = pow(*(LAV + i), 0.5);
	}
}


SWEAbstract2d::~SWEAbstract2d()
{
	freememory(&dx);
}

void SWEAbstract2d::EvaluateSurfFlux(double *nx, double *ny, double *fm, double *fluxM, int *Nfp, int *Ne)
{
	swefacefluxsolver2d.surfluxSolver_evaluate(hmin, gra, nx, ny, fm, fluxM, Nfp, Ne);
};

void SWEAbstract2d::EvaluateSurfNumFlux(double *nx, double *ny, double *fm, double *fp, double *fluxS, int *Nfp, int *Ne)
{
	swehllnumfluxsolver2d.numfluxSolver_evaluate(hmin, gra, nx, ny, fm, fp, fluxS, Nfp, Ne);
};

//计算时间步长
double SWEAbstract2d::UpdateTimeInterval(double *fphys)
{
	double dt = 1e1;
	int N = *meshunion->cell_p->N;
	signed char *status = meshunion->status;
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;

	double dtm = c_UpdateTimeInterval2d(hmin, gra, N, dx, status, fphys, Np, K, Nfield);

	if (dtm > 0) {
		dt = (dt < dtm*cfl) ? dt : dtm * cfl;
	}

	return dt;
};

void SWEAbstract2d::ImposeBoundaryCondition(double *nx, double *ny, double *fm, double *fp, double *fext)
{
	signed char *ftype = meshunion->boundarydge_p->ftype;
	int *Nfp = meshunion->boundarydge_p->Nfp;
	int *Ne = meshunion->boundarydge_p->Ne;
	int Nfield = meshunion->Nfield;

	c_ImposeBoundaryCondition(gra, nx, ny, fp, fext, ftype, Nfp, Ne, Nfield);
	//fP(:,:,6) = fM(:,:,6);
	const int dis = (*Nfp)*(*Ne);
	double *fp_6 = fp + dis * 5;
	double *fm_6 = fm + dis * 5;
	cblas_dcopy(dis, fm_6, 1, fp_6, 1);

	c_HydrostaticReconstruction(hmin, fm, fp, Nfp, Ne, Nfield);

}

void SWEAbstract2d::EvaluateSourceTerm(double *fphys, double *frhs, double *zGrad, double time)
{
	//function matEvaluateSourceTerm(obj, fphys)
	//	% frhs = frhs + BottomTerm
	//	obj.matEvaluateTopographySourceTerm(fphys);
	//% frhs = frhs + CoriolisTerm
	//	obj.coriolisSolver.evaluateCoriolisTermRHS(obj, fphys);
	//% frhs = frhs + FrictionTerm
	//	obj.frictionSolver.evaluateFrictionTermRHS(obj, fphys);
	//% frhs = frhs + WindTerm
	//	obj.windSolver.evaluateWindTermRHS(obj, fphys);

	//obj.NonhydrostaticSolver.evaluateNonhydroRHS(obj, fphys);
	//end

	//sweprebalancevolumeflux2d.evaluate();
	/*ndgwavecurrentvissolver2d.evaluateViscosityRHS(fphys,time,hmin);*/

	//double *Wave;///注意
	//requestmemory(&Wave, meshunion->cell_p->Np, meshunion->K);
	swetopographysourceterm2d.EvaluateTopographySourceTerm(gra, fphys, zGrad, frhs);
	frictiontermsolver2d.evaluateFrictionTermRHS(gra, hmin, fphys, frhs);
	//rollerwaveradiationsolver.valuateWaveRadiationRHS(time, frhs, fphys, Wave, hmin, gra);

	//freememory(&Wave);
};

//void SWEAbstract2d::EvaluateSurfaceValueOld(double *fp, double *fm, double *fphys, double *fext) {
//
//	c_EvaluateSurfaceValue(fp, fm, hmin, gra, double *eidM_, double *eidP_, double *nx_, double *ny_, signed char *eidtype_, fphys, fext, int *Np_, int *K_, double *TNfp_);
//};

//void SWEAbstract2d::EvaluateDerivativeX(double *qx, double *fm, double *fp, double hmin, double *field, int n)
//{
//	int *const TNfp = meshunion->cell_p->TNfp;
//	int *const K = meshunion->K;
//	int *const Np = meshunion->cell_p->Np;
//	signed char *status = meshunion->status;
//
//	double *rx = meshunion->rx;
//	double *sx = meshunion->sx;
//	double *Dr = meshunion->cell_p->Dr;
//	double *Ds = meshunion->cell_p->Ds;
//	double *Js = meshunion->Js;  //
//	double *J = meshunion->J;
//	double *nx = meshunion->nx;
//
//	double *LIFT;
//
//
//
//	double *var_m, *var_p, *varflux, *deltaflux;
//
//	requestmemory(&var_m, TNfp, K);
//	requestmemory(&var_p, TNfp, K);
//
//	const int num = (*TNfp)*(*K);
//
//	double *fm_n, *fp_n;
//	fm_n = fm + (n - 1)*num;
//	fp_n = fp + (n - 1)*num;
//
//	for (int i = 0; i < *K; i++)
//	{
//		if (status[i] == (signed char)enumSWERegion::Wet)
//		{
//			for (int j = 0; j < *TNfp; j++)
//			{
//				var_m[i*(*TNfp) + j] = fm_n[i*(*TNfp) + j] / fm[i*(*TNfp) + j];
//				var_p[i*(*TNfp) + j] = fp_n[i*(*TNfp) + j] / fp[i*(*TNfp) + j];
//			}
//		}
//	}
//
//	cblas_dscal(num, 0.5, var_p, 1);
//	cblas_daxpy(num, -0.5, var_m, 1, var_p, 1);
//
//	double *rx_dr_field, *sx_ds_field, *nx_js_temp1, *nx_js_temp2;
//
//	requestmemory(&rx_dr_field, Np, K);
//	requestmemory(&sx_ds_field, Np, K);
//	requestmemory(&nx_js_temp1, TNfp, K);
//	requestmemory(&nx_js_temp2, Np, K);
//
//	const int num1 = (*Np)*(*K);
//
//	multiply(Dr, field, rx_dr_field, *Np, *K, *Np);
//	dotmul(num1, rx, rx_dr_field, rx_dr_field);
//
//	multiply(Ds, field, sx_ds_field, *Np, *K, *Np);
//	dotmul(num1, sx, sx_ds_field, sx_ds_field);
//
//	dotmul(num, nx, Js, nx_js_temp1);
//	dotmul(num, nx_js_temp1, var_p, nx_js_temp1);
//	dotdiv(num, nx_js_temp1, J, nx_js_temp1);
//	multiply(LIFT, nx_js_temp1, nx_js_temp2, *Np, *K, *TNfp);
//
//	dotmul(num1, rx_dr_field, sx_ds_field, qx);
//	dotmul(num1, nx_js_temp2, qx, qx);
//
//	freememory(&rx_dr_field);
//	freememory(&sx_ds_field);
//	freememory(&nx_js_temp1);
//	freememory(&nx_js_temp2);
//
//	freememory(&var_m);
//	freememory(&var_p);
//
//
//};
//
//
//void SWEAbstract2d::EvaluateDerivativeY(double *qx, double *fm, double *fp, double hmin, double *field, int n)
//{
//	int *const TNfp = meshunion->cell_p->TNfp;
//	int *const K = meshunion->K;
//	int *const Np = meshunion->cell_p->Np;
//	signed char *status = meshunion->status;
//
//	double *ry = meshunion->ry;
//	double *sy = meshunion->sy;
//	double *Dr = meshunion->cell_p->Dr;
//	double *Ds = meshunion->cell_p->Ds;
//	double *Js = meshunion->Js;  //
//	double *J = meshunion->J;
//	double *ny = meshunion->ny;
//
//	double *LIFT;
//
//
//
//	double *var_m, *var_p, *varflux, *deltaflux;
//
//	requestmemory(&var_m, TNfp, K);
//	requestmemory(&var_p, TNfp, K);
//
//	const int num = (*TNfp)*(*K);
//
//	double *fm_n, *fp_n;
//	fm_n = fm + (n - 1)*num;
//	fp_n = fp + (n - 1)*num;
//
//	for (int i = 0; i < *K; i++)
//	{
//		if (status[i] == (signed char)enumSWERegion::Wet)
//		{
//			for (int j = 0; j < *TNfp; j++)
//			{
//				var_m[i*(*TNfp) + j] = fm_n[i*(*TNfp) + j] / fm[i*(*TNfp) + j];
//				var_p[i*(*TNfp) + j] = fp_n[i*(*TNfp) + j] / fp[i*(*TNfp) + j];
//			}
//		}
//	}
//
//	cblas_dscal(num, 0.5, var_p, 1);
//	cblas_daxpy(num, -0.5, var_m, 1, var_p, 1);
//
//	double *ry_dr_field, *sy_ds_field, *ny_js_temp1, *ny_js_temp2;
//
//	requestmemory(&ry_dr_field, Np, K);
//	requestmemory(&sy_ds_field, Np, K);
//	requestmemory(&ny_js_temp1, TNfp, K);
//	requestmemory(&ny_js_temp2, Np, K);
//
//	const int num1 = (*Np)*(*K);
//
//	multiply(Dr, field, ry_dr_field, *Np, *K, *Np);
//	dotmul(num1, ry, ry_dr_field, ry_dr_field);
//
//	multiply(Ds, field, sy_ds_field, *Np, *K, *Np);
//	dotmul(num1, sy, sy_ds_field, sy_ds_field);
//
//	dotmul(num, ny, Js, ny_js_temp1);
//	dotmul(num, ny_js_temp1, var_p, ny_js_temp1);
//	dotdiv(num, ny_js_temp1, J, ny_js_temp1);
//	multiply(LIFT, ny_js_temp1, ny_js_temp2, *Np, *K, *TNfp);
//
//	dotmul(num1, ry_dr_field, sy_ds_field, qx);
//	dotmul(num1, ny_js_temp2, qx, qx);
//
//	freememory(&ry_dr_field);
//	freememory(&sy_ds_field);
//	freememory(&ny_js_temp1);
//	freememory(&ny_js_temp2);
//
//	freememory(&var_m);
//	freememory(&var_p);
//
//
//};
////void dotmul(int num, double *const matrix1, double *matrix2)
////{
////
////	for (int i = 0; i < num; i++)
////	{
////		matrix2[i] = matrix1[i] * matrix2[i];
////	}
////	return;
////}
//
////void multiply(double *const matrix1, double *const matrix2, double *result, int M_, int N_, int L_)
////{
////	const enum CBLAS_ORDER Order = CblasColMajor;
////	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
////	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
////	const int M = M_;//A的行数，C的行数
////	const int N = N_;//B的列数，C的列数
////	const int L = L_;//A的列数，B的行数
////	const double alpha = 1.0;
////	const double beta = 0.0;
////	const int lda = M;//A的行        
////	const int ldb = L;//B的行
////	const int ldc = M;//C的行   //如果列优先，分别写ABC的行
////
////	cblas_dgemm(Order, TransA, TransB, M, N, L, alpha, matrix1, lda, matrix2, ldb, beta, result, ldc);
////	//std::cout << "NdgQuadFreeStrongFormAdvSolver2d.cpp" << std::endl;
////}
////
////void dotmul(int num, double *const matrix1, double *matrix2)
////{
////
////	for (int i = 0; i < num; i++)
////	{
////		matrix2[i] = matrix1[i] * matrix2[i];
////	}
////	return;
////}