#include "NdgSWEHorizSmagrinskyDiffSolver.h"

NdgSWEHorizSmagrinskyDiffSolver::NdgSWEHorizSmagrinskyDiffSolver(double c_) :C(c_)
{ 
}

NdgSWEHorizSmagrinskyDiffSolver::~NdgSWEHorizSmagrinskyDiffSolver()
{
};

void NdgSWEHorizSmagrinskyDiffSolver::EvaluateDiffRHS(double *fphys, double *frhs, double *InnerEdgefm2d, double * InnerEdgefp2d, double * BoundaryEdgefm2d, double * BoundaryEdgefp2d)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	double *fphys_1 = fphys;
	double *fphys_2 = fphys + num;
	double *fphys_3 = fphys + 2 * num;

	double *frhs_1 = frhs;
	double *frhs_2 = frhs + num;
	double *frhs_3 = frhs + 2 * num;

	UpdateViscosity(fphys_2, fphys_3, fphys_1);
	UpdatePenaltyParameter(fphys_1);

	double *uv;
	requestmemory(&uv, Np, K);

	double *InnerEdgefm2d_temp, *InnerEdgefp2d_temp, *BoundaryEdgefm2d_temp, *BoundaryEdgefp2d_temp;
	requestmemory(&InnerEdgefm2d_temp, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&InnerEdgefp2d_temp, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&BoundaryEdgefm2d_temp, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&BoundaryEdgefp2d_temp, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);

	int num_inner = (*meshunion->inneredge_p->Nfp)*(*meshunion->inneredge_p->Ne);
	int num_bound = (*meshunion->boundarydge_p->Nfp) *(*meshunion->boundarydge_p->Ne);

	double *u, *v,
		*InnerEdgefm2d21, *InnerEdgefp2d21, *BoundaryEdgefm2d21, *BoundaryEdgefp2d21,
		*InnerEdgefm2d31, *InnerEdgefp2d31, *BoundaryEdgefm2d31, *BoundaryEdgefp2d31;

	requestmemory(&u, Np, K);
	requestmemory(&v, Np, K);

	requestmemory(&InnerEdgefm2d21, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&InnerEdgefp2d21, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&InnerEdgefm2d31, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&InnerEdgefp2d31, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&BoundaryEdgefm2d21, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&BoundaryEdgefp2d21, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&BoundaryEdgefm2d31, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&BoundaryEdgefp2d31, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);

	for (int i = 1; i < 3; i++)
	{
		double *InnerEdgefm2d_n = InnerEdgefm2d + i * num_inner;
		double *InnerEdgefp2d_n = InnerEdgefp2d + i * num_inner;
		double *BoundaryEdgefm2d_n = BoundaryEdgefm2d + i * num_bound;
		double *BoundaryEdgefp2d_n = BoundaryEdgefp2d + i * num_bound;



		dotdiv(num_inner, InnerEdgefm2d_n, InnerEdgefm2d, InnerEdgefm2d_temp);
		dotdiv(num_inner, InnerEdgefp2d_n, InnerEdgefp2d, InnerEdgefp2d_temp);
		dotdiv(num_bound, BoundaryEdgefm2d_n, BoundaryEdgefm2d, BoundaryEdgefm2d_temp);
		dotdiv(num_bound, BoundaryEdgefp2d_n, BoundaryEdgefp2d, BoundaryEdgefp2d_temp);

		double *fphys_n = fphys + i * num;
		dotdiv(num, fphys_n, fphys_1, uv);
		CalculateAuxialaryVariable(uv, Kappa, i, InnerEdgefm2d_temp, InnerEdgefp2d_temp, BoundaryEdgefm2d_temp, BoundaryEdgefp2d_temp);

		switch (i)
		{
		case 1:

			cblas_dcopy(num, uv, 1, u, 1);
			cblas_dcopy(num_inner, InnerEdgefm2d_temp, 1, InnerEdgefm2d21, 1);
			cblas_dcopy(num_inner, InnerEdgefp2d_temp, 1, InnerEdgefp2d21, 1);
			cblas_dcopy(num_bound, BoundaryEdgefm2d_temp, 1, BoundaryEdgefm2d21, 1);
			cblas_dcopy(num_bound, BoundaryEdgefp2d_temp, 1, BoundaryEdgefp2d21, 1);
			break;
		case 2:
			cblas_dcopy(num, uv, 1, v, 1);
			cblas_dcopy(num_inner, InnerEdgefm2d_temp, 1, InnerEdgefm2d31, 1);
			cblas_dcopy(num_inner, InnerEdgefp2d_temp, 1, InnerEdgefp2d31, 1);
			cblas_dcopy(num_bound, BoundaryEdgefm2d_temp, 1, BoundaryEdgefm2d31, 1);
			cblas_dcopy(num_bound, BoundaryEdgefp2d_temp, 1, BoundaryEdgefp2d31, 1);
			break;
		}

	}

	double *px_1 = px;
	double *px_2 = px + num;
	double *py_1 = py;
	double *py_2 = py + num;

	double *frhs_temp1, *frhs_temp2, *frhs_temp3;
	requestmemory(&frhs_temp1, Np, K);
	requestmemory(&frhs_temp2, Np, K);
	requestmemory(&frhs_temp3, Np, K);


	CalculatePartDerivTermXY(frhs_temp1, px_1, Kappa, u, 1, InnerEdgefm2d21, InnerEdgefp2d21, BoundaryEdgefm2d21, BoundaryEdgefp2d21, meshunion->rx, meshunion->sx, meshunion->inneredge_p->nx, meshunion->boundarydge_p->nx);
	CalculatePartDerivTermXY(frhs_temp2, py_1, Kappa, u, 1, InnerEdgefm2d21, InnerEdgefp2d21, BoundaryEdgefm2d21, BoundaryEdgefp2d21, meshunion->ry, meshunion->sy, meshunion->inneredge_p->ny, meshunion->boundarydge_p->ny);
	CalculatePartDerivTermXY(frhs_temp3, px_2, Kappa, v, 1, InnerEdgefm2d31, InnerEdgefp2d31, BoundaryEdgefm2d31, BoundaryEdgefp2d31, meshunion->ry, meshunion->sy, meshunion->inneredge_p->ny, meshunion->boundarydge_p->ny);

	for (size_t i = 0; i < num; i++)
	{
		frhs_2[i] = frhs_2[i] + 2 * frhs_temp1[i] + frhs_temp2[i] + frhs_temp3[i];
	}


	CalculatePartDerivTermXY(frhs_temp1, py_1, Kappa, u, 1, InnerEdgefm2d21, InnerEdgefp2d21, BoundaryEdgefm2d21, BoundaryEdgefp2d21, meshunion->rx, meshunion->sx, meshunion->inneredge_p->nx, meshunion->boundarydge_p->nx);
	CalculatePartDerivTermXY(frhs_temp2, px_2, Kappa, v, 1, InnerEdgefm2d31, InnerEdgefp2d31, BoundaryEdgefm2d31, BoundaryEdgefp2d31, meshunion->rx, meshunion->sx, meshunion->inneredge_p->nx, meshunion->boundarydge_p->nx);
	CalculatePartDerivTermXY(frhs_temp3, py_2, Kappa, v, 1, InnerEdgefm2d31, InnerEdgefp2d31, BoundaryEdgefm2d31, BoundaryEdgefp2d31, meshunion->ry, meshunion->sy, meshunion->inneredge_p->ny, meshunion->boundarydge_p->ny);

	for (size_t i = 0; i < num; i++)
	{
		frhs_3[i] = frhs_3[i] + frhs_temp1[i] + frhs_temp2[i] + 2 * frhs_temp3[i];
	}

	freememory(&frhs_temp1);
	freememory(&frhs_temp2);
	freememory(&frhs_temp3);

	freememory(&u);
	freememory(&v);
	freememory(&uv);

	freememory(&InnerEdgefm2d_temp);
	freememory(&InnerEdgefp2d_temp);
	freememory(&BoundaryEdgefm2d_temp);
	freememory(&BoundaryEdgefp2d_temp);

	freememory(&InnerEdgefm2d21);
	freememory(&InnerEdgefp2d21);
	freememory(&InnerEdgefm2d31);
	freememory(&InnerEdgefp2d31);
	freememory(&BoundaryEdgefm2d21);
	freememory(&BoundaryEdgefp2d21);
	freememory(&BoundaryEdgefm2d31);
	freememory(&BoundaryEdgefp2d31);

};

void NdgSWEHorizSmagrinskyDiffSolver::UpdateViscosity(double *hu, double *hv, double *h)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	double *LAV = meshunion->LAV;
	double *rx = meshunion->rx;
	double *ry = meshunion->ry;
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *sx = meshunion->sx;
	double *sy = meshunion->sy;

	double *hu_, *hv_;
	double *rx_dr_hu, *sx_ds_hu, *ry_dr_hu, *sy_ds_hu;
	double *rx_dr_hv, *sx_ds_hv, *ry_dr_hv, *sy_ds_hv;
	//double *temp_result;

	requestmemory(&hu_, Np, K);
	requestmemory(&hv_, Np, K);
	requestmemory(&rx_dr_hu, Np, K);
	requestmemory(&sx_ds_hu, Np, K);
	requestmemory(&ry_dr_hu, Np, K);
	requestmemory(&sy_ds_hu, Np, K);
	requestmemory(&rx_dr_hv, Np, K);
	requestmemory(&sx_ds_hv, Np, K);
	requestmemory(&ry_dr_hv, Np, K);
	requestmemory(&sy_ds_hv, Np, K);
	//requestmemory(&temp_result, Np, K);

	cblas_dcopy(num, hu, 1, hu_, 1);
	cblas_dcopy(num, hv, 1, hv_, 1);
	dotdiv(num, hu_, h, hu_);
	dotdiv(num, hv_, h, hv_);

	Evaluate_rdhuv(rx, Dr, hu_, rx_dr_hu);
	Evaluate_rdhuv(sx, Ds, hu_, sx_ds_hu);
	Evaluate_rdhuv(ry, Dr, hu_, ry_dr_hu);
	Evaluate_rdhuv(sy, Ds, hu_, sy_ds_hu);
	Evaluate_rdhuv(rx, Dr, hv_, rx_dr_hv);
	Evaluate_rdhuv(sx, Ds, hv_, sx_ds_hv);
	Evaluate_rdhuv(ry, Dr, hv_, ry_dr_hv);
	Evaluate_rdhuv(sy, Ds, hv_, sy_ds_hv);


	dotplus(num, rx_dr_hu, sx_ds_hu, rx_dr_hu);

	dotplus(num, ry_dr_hu, sy_ds_hu, ry_dr_hu);
	dotplus(num, ry_dr_hu, rx_dr_hv, ry_dr_hu);
	dotplus(num, ry_dr_hu, sx_ds_hv, ry_dr_hu);

	dotplus(num, ry_dr_hv, sy_ds_hv, ry_dr_hv);


	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < Np; j++)
		{
			nv[i*Np + j] = C * LAV[i] * (sqrt(pow(rx_dr_hu[i*Np + j], 2) + 0.5* pow(ry_dr_hu[i*Np + j], 2) + pow(ry_dr_hv[i*Np + j], 2)));
		}
	}


	freememory(&hu_);
	freememory(&hv_);
	freememory(&rx_dr_hu);
	freememory(&sx_ds_hu);
	freememory(&ry_dr_hu);
	freememory(&sy_ds_hu);
	freememory(&rx_dr_hv);
	freememory(&sx_ds_hv);
	freememory(&ry_dr_hv);
	freememory(&sy_ds_hv);
	//freememory(&temp_result);

}


