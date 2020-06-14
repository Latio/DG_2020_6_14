#include "NdgHorizDiffSolver.h"
#define MAX(a,b)  (((a)>(b))?(a):(b))

NdgHorizDiffSolver::NdgHorizDiffSolver() {

	requestmemory(&InnerEdgeTau, meshunion->cell_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&BoundaryEdgeTau, meshunion->cell_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&Kappa, meshunion->K, meshunion->cell_p->Np);
	requestmemory(&px, meshunion->cell_p->Np, meshunion->K, 2);
	requestmemory(&py, meshunion->cell_p->Np, meshunion->K, 2);

	requestmemory(&M, meshunion->cell_p->Np, meshunion->cell_p->Np, meshunion->K);
	requestmemory(&invM, meshunion->cell_p->Np, meshunion->cell_p->Np, meshunion->K);

	assembleMassMatrix();
}

NdgHorizDiffSolver::~NdgHorizDiffSolver() {

	freememory(&InnerEdgeTau);
	freememory(&InnerEdgeTau);
	freememory(&Kappa);
	freememory(&px);
	freememory(&py);
	freememory(&M);
	freememory(&invM);

}

void NdgHorizDiffSolver::assembleMassMatrix()
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int Nq = *(meshunion->cell_p->Nq);
	double *Vq = meshunion->cell_p->Vq;
	double *J = meshunion->J;
	double *wq = meshunion->cell_p->wq;

	double *Jq, *temp;
	requestmemory(&Jq, Nq);
	requestmemory(&temp, Np, Nq);

	for (size_t i = 0; i < K; i++)
	{
		double *M_n = M + i * Np*Np;
		double *invM_n = invM + i * Np*Np;
		double *J_n = J + i * Np;
		multiply(Vq, J_n, Jq, Nq, 1, Np);

		for (size_t j = 0; j < Nq; j++)
		{
			for (size_t k = 0; k < Np; k++)
			{
				temp[j*Nq + k] = Vq[j + k * Nq] * Jq[j] * wq[j];
			}
		}

		multiply(temp, Vq, M_n, Np, Np, Nq);

		GetMatrixInverse(Np, M_n, invM_n);//¼ÆËãÄæ¾ØÕó
	}

	freememory(&Jq);
	freememory(&temp);
};


void NdgHorizDiffSolver::UpdatePenaltyParameter(/*double *HnvM, double *HnVP,*/ double *Height)
{
	// this penalty parameter is calculated as $\tau = \frac{ (D_p + 1)(D_p + d) }{d}\frac{ n_0 }{2}\frac{ A }{V}\miu$
	int Nfp = *(meshunion->inneredge_p->Nfp);
	int Ne = *(meshunion->inneredge_p->Ne);
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	double type = *(meshunion->type);
	double N = (double)*(meshunion->cell_p->N);
	double Nface = *(meshunion->cell_p->Nface);

	double *FToE = meshunion->inneredge_p->FToE;
	double *FToN1 = meshunion->inneredge_p->FToN1;
	double *FToN2 = meshunion->inneredge_p->FToN2;

	double *HnvM, *HnvP;
	requestmemory(&HnvM, Nfp, Ne);
	requestmemory(&HnvP, Nfp, Ne);

	//double *Heigh_temp;
	//requestmemory(&Heigh_temp, Np, K);

	//cblas_dcopy(num, Height, 1, Heigh_temp, 1);
	//dotmul(num, nv, Heigh_temp, Heigh_temp);

	cblas_dcopy(num, Height, 1, Kappa, 1);
	dotmul(num, nv, Kappa, Kappa);

	EvaluateSurfValue(HnvM, HnvP, FToE, FToN1, FToN2, Kappa, Nfp, Ne, Np, K);

	double *InnerEdgeA_fm, *InnerEdgeA_fp;
	requestmemory(&InnerEdgeA_fm, Ne);
	requestmemory(&InnerEdgeA_fp, Ne);

	double *inner_LAV = meshunion->inneredge_p->LAV;
	double *LAV = meshunion->LAV;
	for (size_t i = 0; i < Ne; i++)
	{
		InnerEdgeA_fm[i] = inner_LAV[i] / LAV[(int)FToE[2 * i]];
		InnerEdgeA_fp[i] = inner_LAV[i] / LAV[(int)FToE[2 * i + 1]];
	}

	double *InnerEdgeTau_fm, *InnerEdgeTau_fp;
	requestmemory(&InnerEdgeTau_fm, Nfp, Ne);
	requestmemory(&InnerEdgeTau_fp, Nfp, Ne);


	const double k = (N + 1)*(N + type) / type * Nface / 2;

	for (size_t i = 0; i < Ne; i++)
	{
		for (size_t j = 0; j < Nfp; j++)
		{
			InnerEdgeTau_fm[i*Nfp + j] = InnerEdgeA_fm[i] * k* HnvM[i*Nfp + j];
			InnerEdgeTau_fp[i*Nfp + j] = InnerEdgeA_fp[i] * k* HnvP[i*Nfp + j];
			InnerEdgeTau[i*Nfp + j] = MAX(InnerEdgeTau_fm[i*Nfp + j], InnerEdgeTau_fp[i*Nfp + j]);
		}
	}

	FToE = meshunion->boundarydge_p->FToE;
	FToN1 = meshunion->boundarydge_p->FToN1;
	FToN2 = meshunion->boundarydge_p->FToN2;
	Ne = *(meshunion->boundarydge_p->Ne);
	Nfp = *(meshunion->boundarydge_p->Nfp);

	double *bound_LAV = meshunion->boundarydge_p->LAV;

	double *BoundaryEdgeA_fm;
	requestmemory(&BoundaryEdgeA_fm, Ne);

	for (size_t i = 0; i < Ne; i++)
	{
		BoundaryEdgeA_fm[i] = bound_LAV[i] / LAV[(int)FToE[2 * i]];

	}

	double *Hnv, *Hnu;
	requestmemory(&Hnv, Nfp, Ne);
	requestmemory(&Hnu, Nfp, Ne);

	EvaluateSurfValue(Hnv, Hnu, FToE, FToN1, FToN2, Kappa, Nfp, Ne, Np, K);

	for (size_t i = 0; i < Ne; i++)
	{
		for (size_t j = 0; j < Nfp; j++) {
			BoundaryEdgeTau[i*Nfp + j] = BoundaryEdgeA_fm[i] * k*Hnv[i*Nfp + j];
		}
	}



	freememory(&HnvM);
	freememory(&HnvP);
	//freememory(&Heigh_temp);
	freememory(&InnerEdgeA_fm);
	freememory(&InnerEdgeA_fp);
	freememory(&InnerEdgeTau_fm);
	freememory(&InnerEdgeTau_fp);
	freememory(&BoundaryEdgeA_fm);
	freememory(&Hnv);
	freememory(&Hnu);
};




void NdgHorizDiffSolver::CalculateAuxialaryVariable(double *fphys, double *Kappa, int VarIndex, double *InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp)
{
	int K = *(meshunion->K);
	int Np = *(meshunion->cell_p->Np);
	int num = Np * K;
	int Ne_inner = *meshunion->inneredge_p->Ne;
	int Nfp_inner = *meshunion->inneredge_p->Nfp;
	int Ne_bound = *meshunion->boundarydge_p->Ne;
	int Nfp_bound = *meshunion->boundarydge_p->Nfp;

	double *rx = meshunion->rx;
	double *ry = meshunion->ry;
	double *sx = meshunion->sx;
	double *sy = meshunion->sy;


	double *px_n = px + (VarIndex - 1)*Np*K;
	double *py_n = py + (VarIndex - 1)*Np*K;

	evaluate_pxy(rx, sx, fphys, px_n);
	evaluate_pxy(ry, sy, fphys, py_n);

	double *KappaM, *KappaP;
	requestmemory(&KappaM, Ne_inner, Nfp_inner);
	requestmemory(&KappaP, Ne_inner, Nfp_inner);

	double *FToE = meshunion->inneredge_p->FToE;
	double *FToN1 = meshunion->inneredge_p->FToN1;
	double *FToN2 = meshunion->inneredge_p->FToN2;

	EvaluateSurfValue(KappaM, KappaP, FToE, FToN1, FToN2, Kappa, Nfp_inner, Ne_inner, Np, K);

	evaluate_inner_pxy(meshunion->inneredge_p->nx, KappaM, KappaP, InnerEdgefm, InnerEdgefp, px_n, Ne_inner, Nfp_inner);
	evaluate_inner_pxy(meshunion->inneredge_p->ny, KappaM, KappaP, InnerEdgefm, InnerEdgefp, py_n, Ne_inner, Nfp_inner);

	freememory(&KappaM);
	freememory(&KappaP);

	requestmemory(&KappaM, Ne_bound, Nfp_bound);
	requestmemory(&KappaP, Ne_bound, Nfp_bound);

	FToE = meshunion->boundarydge_p->FToE;
	FToN1 = meshunion->boundarydge_p->FToN1;
	FToN2 = meshunion->boundarydge_p->FToN2;

	EvaluateSurfValue(KappaM, KappaP, FToE, FToN1, FToN2, Kappa, Nfp_bound, Ne_bound, Np, K);

	evaluate_bound_pxy(meshunion->boundarydge_p->nx, KappaM, BoundaryEdgefm, BoundaryEdgefp, px_n, Ne_bound, Nfp_bound);
	evaluate_bound_pxy(meshunion->boundarydge_p->ny, KappaM, BoundaryEdgefm, BoundaryEdgefp, py_n, Ne_bound, Nfp_bound);


	freememory(&KappaM);
	freememory(&KappaP);
};

//obj.px(:, k, VarIndex) = obj.invM(:, : , k) * diag(Kappa(:, k)) *...
//obj.M(:, : , k) * (physClass.meshUnion(1).rx(:, k) .* (physClass.meshUnion(1).cell.Dr *  fphys(:, k)) + ...
//physClass.meshUnion(1).sx(:, k) .* (physClass.meshUnion(1).cell.Ds *  fphys(:, k)));


void NdgHorizDiffSolver::evaluate_pxy(double *r, double *s, double *fphys, double *pxy) {

	int K = *(meshunion->K);
	int Np = *(meshunion->cell_p->Np);
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *invM_n;
	double *Kappa_n;
	double *M_n;
	double *r_n;
	double *s_n;
	double *p;
	double *fphys_n;
	double *temp1, *temp2, *temp3;
	requestmemory(&temp1, Np, Np);
	requestmemory(&temp2, Np);
	requestmemory(&temp3, Np);


	for (size_t i = 0; i < K; i++)
	{
		p = pxy + i * Np;
		invM_n = invM + i * Np*Np;
		Kappa_n = Kappa + i * Np;
		M_n = M + i * Np*Np;
		r_n = r + i * Np;
		s_n = s + i * Np;
		fphys_n = fphys + i * Np;

		multiply(Dr, fphys_n, temp2, Np, 1, Np);
		multiply(Ds, fphys_n, temp3, Np, 1, Np);

		for (size_t j = 0; j < Np; j++)
		{
			for (size_t k = 0; k < Np; k++)
			{
				temp1[j*Np + k] = invM_n[j*Np + k] * Kappa_n[j] * M_n[j*Np + k];
			}
			temp2[j] = r_n[j] * temp2[j] + s_n[j] * temp3[j];
		}

		multiply(temp1, temp2, p, Np, 1, Np);

	}

	freememory(&temp1);
	freememory(&temp2);
	freememory(&temp3);

};


void NdgHorizDiffSolver::evaluate_inner_pxy(double *n, double *KappaM, double *KappaP, double *fm, double *fp, double *pxy_n, int Ne, int Nfp) {

	int K = *(meshunion->K);
	int Np = *(meshunion->cell_p->Np);
	int num = Np * K;
	int Ne_inner = Ne;
	int Nfp_inner = Nfp;

	double *frhs_temp;
	requestmemory(&frhs_temp, Np, K);
	double *KappaM_inner_nx, *KappaP_inner_nx, *KappaMP_inner_nx;
	requestmemory(&KappaM_inner_nx, Ne_inner, Nfp_inner);
	requestmemory(&KappaP_inner_nx, Ne_inner, Nfp_inner);
	requestmemory(&KappaMP_inner_nx, Ne_inner, Nfp_inner);

	int num_inner = Ne_inner * Nfp_inner;
	dotmul(num_inner, KappaM, fm, KappaM_inner_nx);
	dotmul(num_inner, KappaP, fp, KappaP_inner_nx);
	dotplus(num_inner, KappaM_inner_nx, KappaP_inner_nx, KappaMP_inner_nx);

	double *nx_inner = n;
	//double *ny_inner = meshunion->inneredge_p->ny;

	dotmul(num_inner, KappaM_inner_nx, nx_inner, KappaM_inner_nx);
	dotmul(num_inner, KappaP_inner_nx, nx_inner, KappaP_inner_nx);
	dotmul(num_inner, KappaMP_inner_nx, nx_inner, KappaMP_inner_nx);
	cblas_dscal(num_inner, 0.5, KappaMP_inner_nx, 1);

	mesh.inneredge.EvaluateStrongFromEdgeRHS(KappaM_inner_nx, KappaP_inner_nx, KappaMP_inner_nx, frhs_temp, meshunion->cell_p->invM, meshunion->J, meshunion->cell_p->Np, meshunion->K, 1);


	for (size_t i = 0; i < num; i++)
	{
		pxy_n[i] = pxy_n[i] - frhs_temp[i];
	}

	freememory(&frhs_temp);
	freememory(&KappaM_inner_nx);
	freememory(&KappaP_inner_nx);
	freememory(&KappaMP_inner_nx);
}


void NdgHorizDiffSolver::evaluate_bound_pxy(double *n, double *KappaM, double *fm, double *fp, double *pxy_n, int Ne, int Nfp) {

	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	int num = (*Np) * (*K);
	int Ne_bound = Ne;
	int Nfp_bound = Nfp;

	double *frhs_temp;
	requestmemory(&frhs_temp, Np, K);

	double *KappaM_fm, *KappaM_fp, *KappaM_fmp;

	requestmemory(&KappaM_fm, Ne_bound, Nfp_bound);
	requestmemory(&KappaM_fp, Ne_bound, Nfp_bound);
	requestmemory(&KappaM_fmp, Ne_bound, Nfp_bound);


	int num_bound = Ne_bound * Nfp_bound;
	dotmul(num_bound, KappaM, fm, KappaM_fm);
	dotmul(num_bound, KappaM, fp, KappaM_fp);
	dotplus(num_bound, KappaM_fm, KappaM_fp, KappaM_fmp);

	double *nx_bound = n;

	dotmul(num_bound, KappaM_fm, nx_bound, KappaM_fm);
	dotmul(num_bound, KappaM_fmp, nx_bound, KappaM_fmp);

	cblas_dscal(num_bound, 0.5, KappaM_fmp, 1);

	mesh.boundarydge.EvaluateStrongFromEdgeRHS(meshunion->cell_p->invM, meshunion->J, KappaM_fm, KappaM_fmp, Np, K, 1, frhs_temp);


	for (size_t i = 0; i < num; i++)
	{
		pxy_n[i] = pxy_n[i] - frhs_temp[i];
	}

	freememory(&frhs_temp);
	freememory(&KappaM_fm);
	freememory(&KappaM_fp);
	freememory(&KappaM_fmp);
}


void NdgHorizDiffSolver::CalculatePartDerivTermXY(double *frhs_tmp, double * pfield, double * Kappa, double * fphys, int Prantl, double * InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp, double *r, double *s, double *n_inner, double *n_bound)
{
	int K = *meshunion->K;
	int Np = *meshunion->cell_p->Np;
	int num = Np * K;
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	//double *rx = meshunion->rx;
	//double *sx = meshunion->sx;

	double *rx_dr_pfield, *sx_ds_pfield;
	requestmemory(&rx_dr_pfield, Np, K);
	requestmemory(&sx_ds_pfield, Np, K);

	multiply(Dr, pfield, rx_dr_pfield, Np, K, Np);
	multiply(Ds, pfield, sx_ds_pfield, Np, K, Np);
	dotmul(num, r, rx_dr_pfield, rx_dr_pfield);
	dotmul(num, s, sx_ds_pfield, sx_ds_pfield);
	dotplus(num, rx_dr_pfield, sx_ds_pfield, frhs_tmp);

	double *LocalVariable;
	requestmemory(&LocalVariable, Np, K);

	double *temp1, *temp2, *temp3;
	requestmemory(&temp1, Np, Np);
	requestmemory(&temp2, Np);
	requestmemory(&temp3, Np);


	for (size_t i = 0; i < K; i++)
	{
		double * LocalVariable_n = LocalVariable + i * Np;
		double *invM_n = invM + i * Np*Np;
		double *Kappa_n = Kappa + i * Np;
		double *M_n = M + i * Np*Np;

		double *r_n = r + i * Np;
		double *s_n = s + i * Np;
		double *fphys_n = fphys + i * Np;

		multiply(Dr, fphys_n, temp2, Np, 1, Np);
		multiply(Ds, fphys_n, temp3, Np, 1, Np);

		for (size_t j = 0; j < Np; j++)
		{
			for (size_t k = 0; k < Np; k++)
			{
				temp1[j*Np + k] = invM_n[j*Np + k] * Kappa_n[j] * M_n[j*Np + k];
			}
			temp2[j] = r_n[j] * temp2[j] + s_n[j] * temp3[j];
		}

		multiply(temp1, temp2, LocalVariable_n, Np, 1, Np);

	}

	freememory(&temp1);
	freememory(&temp2);
	freememory(&temp3);

	double *pfieldM, *pfieldP, *fluxS;
	requestmemory(&pfieldM, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&pfieldP, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&fluxS, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);

	EvaluateSurfValue(pfieldM, pfieldP, meshunion->inneredge_p->FToE, meshunion->inneredge_p->FToN1, meshunion->inneredge_p->FToN2, pfield, *meshunion->inneredge_p->Nfp, *meshunion->inneredge_p->Ne, Np, K);
	int num_inner = (*meshunion->inneredge_p->Nfp)*(*meshunion->inneredge_p->Ne);
	dotmul(num_inner, pfieldM, n_inner, pfieldM);
	dotmul(num_inner, pfieldP, n_inner, pfieldP);

	double *LocalVariableM, *LocalVariableP;
	requestmemory(&LocalVariableM, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	requestmemory(&LocalVariableP, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne);
	EvaluateSurfValue(LocalVariableM, LocalVariableP, meshunion->inneredge_p->FToE, meshunion->inneredge_p->FToN1, meshunion->inneredge_p->FToN2, LocalVariable, *meshunion->inneredge_p->Nfp, *meshunion->inneredge_p->Ne, *meshunion->cell_p->Np, *meshunion->K);

	for (size_t i = 0; i < num_inner; i++)
	{
		fluxS[i] = n_inner[i] * ((LocalVariableM[i] + LocalVariableP[i]) / 2 - InnerEdgeTau[i] / Prantl * n_inner[i] * (InnerEdgefm[i] - InnerEdgefp[i]));
	}

	double *frhs_temp1;
	requestmemory(&frhs_temp1, Np, K);

	mesh.inneredge.EvaluateStrongFromEdgeRHS(pfieldM, pfieldP, fluxS, frhs_temp1, meshunion->cell_p->invM, meshunion->J, meshunion->cell_p->Np, meshunion->K, 1);

	//for (size_t i = 0; i < num; i++)
	//{
	//	frhs_tmp[i] = frhs_tmp[i] - frhs_temp1[i];
	//}

	freememory(&pfieldM);
	freememory(&pfieldP);
	freememory(&LocalVariableM);
	freememory(&LocalVariableP);
	freememory(&fluxS);

	requestmemory(&pfieldM, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&pfieldP, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);

	EvaluateSurfValue(pfieldM, pfieldP, meshunion->boundarydge_p->FToE, meshunion->boundarydge_p->FToN1, meshunion->boundarydge_p->FToN2, pfield, *meshunion->boundarydge_p->Nfp, *meshunion->boundarydge_p->Ne, Np, K);
	int num_bound = (*meshunion->boundarydge_p->Nfp)*(*meshunion->boundarydge_p->Ne);
	dotmul(num_bound, pfieldM, n_bound, pfieldM);

	requestmemory(&LocalVariableM, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&LocalVariableP, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	requestmemory(&fluxS, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne);
	EvaluateSurfValue(LocalVariableM, LocalVariableP, meshunion->boundarydge_p->FToE, meshunion->boundarydge_p->FToN1, meshunion->boundarydge_p->FToN2, LocalVariable, *meshunion->boundarydge_p->Nfp, *meshunion->boundarydge_p->Ne, Np, K);

	for (size_t i = 0; i < num_bound; i++)
	{
		fluxS[i] = n_bound[i] * (LocalVariableM[i] - BoundaryEdgeTau[i] / Prantl * n_bound[i] * (BoundaryEdgefm[i] - BoundaryEdgefp[i]));
	}
	double *frhs_temp2;
	requestmemory(&frhs_temp2, Np, K);

	mesh.boundarydge.EvaluateStrongFromEdgeRHS(meshunion->cell_p->invM, meshunion->J, pfieldM, fluxS, meshunion->cell_p->Np, meshunion->K, 1, frhs_temp2);

	for (size_t i = 0; i < num; i++)
	{
		frhs_tmp[i] = frhs_tmp[i] - frhs_temp1[i] - frhs_temp2[i];
	}


	freememory(&frhs_temp2);
	freememory(&frhs_temp1);
	freememory(&pfieldM);
	freememory(&pfieldP);
	freememory(&fluxS);
	freememory(&LocalVariableM);
	freememory(&LocalVariableP);

	freememory(&LocalVariable);
	freememory(&rx_dr_pfield);
	freememory(&sx_ds_pfield);

}
//void NdgHorizDiffSolver::CalculatePartDerivTermY(double *frhs_tmp, double * pfield, double * Kappa, double * fphys, int Prantl, double * InnerEdgefm, double * InnerEdgefp, double * BoundaryEdgefm, double * BoundaryEdgefp)
//{
//
//}

