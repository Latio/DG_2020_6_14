#include "BoundaryEdge.h"
#include<iostream>

double *BoundaryEdge::FToE = NULL;
double *BoundaryEdge::FToF = NULL;
double *BoundaryEdge::FToM = NULL;
double *BoundaryEdge::FToN1 = NULL;
double *BoundaryEdge::FToN2 = NULL;
double *BoundaryEdge::FToV = NULL;
signed char *BoundaryEdge::ftype = NULL;
double *BoundaryEdge::Js = NULL;
double *BoundaryEdge::LAV = NULL;
double *BoundaryEdge::M = NULL;
int *BoundaryEdge::Ne = NULL;
int *BoundaryEdge::Nfp = NULL;
double *BoundaryEdge::nx = NULL;
double *BoundaryEdge::ny = NULL;
double *BoundaryEdge::nz = NULL;
double *BoundaryEdge::r = NULL;
double *BoundaryEdge::xb = NULL;
double *BoundaryEdge::yb = NULL;

BoundaryEdge::BoundaryEdge() :bcell()
{
	MeshUnion_dim::ncvar_read(FToE, "BoundaryEdge_FToE", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "BoundaryEdge_FToF", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "BoundaryEdge_FToM", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToN1, "BoundaryEdge_FToN1", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToN2, "BoundaryEdge_FToN2", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToV, "BoundaryEdge_FToV", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(ftype, "BoundaryEdge_ftype", MeshUnion_dim::one, MeshUnion_dim::Ne_boundary);
	MeshUnion_dim::ncvar_read(Js, "BoundaryEdge_Js", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(LAV, "BoundaryEdge_LAV", MeshUnion_dim::Ne_boundary, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "BoundaryEdge_M", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ne, "BoundaryEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "BoundaryEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "BoundaryEdge_nx", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(ny, "BoundaryEdge_ny", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(nz, "BoundaryEdge_nz", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(r, "BoundaryEdge_r", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(xb, "BoundaryEdge_xb", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(yb, "BoundaryEdge_yb", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
}


BoundaryEdge::~BoundaryEdge()
{
	freememory(&FToE);
	freememory(&FToF);
	freememory(&FToM);
	freememory(&FToN1);
	freememory(&FToN2);
	freememory(&FToV);
	freememory(&ftype);
	freememory(&Js);
	freememory(&LAV);
	freememory(&M);
	freememory(&Ne);
	freememory(&Nfp);
	freememory(&nx);
	freememory(&ny);
	freememory(&nz);
	freememory(&r);
	freememory(&xb);
	freememory(&yb);

	std::cout << "Îö¹¹MeshUnion_BoundaryEdge" << std::endl;
}
void BoundaryEdge::EvaluateSurfValue(double *fphys, double *fm, double *fp, int *Np, int *K, int Nfield)
{
	c_boundary_EvaluateSurfValue(FToM, FToE, FToN1, FToN2, fphys, fm, fp, Nfp, Ne, Np, K, Nfield);
}

void BoundaryEdge::EvaluateStrongFromEdgeRHS(double *invM_,  double *J, double *fluxM, double *fluxS, int *Np, int *K,  int Nfield,double *const frhs_temp)
{
	c_boundary_EvaluateStrongFormEdgeRHS(invM_, M, FToE, FToN1, Js, J, fluxM, fluxS, Np, K, Nfp, Ne, Nfield, frhs_temp);

};

