#include "InnerEdge.h"



double *InnerEdge::FToE = NULL;
double *InnerEdge::FToF = NULL;
double *InnerEdge::FToM = NULL;
double *InnerEdge::FToN1 = NULL;
double *InnerEdge::FToN2 = NULL;
double *InnerEdge::FToV = NULL;
double *InnerEdge::Js = NULL;
double *InnerEdge::LAV = NULL;
double *InnerEdge::M = NULL;
int *InnerEdge::Ne = NULL;
int *InnerEdge::Nfp = NULL;
double *InnerEdge::nx = NULL;
double *InnerEdge::ny = NULL;
double *InnerEdge::nz = NULL;
double *InnerEdge::r = NULL;

InnerEdge::InnerEdge() :icell()
{
	MeshUnion_dim::ncvar_read(FToE, "InnerEdge_FToE", MeshUnion_dim::Ne_inner, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "InnerEdge_FToF", MeshUnion_dim::Ne_inner, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "InnerEdge_FToM", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN1, "InnerEdge_FToN1", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToN2, "InnerEdge_FToN2", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToV, "InnerEdge_FToV", MeshUnion_dim::Ne_inner, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Js, "InnerEdge_Js", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(LAV, "InnerEdge_LAV", MeshUnion_dim::Ne_inner, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "InnerEdge_M", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ne, "InnerEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "InnerEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "InnerEdge_nx", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(ny, "InnerEdge_ny", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(nz, "InnerEdge_nz", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(r, "InnerEdge_r", MeshUnion_dim::one, MeshUnion_dim::Nfp);
}


InnerEdge::~InnerEdge()
{
	freememory(&FToE);
	freememory(&FToF);
	freememory(&FToM);
	freememory(&FToN1);
	freememory(&FToN2);
	freememory(&FToV);
	freememory(&Js);
	freememory(&LAV);
	freememory(&M);
	freememory(&Ne);
	freememory(&Nfp);
	freememory(&nx);
	freememory(&ny);
	freememory(&nz);
	freememory(&r);
	std::cout << "Îö¹¹MeshUnion_InnerEdge" << std::endl;
}



void InnerEdge::EvaluateSurfValue(double *fphys, double *fm, double *fp, int *Np, int *K, int Nfield)
{
	c_inner_EvaluateSurfValue(FToE, FToN1, FToN2, fphys, fm, fp, Nfp, Ne, Np, K, Nfield);
}

void InnerEdge::EvaluateStrongFromEdgeRHS(double *fluxM, double *fluxP, double *fluxS, double *frhs, double *invM, double *J, int *Np, int *K, int Nfield)

{
	c_inner_EvaluateStrongFromEdgeRHS(invM, M, FToE, FToN1, FToN2, Js, J, fluxM, fluxP, fluxS, frhs, Np, K, Nfp, Ne, Nfield);
};