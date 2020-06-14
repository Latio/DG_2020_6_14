#include "Bcell.h"
#include<iostream>

double* Bcell::Dr = NULL;
double* Bcell::Ds = NULL;
double* Bcell::Dt = NULL;
double* Bcell::faceType = NULL;
double* Bcell::Fmask = NULL;
double* Bcell::FToV = NULL;
double* Bcell::invM = NULL;
double* Bcell::LAV = NULL;
double* Bcell::M = NULL;
double* Bcell::N = NULL;
double* Bcell::Nface = NULL;
double* Bcell::Nfp = NULL;
double* Bcell::Nfv = NULL;
double* Bcell::Np = NULL;
double* Bcell::Nq = NULL;
double* Bcell::Nv = NULL;
double* Bcell::r = NULL;
double* Bcell::rq = NULL;
double* Bcell::s = NULL;
double* Bcell::sq = NULL;
double* Bcell::t = NULL;
double* Bcell::TNfp = NULL;
double* Bcell::tq = NULL;
double* Bcell::type = NULL;
double* Bcell::V = NULL;
double* Bcell::Vq = NULL;
double* Bcell::vr = NULL;
double* Bcell::vs = NULL;
double* Bcell::vt = NULL;
double* Bcell::wq = NULL;


Bcell::Bcell()
{
	MeshUnion_dim::ncvar_read(Dr, "InnerEdge_cell_Dr", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ds, "InnerEdge_cell_Ds", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Dt, "InnerEdge_cell_Dt", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(faceType, "InnerEdge_cell_faceType", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Fmask, "InnerEdge_cell_Fmask", MeshUnion_dim::two, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToV, "InnerEdge_cell_FToV", MeshUnion_dim::two, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(invM, "InnerEdge_cell_invM", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(LAV, "InnerEdge_cell_LAV", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "InnerEdge_cell_M", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(N, "InnerEdge_cell_N", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nface, "InnerEdge_cell_Nface", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "InnerEdge_cell_Nfp", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Nfv, "InnerEdge_cell_Nfv", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Np, "InnerEdge_cell_Np", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nq, "InnerEdge_cell_Nq", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv, "InnerEdge_cell_Nv", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(r, "InnerEdge_cell_r", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(rq, "InnerEdge_cell_rq", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(s, "InnerEdge_cell_s", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(sq, "InnerEdge_cell_sq", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(t, "InnerEdge_cell_t", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(TNfp, "InnerEdge_cell_TNfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tq, "InnerEdge_cell_tq", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(type, "InnerEdge_cell_type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V, "InnerEdge_cell_V", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Vq, "InnerEdge_cell_Vq", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(vr, "InnerEdge_cell_vr", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(vs, "InnerEdge_cell_vs", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(vt, "InnerEdge_cell_vt", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(wq, "InnerEdge_cell_wq", MeshUnion_dim::one, MeshUnion_dim::Nfp);

}


Bcell::~Bcell()
{
	freememory(&Dr);
	freememory(&Ds);
	freememory(&Dt);
	freememory(&faceType);
	freememory(&Fmask);
	freememory(&FToV);
	freememory(&invM);
	freememory(&LAV);
	freememory(&M);
	freememory(&N);
	freememory(&Nface);
	freememory(&Nfp);
	freememory(&Nfv);
	freememory(&Np);
	freememory(&Nq);
	freememory(&Nv);
	freememory(&r);
	freememory(&rq);
	freememory(&s);
	freememory(&sq);
	freememory(&t);
	freememory(&TNfp);
	freememory(&tq);
	freememory(&type);
	freememory(&V);
	freememory(&Vq);
	freememory(&vr);
	freememory(&vs);
	freememory(&vt);
	freememory(&wq);

	std::cout << "Îö¹¹MeshUnion_BoundaryEdge_Bcell" << std::endl;
}
