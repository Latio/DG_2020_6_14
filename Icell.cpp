#include "Icell.h"
#include<iostream>

double* Icell::Dr = NULL;
double* Icell::Ds = NULL;
double* Icell::Dt = NULL;
double* Icell::faceType = NULL;
double* Icell::Fmask = NULL;
double* Icell::FToV = NULL;
double* Icell::invM = NULL;
double* Icell::LAV = NULL;
double* Icell::M = NULL;
double* Icell::N = NULL;
double* Icell::Nface = NULL;
double* Icell::Nfp = NULL;
double* Icell::Nfv = NULL;
double* Icell::Np = NULL;
double* Icell::Nq = NULL;
double* Icell::Nv = NULL;
double* Icell::r = NULL;
double* Icell::rq = NULL;
double* Icell::s = NULL;
double* Icell::sq = NULL;
double* Icell::t = NULL;
double* Icell::TNfp = NULL;
double* Icell::tq = NULL;
double* Icell::type = NULL;
double* Icell::V = NULL;
double* Icell::Vq = NULL;
double* Icell::vr = NULL;
double* Icell::vs = NULL;
double* Icell::vt = NULL;
double* Icell::wq = NULL;

Icell::Icell():meshunion_inneredge()
{
	MeshUnion_dim::ncvar_read(Dr, "InnerEdge_cell_Dr", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ds, "InnerEdge_cell_Ds", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Dt, "InnerEdge_cell_Dt", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(faceType,"InnerEdge_cell_faceType", MeshUnion_dim::one, MeshUnion_dim::two);
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


Icell::~Icell()
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

	std::cout << "Îö¹¹MeshUnion_InnerEdge_Icell" << std::endl;
}
