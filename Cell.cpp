#include "Cell.h"
#include<netcdfcpp.h>

using namespace netCDF;


int *Cell::Nq = NULL;
int *Cell::Nv = NULL;

double *Cell::Dr = NULL;
double *Cell::Ds = NULL;
double *Cell::Dt = NULL;
double *Cell::FToV = NULL;
double *Cell::faceType = NULL;
double *Cell::Fmask = NULL;
double *Cell::invM = NULL;
double *Cell::LAV = NULL;
double *Cell::M = NULL;
int *Cell::N = NULL;
double *Cell::Nface = NULL;
int *Cell::Nfp = NULL;
double *Cell::Nfv = NULL;
int *Cell::Np = NULL;
//double *cell:: Nq = NULL;
//double *cell:: Nv = NULL;
double *Cell::r = NULL;
double *Cell::rq = NULL;
double *Cell::s = NULL;
double *Cell::sq = NULL;
double *Cell::t = NULL;
int *Cell::TNfp = NULL;
double *Cell::tq = NULL;
double *Cell::type = NULL;
double *Cell::V = NULL;
double *Cell::Vq = NULL;
double *Cell::vr = NULL;
double *Cell::vs = NULL;
double *Cell::vt = NULL;
double *Cell::wq = NULL;

Cell::Cell()
{
	MeshUnion_dim::ncvar_read(Nq, "cell_Nq", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv, "cell_Nv", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(Dr, "cell_Dr", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Ds, "cell_Ds", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Dt, "cell_Dt", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(FToV, "cell_FToV", MeshUnion_dim::cell_Nv, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(faceType, "cell_faceType", MeshUnion_dim::cell_Nv, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Fmask, "cell_Fmask", MeshUnion_dim::cell_Nv, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(invM, "cell_invM", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(LAV, "cell_LAV", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "cell_M", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(N, "cell_N", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nface, "cell_Nface", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "cell_Nfp", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(Nfv, "cell_Nfv", MeshUnion_dim::cell_Nv, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Np, "cell_Np", MeshUnion_dim::one);
	//MeshUnion_dim::ncvar_read(cell_Nq, "cell_Nq", one);
	//MeshUnion_dim::ncvar_read(cell_Nv, "cell_Nv", one);
	MeshUnion_dim::ncvar_read(r, "cell_r", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(rq, "cell_rq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(s, "cell_s", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sq, "cell_sq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(t, "cell_t", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(TNfp, "cell_TNfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tq, "cell_tq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(type, "cell_type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V, "cell_V", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Vq, "cell_Vq", MeshUnion_dim::Np, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(vr, "cell_vr", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(vs, "cell_vs", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(vt, "cell_vt", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(wq, "cell_wq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
}


Cell::~Cell()
{
	freememory(&Nq);
	freememory(&Nv);
	freememory(&Dr);
	freememory(&Ds);
	freememory(&Dt);
	freememory(&FToV);
	freememory(&faceType);
	freememory(&Fmask);
	freememory(&invM);
	freememory(&LAV);
	freememory(&M);
	freememory(&N);
	freememory(&Nface);
	freememory(&Nfp);
	freememory(&Nfv);
	freememory(&Np);
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
	std::cout << "Îö¹¹MeshUnion_Cell" << std::endl;
}
