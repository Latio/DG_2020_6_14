#include "MeshUnion.h"


using namespace netCDF;

int *MeshUnion::K = NULL;
int *MeshUnion::Nv = NULL;
int MeshUnion::Nfield = 7;


double *MeshUnion::charLength = NULL;
double *MeshUnion::EToE = NULL;
double *MeshUnion::EToF = NULL;
double *MeshUnion::EToM = NULL;
signed char *MeshUnion::EToR = NULL;
double *MeshUnion::EToV = NULL;
double *MeshUnion::ind = NULL;
double *MeshUnion::J = NULL;
//double *MeshUnion::K = NULL;
double *MeshUnion::LAV = NULL;
//double *MeshUnion::Nv = NULL;
double *MeshUnion::rx = NULL;
double *MeshUnion::ry = NULL;
double *MeshUnion::rz = NULL;
signed char *MeshUnion::status = NULL;
double *MeshUnion::sx = NULL;
double *MeshUnion::sy = NULL;
double *MeshUnion::sz = NULL;
double *MeshUnion::tx = NULL;
double *MeshUnion::ty = NULL;
double *MeshUnion::type = NULL;
double *MeshUnion::tz = NULL;
double *MeshUnion::vx = NULL;
double *MeshUnion::vy = NULL;
double *MeshUnion::vz = NULL;
double *MeshUnion::x = NULL;
double *MeshUnion::xc = NULL;
double *MeshUnion::y = NULL;
double *MeshUnion::yc = NULL;
double *MeshUnion::z = NULL;
double *MeshUnion::zc = NULL;
Cell *MeshUnion::cell_p = NULL;
BoundaryEdge *MeshUnion::boundarydge_p = NULL;
InnerEdge *MeshUnion::inneredge_p = NULL;


MeshUnion::MeshUnion() :boundarydge(), inneredge(), cell()
{
	boundarydge_p = &boundarydge;
	inneredge_p = &inneredge;
	cell_p = &cell;

	MeshUnion_dim::ncvar_read(K, "K", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv, "Nv", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(charLength, "charLength", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(EToE, "EToE", MeshUnion_dim::K, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(EToF, "EToF", MeshUnion_dim::K, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(EToM, "EToM", MeshUnion_dim::K, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(EToR, "EToR", MeshUnion_dim::one, MeshUnion_dim::K);
	MeshUnion_dim::ncvar_read(EToV, "EToV", MeshUnion_dim::K, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(ind, "ind", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(J, "J", MeshUnion_dim::K, MeshUnion_dim::Np);
	//ncvar_read(K, "K", one);
	MeshUnion_dim::ncvar_read(LAV, "LAV", MeshUnion_dim::K, MeshUnion_dim::one);
	//ncvar_read(Nv, "Nv", one);
	MeshUnion_dim::ncvar_read(rx, "rx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(ry, "ry", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(rz, "rz", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(status, "status", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(sx, "sx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sy, "sy", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sz, "sz", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(tx, "tx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(ty, "ty", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(type, "type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tz, "tz", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(vx, "vx", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(vy, "vy", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(vz, "vz", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(x, "x", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(xc, "xc", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(y, "y", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(yc, "yc", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(z, "z", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(zc, "zc", MeshUnion_dim::K, MeshUnion_dim::one);

}


MeshUnion::~MeshUnion()
{
	freememory(&charLength);
	freememory(&EToE);
	freememory(&EToF);
	freememory(&EToM);
	freememory(&EToR);
	freememory(&EToV);
	freememory(&ind);
	freememory(&J);
	freememory(&LAV);
	freememory(&rx);
	freememory(&ry);
	freememory(&rz);
	freememory(&status);
	freememory(&sx);
	freememory(&sy);
	freememory(&sz);
	freememory(&tx);
	freememory(&ty);
	freememory(&type);
	freememory(&tz);
	freememory(&vx);
	freememory(&vy);
	freememory(&vz);
	freememory(&x);
	freememory(&xc);
	freememory(&y);
	freememory(&yc);
	freememory(&z);
	freememory(&zc);

	std::cout << "Îö¹¹MeshUnion" << std::endl;
}

void MeshUnion::GetMeshAverageValue(double *nodeVal, double *averageValue)
{
	double *wq = cell_p->wq;
	double *Vq = cell_p->Vq;
	int *Np = cell_p->Np;
	int *Nq = cell_p->Nq;

	double *integralValue;
	requestmemory(&integralValue, K);
	c_GetMeshIntegralValue(nodeVal, wq, J, Vq, Np, K, Nq, integralValue);

	for (int i = 0; i < (*K); i++)
	{
		averageValue[i] = integralValue[i] / LAV[i];
	}


	freememory(&integralValue);

};