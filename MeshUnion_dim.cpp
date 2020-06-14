#include "MeshUnion_dim.h"
using namespace netCDF;
using namespace netCDF::exceptions;

int MeshUnion_dim::K = 0;
int MeshUnion_dim::Nv = 0;
int	MeshUnion_dim::Ne_inner = 0;
int MeshUnion_dim::Ne_boundary = 0;
int MeshUnion_dim::Nfp = 0;
int MeshUnion_dim::Np = 0;
int MeshUnion_dim::cell_Nv = 0;
int MeshUnion_dim::cell_Nq = 0;
int MeshUnion_dim::two = 2;
int MeshUnion_dim::one = 1;

MeshUnion_dim::MeshUnion_dim()
{
	ncdim_read();
}

MeshUnion_dim::~MeshUnion_dim()
{
}

void MeshUnion_dim::ncdim_read()
{
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar K_v = dataFile.getVar("K");
	NcVar Nv_v = dataFile.getVar("Nv");
	NcVar Ne_inner_v = dataFile.getVar("InnerEdge_Ne");
	NcVar Ne_boundary_v = dataFile.getVar("BoundaryEdge_Ne");
	NcVar Nfp_v = dataFile.getVar("InnerEdge_Nfp");
	NcVar Np_v = dataFile.getVar("cell_Np");
	NcVar cell_Nv_v = dataFile.getVar("cell_Nv");
	NcVar cell_Nq_v = dataFile.getVar("cell_Nq");

	
	K_v.getVar(&K);
	Nv_v.getVar(&Nv);
	Ne_inner_v.getVar(&Ne_inner);
	Ne_boundary_v.getVar(&Ne_boundary);
	Nfp_v.getVar(&Nfp);
	Np_v.getVar(&Np);
	cell_Nv_v.getVar(&cell_Nv);
	cell_Nq_v.getVar(&cell_Nq);

}



void MeshUnion_dim::ncvar_read(double *&meshunion_data, std::string ncvarname, int &dim1, int &dim2)
{
	//meshunion_data = new double[Ne_boundary*two];
	meshunion_data = new double[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar temp_v = dataFile.getVar(ncvarname);

	temp_v.getVar(meshunion_data);
}

void MeshUnion_dim::ncvar_read(double *&meshunion_data, std::string ncvarname, int &dim1)
{
	meshunion_data = new double[dim1];
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar temp_v = dataFile.getVar(ncvarname);

	temp_v.getVar(meshunion_data);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, std::string ncvarname, int &dim1, int &dim2)
{
	//meshunion_data = new double[Ne_boundary*two];
	meshunion_data = new int[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar temp_v = dataFile.getVar(ncvarname);

	temp_v.getVar(meshunion_data);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, std::string ncvarname, int &dim1)
{
	meshunion_data = new int[dim1];
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar temp_v = dataFile.getVar(ncvarname);

	temp_v.getVar(meshunion_data);
}

void MeshUnion_dim::ncvar_read(int8_t *&meshunion_data, std::string ncvarname, int &dim1, int &dim2)
{
	//meshunion_data = new double[Ne_boundary*two];
	meshunion_data = new int8_t[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::read);

	NcVar temp_v = dataFile.getVar(ncvarname);

	temp_v.getVar(meshunion_data);
}