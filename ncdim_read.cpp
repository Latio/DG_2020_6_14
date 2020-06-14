//#include <netcdfcpp.h>
////#include<iostream>
////using namespace std;
//
//using namespace netCDF;
////using namespace netCDF::exceptions;
//
//int K;
//int Nv;
//int	Ne_inner;
//int Ne_boundary;
//int Nfp;
//int Np;
//int cell_Nv;
//int cell_Nq;
//int two = 2;
//int one = 1;
//
//
//void ncdim_read()
//{
//
//	NcFile dataFile("meshUnion.nc", NcFile::read);
//
//	NcVar K_v = dataFile.getVar("K");
//	NcVar Nv_v = dataFile.getVar("Nv");
//	NcVar Ne_inner_v = dataFile.getVar("InnerEdge_Ne");
//	NcVar Ne_boundary_v = dataFile.getVar("BoundaryEdge_Ne");
//	NcVar Nfp_v = dataFile.getVar("InnerEdge_Nfp");
//	NcVar Np_v = dataFile.getVar("cell_Np");
//	NcVar cell_Nv_v = dataFile.getVar("cell_Nv");
//	NcVar cell_Nq_v = dataFile.getVar("cell_Nq");
//
//	Ne_boundary_v.getVar(&Ne_boundary);
//	K_v.getVar(&K);
//	Nv_v.getVar(&Nv);
//	Ne_inner_v.getVar(&Ne_inner);
//	Ne_boundary_v.getVar(&Ne_boundary);
//	Nfp_v.getVar(&Nfp);
//	Np_v.getVar(&Np);
//	cell_Nv_v.getVar(&cell_Nv);
//	cell_Nq_v.getVar(&cell_Nq);
//
//
//	//cout << K << "," << Nv << "," << Ne_inner << "," << Ne_boundary << "," << Np << "," << cell_Nv << "," << cell_Nq << "," << two << "," << one << endl;
//
//
//	//NcVar Nfp_v = dataFile.getVar("InnerEdge_Nfp");
//	//Nfp_v.getVar(&Nfp);
//
//	//cout << Nfp << endl;
//	//double *BoundaryEdge_M = new double[Nfp*Nfp];
//
//	//NcVar BoundaryEdge_M_v = dataFile.getVar("BoundaryEdge_M");
//
//	//BoundaryEdge_M_v.getVar(BoundaryEdge_M);
//
//	//double *BoundaryEdge_FToF = new double[2 * Ne_boundary];
//
//
//	//NcVar BoundaryEdge_FToF_v = dataFile.getVar("BoundaryEdge_FToF");
//
//	//BoundaryEdge_FToF_v.getVar(BoundaryEdge_FToF);
//	
//	//system("pause");
//	//return 0;
//}
//
//
