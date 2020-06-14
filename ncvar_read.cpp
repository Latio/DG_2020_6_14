//#include<iostream>
//#include <netcdfcpp.h>
////#include"ncread_dim.h"
//
//using namespace std;
//using namespace netCDF;
//using namespace netCDF::exceptions;
//
//int ncvar_read(double *meshunion_data,string ncvarname)
//{
//	
//
//	NcFile dataFile("meshUnion.nc", NcFile::read);
//
//	NcVar temp_v = dataFile.getVar(ncvarname);
//
//	temp_v.getVar(meshunion_data);
//
//	
//	//NcVar BoundaryEdge_FToF_v = dataFile.getVar("BoundaryEdge_FToF");
//	//double *BoundaryEdge_FToF = new double[two * Ne_boundary];
//	//BoundaryEdge_FToF_v.getVar(BoundaryEdge_FToF);
//
//
//	//for (int i = 0; i < Nfp*Nfp; i++)
//	//{
//	//	cout << *(BoundaryEdge_M++) << " , ";
//	//}
//	//cout << endl;
//
//	//for (int i = 0; i < two * Ne_boundary; i++)
//	//{
//	//	cout << *(BoundaryEdge_FToF++) << " , ";
//	//}
//	//cout << endl;
//
//	//system("pause");
//	//return 0;
//}