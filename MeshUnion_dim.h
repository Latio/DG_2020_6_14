#pragma once
#include<iostream>
#include<netcdfcpp.h>
#include"new_delete.h"
class MeshUnion_dim
{
public:
	MeshUnion_dim();
	~MeshUnion_dim();



	static void ncdim_read();
	static void ncvar_read(double *&meshunion_data, std::string ncvarname, int &dim1, int &dim2);
	static void ncvar_read(double *&meshunion_data, std::string ncvarname, int &dim1);
	static void ncvar_read(int *&meshunion_data, std::string ncvarname, int &dim1, int &dim2);
	static void ncvar_read(int *&meshunion_data, std::string ncvarname, int &dim1);
	static void ncvar_read(int8_t *&meshunion_data, std::string ncvarname, int &dim1, int &dim2);
	static int K;
	static int Nv;
	static int Ne_inner;
	static int Ne_boundary;
	static int Nfp;
	static int Np;
	static int cell_Nv;
	static int cell_Nq;
	static int two;
	static int one;
};

