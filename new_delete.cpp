/*共有方法创建和释放内存*/
#include<iostream>
#include"cblas.h"
#include"new_delete.h"

void requestmemory(bool **meshpoint, int*dim1)
{
	*meshpoint = new bool[*dim1]();
};

void requestmemory(bool **meshpoint, int&dim1)
{
	*meshpoint = new bool[dim1]();
};

void requestmemory(bool **meshpoint, int&dim1, int&dim2)
{
	*meshpoint = new bool[dim1*dim2]();
};


void requestmemory(double **meshpoint, int&dim1, int&dim2)
{
	*meshpoint = new double[dim1*dim2]();
};

void requestmemory(double **meshpoint, int&dim1, int&dim2, int dim3)
{
	*meshpoint = new double[dim1*dim2*dim3]();
};


void requestmemory(double **meshpoint, int*dim1, int*dim2)
{
	*meshpoint = new double[(*dim1)*(*dim2)]();
};


void requestmemory(double **meshpoint, int*dim1, int*dim2, int*dim3)
{
	*meshpoint = new double[(*dim1)*(*dim2)*(*dim3)]();
};

void requestmemory(double **meshpoint, int*dim1, int dim2, int dim3)
{
	*meshpoint = new double[(*dim1)*(dim2)*(dim3)]();
};

void requestmemory(double **meshpoint, int*dim1, int*dim2, int dim3)
{
	*meshpoint = new double[(*dim1)*(*dim2)*dim3]();
};


void requestmemory(double **meshpoint, int&dim1)
{
	*meshpoint = new double[dim1]();
};

void requestmemory(double **meshpoint, int*dim1)
{
	*meshpoint = new double[*dim1]();
};


void requestmemory(int **meshpoint, int*dim1, int*dim2)
{
	*meshpoint = new int[(*dim1)*(*dim2)]();
};

void requestmemory(int **meshpoint, int&dim1, int&dim2)
{
	*meshpoint = new int[dim1*dim2];
};
void requestmemory(int **meshpoint, int&dim1)
{
	*meshpoint = new int[dim1]();
};



void freememory(double **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[](*meshpoint);
		*meshpoint = NULL;

	}

};

void freememory(int **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[](*meshpoint);
		*meshpoint = NULL;

	}

};

void freememory(signed char **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[] * meshpoint;
		*meshpoint = NULL;

	}

};

void freememory(bool **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[] * meshpoint;
		*meshpoint = NULL;

	}

};

void multiply(double *const matrix1, double *const matrix2, double *result, int M_, int N_, int L_)
{
	const enum CBLAS_ORDER Order = CblasColMajor;
	const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	const int M = M_;//A的行数，C的行数
	const int N = N_;//B的列数，C的列数
	const int L = L_;//A的列数，B的行数
	const double alpha = 1.0;
	const double beta = 0.0;
	const int lda = M;//A的行        
	const int ldb = L;//B的行
	const int ldc = M;//C的行   //如果列优先，分别写ABC的行

	cblas_dgemm(Order, TransA, TransB, M, N, L, alpha, matrix1, lda, matrix2, ldb, beta, result, ldc);
	//std::cout << "NdgQuadFreeStrongFormAdvSolver2d.cpp" << std::endl;
}

void dotmul(int num, double *const matrix1, double *matrix2, double *result)
{

	for (int i = 0; i < num; i++)
	{
		result[i] = matrix1[i] * matrix2[i];
	}
	return;
}

void dotdiv(int num, double *const matrix1, double *matrix2, double *result)
{

	for (int i = 0; i < num; i++)
	{
		result[i] = matrix1[i] / matrix2[i];
	}
	return;
}

void dotplus(int num, double *const matrix1, double *matrix2, double *result)
{

	for (int i = 0; i < num; i++)
	{
		result[i] = matrix1[i] + matrix2[i];
	}
	return;
}

enum enumSWERegion {
	Sponge = 3, // % sponge cell
	Wet = 4,		//well cell(SWE)
	Dry = 5,		//dry cell(SWE)
	PartialWet = 6,
	PartialWetFlood = 7,
	PartialWetDamBreak = 8
} enumsweregion;