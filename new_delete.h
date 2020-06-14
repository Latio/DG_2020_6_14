#pragma once

void requestmemory(bool **meshpoint, int*dim1);
void requestmemory(bool **meshpoint, int&dim1);
void requestmemory(bool **meshpoint, int&dim1, int&dim2);

void requestmemory(double **meshpoint, int&dim1, int&dim2);
void requestmemory(double **meshpoint, int&dim1, int&dim2, int dim3);
void requestmemory(double **meshpoint, int*dim1, int*dim2);
void requestmemory(double **meshpoint, int*dim1, int*dim2, int*dim3);
void requestmemory(double **meshpoint, int*dim1, int*dim2, int dim3);
void requestmemory(double **meshpoint, int*dim1, int dim2, int dim3);
void requestmemory(double **meshpoint, int&dim1);
void requestmemory(double **meshpoint, int*dim1);

void requestmemory(int **meshpoint, int*dim1, int*dim2);
void requestmemory(int **meshpoint, int&dim1, int&dim2);
void requestmemory(int **meshpoint, int&dim1);

void freememory(double **meshpoint);

void freememory(int **meshpoint);

void freememory(signed char **meshpoint);
void freememory(bool **meshpoint);
int  GetMatrixInverse(int np, double *array, double *result);

void multiply(double *const matrix1, double *const matrix2, double *result, int M_, int N_, int L_);
void dotmul(int num, double *const matrix1, double *matrix2, double *result);
void dotdiv(int num, double *const matrix1, double *matrix2, double *result);
void dotplus(int num, double *const matrix1, double *matrix2, double *result);

//enum enumSWERegion {
//	Sponge = 3, // % sponge cell
//	Wet = 4,		//well cell(SWE)
//	Dry = 5,		//dry cell(SWE)
//	PartialWet = 6,
//	PartialWetFlood = 7,
//	PartialWetDamBreak = 8
//} enumsweregion;