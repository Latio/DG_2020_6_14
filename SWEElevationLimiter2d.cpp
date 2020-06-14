#include "SWEElevationLimiter2d.h"
#include<fstream>
SWEElevationLimiter2d::SWEElevationLimiter2d()
{
	int Nv = *meshunion->Nv;
	int	Nv_cell = *meshunion->cell_p->Nv;
	int K = *meshunion->K;

	double *xc = meshunion->xc;
	double *yc = meshunion->yc;
	double *vx = meshunion->vx;
	double *vy = meshunion->vy;
	double *EToV = meshunion->EToV;

	requestmemory(&Nvc, Nv);

	double *v;
	requestmemory(&v, Nv_cell);

	for (int k = 0; k < K; k++)
	{
		double *temp_etov = EToV + k * Nv_cell;
		cblas_dcopy(Nv_cell, temp_etov, 1, v, 1);

		for (int i = 0; i < Nv_cell; i++)
		{
			Nvc[(int)v[i] - 1] = Nvc[(int)v[i] - 1] + 1;
		}
	}

	const int maxindex = cblas_idamax(Nv, Nvc, 1);
	Nvcmax = (int)Nvc[maxindex];

	freememory(&Nvc);

	requestmemory(&VToK, Nvcmax, Nv);
	requestmemory(&VToM, Nvcmax, Nv);
	requestmemory(&Nvc, Nv);

	double *ind;
	requestmemory(&ind, Nv_cell);

	for (int k = 0; k < K; k++)
	{
		double *temp_etov = EToV + k * Nv_cell;
		cblas_dcopy(Nv_cell, temp_etov, 1, v, 1);

		for (int n = 0; n < Nv_cell; n++)
		{
			ind[n] = Nvc[(int)v[n] - 1] + 1 + (v[n] - 1)*Nvcmax;
			VToK[(int)ind[n] - 1] = k + 1;
			VToM[(int)ind[n] - 1] = 1;
			Nvc[(int)v[n] - 1] = Nvc[(int)v[n] - 1] + 1;
		}

	}

	requestmemory(&VToW, Nvcmax, Nv);
	for (int n = 0; n < Nv; n++)
	{
		int nvcn = (int)Nvc[n];
		double *w;
		requestmemory(&w, nvcn);
		for (int m = 0; m < nvcn; m++)
		{
			int cellid = VToK[n*Nvcmax + m] - 1;
			int msehid = VToM[n*Nvcmax + m] - 1;
			double xc_ = xc[cellid];
			double yc_ = yc[cellid];

			w[m] = 1.0 / (pow(vx[n] - xc_, 2) + pow(vy[n] - yc_, 2));

		}

		double sum = cblas_dasum(nvcn, w, 1);

		for (int i = 0; i < nvcn; i++)
		{
			VToW[n*Nvcmax + i] = w[i] / sum;
		}


		freememory(&w);
	}

	freememory(&ind);
	freememory(&v);
}


SWEElevationLimiter2d::~SWEElevationLimiter2d()
{
	freememory(&Nvc);
	freememory(&VToK);
	freememory(&VToM);
	freememory(&VToW);

}

void SWEElevationLimiter2d::apply(double *fphys)
{

	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;
	signed char *status = meshunion->status;

	const int num = (*Np)*(*K);
	double *fphys_1 = fphys;
	double *fphys_4 = fphys + 3 * num;
	double *fphys_5 = fphys + 4 * num;
	cblas_dcopy(num, fphys_1, 1, fphys_5, 1);
	cblas_daxpy(num, 1, fphys_4, 1, fphys_5, 1);



	matLimit(fphys, 5);



	bool *ind;
	requestmemory(&ind, K);

	for (int i = 0; i < (*K); i++)
	{
		if (status[i] == (signed char)enumSWERegion::Wet) {
			ind[i] = true;
		}
		else {
			ind[i] = false;
		}
	}

	for (int i = 0; i < (*K); i++)
	{
		if (ind[i] == true)
		{
			for (int j = 0; j < (*Np); j++)
			{
				fphys_1[i*(*Np) + j] = fphys_5[i*(*Np) + j] - fphys_4[i*(*Np) + j];
			}
		}
	}
	//for (size_t j = 0; j < num; j++)
	//{
	//	std::cout << j + 1 << "  : " << fphys_1[j] << std::endl;
	//}

	//for (size_t i = 0; i < 4; i++)
	//{


	//}


	matLimit(fphys, 2);


	matLimit(fphys, 3);
	freememory(&ind);
};


void SWEElevationLimiter2d::matLimit(double *fphys, int fieldId)
{
	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	int *Nv = meshunion->Nv;
	int *Nv_cell = meshunion->cell_p->Nv;
	int *Nfp = meshunion->inneredge_p->Nfp;

	double *x = meshunion->x;
	double *y = meshunion->y;
	double *xc = meshunion->xc;
	double *yc = meshunion->yc;
	double *vx = meshunion->vx;
	double *vy = meshunion->vy;
	double *EToV = meshunion->EToV;
	double *Fmask = meshunion->cell_p->Fmask;

	double  *fvert, *fvmin, *fvmax, *cvar;
	requestmemory(&fvert, Nv);
	requestmemory(&fvmin, Nv);
	requestmemory(&fvmax, Nv);
	requestmemory(&cvar, K);

	EvaluateVertAverage(fphys, fieldId, fvert, fvmin, fvmax, cvar);

	double *fphys_fieldId = fphys + (fieldId - 1)*(*K)*(*Np);

	c_VertLimit2d(fphys_fieldId, x, y, xc, yc, vx, vy, fvert, fvmin, fvmax, cvar, EToV, Fmask, Np, Nv_cell, K, Nfp);

	freememory(&fvert);
	freememory(&fvmin);
	freememory(&fvmax);
	freememory(&cvar);
};


void SWEElevationLimiter2d::EvaluateVertAverage(double *fphys, int fieldId, double *fvert, double *fvmin, double *fvmax, double *cvar)
{
	int *Nv = meshunion->Nv;
	int Np = *meshunion->cell_p->Np;
	int K = *meshunion->K;

	double *fphys_fieldId = fphys + (fieldId - 1) * K*Np;
	mesh.GetMeshAverageValue(fphys_fieldId, cvar);

	c_EvaluateVertAverage(cvar, Nv, Nvc, VToM, VToK, VToW, fvert, fvmin, fvmax, Nvcmax);
};

void  assembleVertexCellConnect()
{

};