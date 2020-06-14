#include "SWEConventional2d.h"



SWEConventional2d::SWEConventional2d()
{

}


SWEConventional2d::~SWEConventional2d()
{

}

void SWEConventional2d::EvaluatePostFunc(double *fphys)
{
	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	double *hc;
	double *qxc;
	double *qyc;

	requestmemory(&hc, K);
	requestmemory(&qxc, K);
	requestmemory(&qyc, K);
	const int dis = (*K)*(*Np);

	double *fphys1 = fphys;
	mesh.GetMeshAverageValue(fphys1, hc);


	double *fphys2 = fphys + dis;
	mesh.GetMeshAverageValue(fphys2, qxc);

	double *fphys3 = fphys + dis * 2;
	mesh.GetMeshAverageValue(fphys3, qyc);

	c_EvaluatePostFunc2d(hmin, fphys, hc, qxc, qyc, K, Np);

	pre_UpdateWetDryState(fphys);

	freememory(&hc);
	freememory(&qxc);
	freememory(&qyc);

	

}

void SWEConventional2d::UpdateWetDryState(double *fphys)
{
	int K = *meshunion->K;
	int Np = *meshunion->cell_p->Np;

	bool *wetflag = new bool[K];

	for (int i = 0; i < K; i++)
	{
		wetflag[i] = true;
		for (int j = 0; j < Np; j++)
		{
			if (fphys[i*Np + j] < hmin)
			{
				wetflag[i] = false;
				break;
			}
		}
	}

	signed char *status = meshunion->status;

	for (int k = 0; k < K; k++)
	{
		if (wetflag)
		{
			status[k] = (signed char)enumSWERegion::Wet;
		}
		else
		{
			status[k] = (signed char)enumSWERegion::Dry;
		}
	}


	delete[] wetflag;
};

void SWEConventional2d::pre_UpdateWetDryState(double *fphys)
{
	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	signed char *status = meshunion->status;

	c_UpdateWDWetDryState(hmin, fphys, status, Np, K, Nfield);
};