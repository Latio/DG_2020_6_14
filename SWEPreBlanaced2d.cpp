#include "SWEPreBlanaced2d.h"



SWEPreBlanaced2d::SWEPreBlanaced2d()
{
}


SWEPreBlanaced2d::~SWEPreBlanaced2d()
{
}

void SWEPreBlanaced2d::EvaluateFlux(double *fphys, double *E, double *G)
{
	sweprebalancevolumeflux2d.evaluate(hmin, gra, fphys, E, G);

};

//void SWEPreBlanaced2d::UpdateWetDryState(double *fphys)
//{
//	int K = *meshunion->K;
//	int Np = *meshunion->cell_p->Np;
//
//	bool *wetflag = new bool[K];
//
//	for (int i = 0; i < K; i++)
//	{
//		wetflag[i] = true;
//		for (int j = 0; j < Np; j++)
//		{
//			if (fphys[i*Np + j] < hmin)
//			{
//				wetflag[i] = false;
//				break;
//			}
//		}
//	}
//
//	signed char *status = meshunion->status;
//
//	for (int k = 0; k < K; k++)
//	{
//		if (wetflag)
//		{
//			status[k] = (signed char)enumSWERegion::Wet;
//		}
//		else
//		{
//			status[k] = (signed char)enumSWERegion::Dry;
//		}
//	}
//
//
//	delete[] wetflag;
//};