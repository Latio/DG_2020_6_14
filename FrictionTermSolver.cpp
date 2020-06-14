#include "FrictionTermSolver.h"



FrictionTermSolver::FrictionTermSolver(double r_, double n_, double ks_, double den_, double TimeInterval_, double Time1_, double Time2_) :r(n_), n(n_), ks(ks_), den(den_), TimeInterval(TimeInterval_)/*待确定*/, Time1(Time1_), Time2(Time2_)
{

}


FrictionTermSolver::~FrictionTermSolver()
{

}

void FrictionTermSolver::EvaluateWaveFrictionCoefficient(double *fphys, double *Wave, double time, double hmin, double gra) {
	double a = Time1;
	double t1 = EvaluateTimeStep(a);
	double t2 = EvaluateTimeStep(time);
	const double pi = 3.141592654;

	if (t1 == t2 && time > 1e-6)
	{
		Time1 = time;
	}
	else
	{
		double t = EvaluateTimeStep(time);
		int Np = *(mesh.cell.Np);
		int K = *(mesh.K);


		double *H = Wave;//从swan得到的有效波高Hs
		double *T = Wave + Np * K;//从swan得到的谱峰周期PER

		int num = Np * K;
		double *Hrms;
		requestmemory(&Hrms, Np, K);
		cblas_zcopy(num, H, 1, Hrms, 1);
		cblas_zdscal(num, 1 / sqrt(2), Hrms, 1);
		double *h = fphys;//模型计算得到的实际水深
		bool *ind;
		requestmemory(&ind, num);

		for (size_t i = 0; i < num; i++)
		{
			if (Hrms[i] > 0 && h[i] > hmin)
			{
				ind[i] = true;
			}
			else
			{
				ind[i] = false;
			}
		}

		double *miu0, *L, *ub, *Ab, *fw;
		requestmemory(&miu0, Np, K);
		requestmemory(&L, Np, K);
		requestmemory(&ub, Np, K);
		requestmemory(&Ab, Np, K);
		requestmemory(&fw, Np, K);

		requestmemory(&miu0, Np, K);
		for (size_t i = 0; i < num; i++)
		{
			if (ind[i] == true)
			{
				//计算纯波浪情况下的底部切应力
				miu0[i] = (pow((2 * pi), 2)*h[i]) / (gra*T[i] * T[i]);
				L[i] = (2 * pi*h[i] * sqrt(tanh(miu0[i]))) / (miu0[i] * (1 + miu0[i] * exp(-1.835 - 1.225*pow(miu0[i], 1.35))));
				ub[i] = pi * Hrms[i] / (T[i] * sinh(2 * pi*h[i] / L[i]));
				Ab[i] = ub[i] * T[i] / (2 * pi);
				fw[i] = 1.39*pow((30 * Ab[i] / ks), -0.52);
				WaveFrictionCoefficient[i] = 0.5*den*fw[i] * ub[i] * ub[i];
			}
		}

		freememory(&miu0);
		freememory(&L);
		freememory(&ub);
		freememory(&Ab);
		freememory(&fw);
		freememory(&fw);
		freememory(&miu0);
		freememory(&Hrms);
		freememory(&ind);
	}

	Time1 = time;

}

int  FrictionTermSolver::EvaluateTimeStep(double t) {
	int step;
	if (t < 1e-3)
	{
		step = 1;
	}
	else
	{
		step = ceil(t / TimeInterval);
	}
	return step;
}

//void FrictionTermSolver::evaluateFrictionTermRHS(double gra, double hmin, double * fphys, double * frhs, double *Wave, double time)
//{
//	signed char *status = meshunion->status;
//	int Np = *(meshunion->cell_p->Np);
//	int K = *(meshunion->K);
//	double *fphys_1 = fphys;
//	double *fphys_2 = fphys + Np * K;
//	double *fphys_3 = fphys + 2 * Np * K;
//	double *frhs_1 = frhs;
//	double *frhs_2 = frhs + Np * K;
//	double *frhs_3 = frhs + 2 * Np*K;
//	int num = K * Np;
//
//	requestmemory(&WaveFrictionCoefficient, Np, K);
//
//	EvaluateWaveFrictionCoefficient(fphys, Wave, time, hmin, gra);
//
//	bool *ind, *ind2;
//	double *Tc, *Tm, *speed;
//	requestmemory(&ind, K);
//	requestmemory(&ind2, K, Np);
//
//	requestmemory(&Tc, Np, K);
//	requestmemory(&Tm, Np, K);
//	requestmemory(&speed, Np, K);
//
//	for (int i = 0; i < K; i++)
//	{
//		if (status[i] == (signed char)enumSWERegion::Wet)
//		{
//			ind[i] = true;
//		}
//		else
//		{
//			ind[i] = false;
//		}
//	}
//
//	for (size_t i = 0; i < num; i++)
//	{
//		speed[i] = pow((fphys_2[i] / fphys_1[i]), 2) + pow((fphys_3[i] / fphys_1[i]), 2);
//	}
//	for (size_t i = 0; i < K; i++)
//	{
//
//		if (ind[i] == true)
//		{
//			for (size_t j = 0; j < Np; j++)
//			{
//				Tc[i*Np + j] = den * gra*(speed[i*Np + j]) / (pow(18 * log10(12 * fphys_1[i*Np + j] / ks), 2));
//			}
//		}
//	}
//
//	for (int i = 0; i < num; i++)
//	{
//		if (abs(Tc[i]) > 0.0000001 && abs(WaveFrictionCoefficient[i]) > 0.0000001)//防止double类型和0的对比出错
//		{
//			ind2[i] = true;
//		}
//		else
//		{
//			ind2[i] = false;
//		}
//	}
//	for (size_t i = 0; i < num; i++)
//	{
//		if (ind2[i] == true)
//		{
//			Tm[i] = Tc[i] * (1 + 1.2*pow((WaveFrictionCoefficient[i] / (WaveFrictionCoefficient[i] + Tc[i])), 3.2));
//		}
//		else
//		{
//			Tm[i] = Tc[i];
//		}
//	}
//
//	bool *ind3, *temp;
//	requestmemory(&temp, K);
//	requestmemory(&ind3, K);
//
//	for (size_t i = 0; i < K; i++)
//	{
//		temp[i] = true;
//	}
//
//
//	for (size_t i = 0; i < K; i++)
//	{
//		for (size_t j = 0; j < Np; j++)
//		{
//			if (speed[i*Np + j] == 0 || std::isnan(speed[i*Np + j]))
//				temp[i] = false;
//		}
//	}
//
//
//
//	for (size_t i = 0; i < K; i++)
//	{
//		if (temp[i] && ind[i])
//		{
//			for (size_t j = 0; j < Np; j++)
//			{
//				// frhs = frhs - Tbx / rou
//				frhs_2[i*Np + j] = frhs_2[i*Np + j] - Tm[i*Np + j] * fphys_2[i*Np + j] / fphys_1[i*Np + j] / (den*sqrt(speed[i*Np + j]));
//				// frhs = frhs - Tby / rou
//				frhs_3[i*Np + j] = frhs_3[i*Np + j] - Tm[i*Np + j] * fphys_3[i*Np + j] / fphys_1[i*Np + j] / (den*sqrt(speed[i*Np + j]));
//
//			}
//		}
//
//	}
//
//
//	//for (size_t i = 0; i < K; i++)
//	//{
//	//	if (temp[i] && ind[i])
//	//	{
//	//		ind3[i] = true;
//	//	}
//	//	else
//	//	{
//	//		ind3[i] = false;
//	//	}
//	//}
//
//
//	freememory(&temp);
//	freememory(&ind3);
//	freememory(&ind);
//	freememory(&ind2);
//	freememory(&Tc);
//	freememory(&Tm);
//	freememory(&speed);
//	freememory(&WaveFrictionCoefficient);
//}





















void FrictionTermSolver::evaluateFrictionTermRHS(double gra, double hmin, double * fphys, double * frhs/*, double *Wave, double time*/)
{
	signed char *status = meshunion->status;
	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	const int num = (*K)*(*Np);
	double g = gra;

	double *fphys_1 = fphys;
	double *fphys_2 = fphys + num;
	double *fphys_3 = fphys + 2 * num;
	double *frhs_2 = frhs + num;
	double *frhs_3 = frhs + 2 * num;

	//bool *ind;
	//requestmemory(&ind, K);
	double *s;
	requestmemory(&s, Np, K);

	for (int i = 0; i < (*K); i++)
	{
		if (status[i] == (signed char)enumSWERegion::Wet) {
			for (int j = 0; j < (*Np); j++) {
				s[i*(*Np) + j] = sqrt(pow(fphys_2[i*(*Np) + j], 2) +
					pow(fphys_3[i*(*Np) + j], 2)) /
					pow(fphys_1[i*(*Np) + j], 2);

				frhs_2[i*(*Np) + j] = frhs_2[i*(*Np) + j]
					- g * r*r*(fphys_2[i*(*Np) + j] * s[i*(*Np) + j]) /
					(pow(fphys_1[i*(*Np) + j], 1.0 / 3.0));

				frhs_3[i*(*Np) + j] = frhs_3[i*(*Np) + j]
					- g * r*r*(fphys_3[i*(*Np) + j] * s[i*(*Np) + j]) /
					(pow(fphys_1[i*(*Np) + j], 1.0 / 3.0));
				//ind[i] = true;
			}
			//else {
			//	ind[i] = false;
			//}
		}



		//for (int i = 0; i < (*K); i++)
		//{
		//	if (ind[i] == true)
		//	{
		//		for (int j = 0; j < (*Np); j++)
		//		{
		//			s[i*(*Np) + j] = (pow(fphys_2[i*(*Np) + j], 2) +
		//				pow(fphys_3[i*(*Np) + j], 2)) /
		//				pow(fphys_1[i*(*Np) + j], 2);

		//			frhs_2[i*(*Np) + j] = frhs_2[i*(*Np) + j]
		//				- g * r*r*(fphys_2[i*(*Np) + j] * s[i*(*Np) + j]) /
		//				(pow(fphys_1[i*(*Np) + j], 1 / 3));

		//			frhs_3[i*(*Np) + j] = frhs_3[i*(*Np) + j]
		//				- g * r*r*(fphys_3[i*(*Np) + j] * s[i*(*Np) + j]) /
		//				(pow(fphys_1[i*(*Np) + j], 1 / 3));
		//		}
		//	}
		//}



		//freememory(&ind);

	};
	freememory(&s);
}

