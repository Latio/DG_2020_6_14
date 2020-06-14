#include "RollerWaveRadiationSolver.h"
RollerWaveRadiationSolver::RollerWaveRadiationSolver() :TimeInterval(5), den(1000), Time1(0), Time2(0)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);


	requestmemory(&WaveField, Np, K, 2);
};
RollerWaveRadiationSolver::~RollerWaveRadiationSolver()
{
	freememory(&WaveField);
};

void RollerWaveRadiationSolver::valuateWaveRadiationRHS(double time, double *frhs, double *fphys, double *Wave, double hmin, double gra) {

	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	signed char *status = meshunion->status;
	double *rx = meshunion->rx;
	double *Dr = meshunion->cell_p->Dr;
	double *sx = meshunion->sx;
	double *Ds = meshunion->cell_p->Ds;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;

	double *frhs_2 = frhs + num;
	double *frhs_3 = frhs + 2 * num;
	double *WaveField_1 = WaveField;
	double *WaveField_2 = WaveField + num;

	double a = Time1;
	double t1 = EvaluateTimeStep(a);
	double t2 = EvaluateTimeStep(time);


	if (t1 == t2 && time > 1e-6) {
		Time1 = time;
	}
	else
	{

		double pi = 3.141592654;

		requestmemory(&WaveField, Np, K, 2);

		double *Sxx, *Sxy, *Syy, *Rxx, *Rxy, *Ryy, *L, *miu0, *n, *Ar, *Er;
		requestmemory(&Sxx, Np, K);
		requestmemory(&Sxy, Np, K);
		requestmemory(&Syy, Np, K);
		requestmemory(&Rxx, Np, K);
		requestmemory(&Rxy, Np, K);
		requestmemory(&Ryy, Np, K);
		requestmemory(&L, Np, K);
		requestmemory(&miu0, Np, K);
		requestmemory(&n, Np, K);
		requestmemory(&Ar, Np, K);
		requestmemory(&Er, Np, K);

		bool *ind1;
		requestmemory(&ind1, num);

		double *H = Wave;
		double *T = Wave + num;
		double *dir = Wave + 2 * num;
		cblas_zdscal(num, pi / 180, dir, 1);//×ª»¯dir
		double *Qb = Wave + 3 * num;
		double *h = fphys;

		for (size_t i = 0; i < num; i++)
		{
			if (H[i] > 0 && h[i] > hmin)
			{
				ind1[i] = true;
			}
			else
			{
				ind1[i] = false;
			}
		}

		for (size_t i = 0; i < num; i++)
		{
			if (ind1[i] == true)
			{
				miu0[i] = (pow((2 * pi), 2)*h[i]) / (gra*T[i] * T[i]);
				L[i] = (2 * pi*h[i] * sqrt(tanh(miu0[i]))) / (miu0[i] * (1 + miu0[i] * exp(-1.835 - 1.225*pow(miu0[i], 1.35))));
				//wave radiation term
				n[i] = 0.5*(1 + 2 * 2 * pi*h[i] / L[i] / sinh(2 * 2 * pi*h[i] / L[i]));
				Sxx[i] = (den*gra*pow(H[i], 2) / 16)*(n[i] * pow(cos(dir[i]), 2) + 0.5*(2 * n[i] - 1));
				Sxy[i] = (den*gra*pow(H[i], 2) / 16)*(0.5*n[i] * sin(2 * dir[i]));
				Sxx[i] = (den*gra*pow(H[i], 2) / 16)*(n[i] * pow(sin(dir[i]), 2) + 0.5*(2 * n[i] - 1));
				//surface roller term
				Ar[i] = 0.06 *H[i] * L[i] * Qb[i] / sqrt(2);
				Er[i] = den * Ar[i] * L[i] / (2 * pow(T[i], 2));
				Rxx[i] = 2 * Er[i] * (pow(cos(dir[i]), 2));
				Rxy[i] = 2 * Er[i] * cos(dir[i])*sin(dir[i]);
				Ryy[i] = 2 * Er[i] * (pow(sin(dir[i]), 2));

			}

		}



		double *rx_dr_sxx, *sx_ds_sxx, *ry_dr_sxy, *sy_ds_sxy, *rx_dr_rxx, *sx_ds_rxx, *ry_dr_rxy, *sy_ds_rxy;
		requestmemory(&rx_dr_sxx, Np, K);
		requestmemory(&sx_ds_sxx, Np, K);
		requestmemory(&ry_dr_sxy, Np, K);
		requestmemory(&sy_ds_sxy, Np, K);
		requestmemory(&rx_dr_rxx, Np, K);
		requestmemory(&sx_ds_rxx, Np, K);
		requestmemory(&ry_dr_rxy, Np, K);
		requestmemory(&sy_ds_rxy, Np, K);

		Evaluate_SR(rx, Dr, Sxx, rx_dr_sxx);
		Evaluate_SR(sx, Ds, Sxx, sx_ds_sxx);
		Evaluate_SR(ry, Dr, Sxy, ry_dr_sxy);
		Evaluate_SR(sy, Ds, Sxy, sy_ds_sxy);
		Evaluate_SR(rx, Dr, Rxx, rx_dr_rxx);
		Evaluate_SR(sx, Ds, Rxx, sx_ds_rxx);
		Evaluate_SR(ry, Dr, Rxy, ry_dr_rxy);
		Evaluate_SR(sy, Ds, Rxy, sy_ds_rxy);

		cblas_dcopy(num, rx_dr_sxx, 1, WaveField_1, 1);
		cblas_dscal(num, -1, WaveField_1, 1);
		cblas_daxpy(num, -1, sx_ds_sxx, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, ry_dr_sxy, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, sy_ds_sxy, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, rx_dr_rxx, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, sx_ds_rxx, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, ry_dr_rxy, 1, WaveField_1, 1);
		cblas_daxpy(num, -1, sy_ds_rxy, 1, WaveField_1, 1);

		freememory(&rx_dr_sxx);
		freememory(&sx_ds_sxx);
		freememory(&ry_dr_sxy);
		freememory(&sy_ds_sxy);
		freememory(&rx_dr_rxx);
		freememory(&sx_ds_rxx);
		freememory(&ry_dr_rxy);
		freememory(&sy_ds_rxy);

		double *rx_dr_sxy, *sx_ds_sxy, *ry_dr_syy, *sy_ds_syy, *rx_dr_rxy, *sx_ds_rxy, *ry_dr_ryy, *sy_ds_ryy;
		requestmemory(&rx_dr_sxy, Np, K);
		requestmemory(&sx_ds_sxy, Np, K);
		requestmemory(&ry_dr_syy, Np, K);
		requestmemory(&sy_ds_syy, Np, K);
		requestmemory(&rx_dr_rxy, Np, K);
		requestmemory(&sx_ds_rxy, Np, K);
		requestmemory(&ry_dr_ryy, Np, K);
		requestmemory(&sy_ds_ryy, Np, K);

		Evaluate_SR(rx, Dr, Sxy, rx_dr_sxy);
		Evaluate_SR(sx, Ds, Sxy, sx_ds_sxy);
		Evaluate_SR(ry, Dr, Syy, ry_dr_syy);
		Evaluate_SR(sy, Ds, Syy, sy_ds_syy);
		Evaluate_SR(rx, Dr, Rxy, rx_dr_rxy);
		Evaluate_SR(sx, Ds, Rxy, sx_ds_rxy);
		Evaluate_SR(ry, Dr, Ryy, ry_dr_ryy);
		Evaluate_SR(sy, Ds, Ryy, sy_ds_ryy);

		cblas_dcopy(num, rx_dr_sxy, 1, WaveField_2, 1);
		cblas_dscal(num, -1, WaveField_2, 1);
		cblas_daxpy(num, -1, sx_ds_sxy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, ry_dr_syy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, sy_ds_syy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, rx_dr_rxy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, sx_ds_rxy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, ry_dr_ryy, 1, WaveField_2, 1);
		cblas_daxpy(num, -1, sy_ds_ryy, 1, WaveField_2, 1);

		freememory(&rx_dr_sxy);
		freememory(&sx_ds_sxy);
		freememory(&ry_dr_syy);
		freememory(&sy_ds_syy);
		freememory(&rx_dr_rxy);
		freememory(&sx_ds_rxy);
		freememory(&ry_dr_ryy);
		freememory(&sy_ds_ryy);

		freememory(&Sxx);
		freememory(&Sxy);
		freememory(&Syy);
		freememory(&Rxx);
		freememory(&Rxy);
		freememory(&Ryy);
		freememory(&L);
		freememory(&miu0);
		freememory(&ind1);
		freememory(&n);
		freememory(&Ar);
		freememory(&Er);
		Time1 = time;
	}

	for (int i = 0; i < K; i++)
	{
		if (status[i] == (signed char)enumSWERegion::Wet)
		{
			for (size_t j = 0; j < Np; j++)
			{
				frhs_2[i*Np + j] = frhs_2[i*Np + j] + WaveField_1[i*Np + j] / den;
				frhs_3[i*Np + j] = frhs_3[i*Np + j] + WaveField_2[i*Np + j] / den;
			}
		}

	}


}


void RollerWaveRadiationSolver::Evaluate_SR(double *r, double *d, double *SR, double *temp)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int num = Np * K;

	multiply(d, SR, temp, Np, K, Np);
	dotmul(num, r, temp, temp);
}




int RollerWaveRadiationSolver::EvaluateTimeStep(double t)
{
	int step;
	if (t < 1e-3)
	{
		step = 1;
	}
	else
	{
		step = ceil(t / TimeInterval);
	}
	return  step;
}