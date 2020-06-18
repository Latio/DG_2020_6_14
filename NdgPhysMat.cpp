#include "NdgPhysMat.h"
#define max(a, b) ((a > b) ? a : b)

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

NdgPhysMat::NdgPhysMat() :frhs(NULL),
startTime(0),
finalTime(259200),
outputIntervalNum(500),
tidalinterval(600)/*潮流数据间隔*/,
abstractoutputfile("20191208jiakuosan.nc", 259200.0 / 500.0, 500),
ndgswehorizsmagrinskydiffsolver(0.25)
{
	Np = meshunion->cell_p->Np;
	K = meshunion->K;
	boundarydge_Nfp = meshunion->boundarydge_p->Nfp;
	boundarydge_Ne = meshunion->boundarydge_p->Ne;
	//Nfield = meshunion->Nfield;
	int Nvar = 3;
	requestmemory(&InnerEdgefm2d, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne, Nvar);
	requestmemory(&InnerEdgefp2d, meshunion->inneredge_p->Nfp, meshunion->inneredge_p->Ne, Nvar);
	requestmemory(&BoundaryEdgefm2d, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne, Nvar);
	requestmemory(&BoundaryEdgefp2d, meshunion->boundarydge_p->Nfp, meshunion->boundarydge_p->Ne, Nvar);

	requestmemory(&fphys, Np, K, Nfield);
	requestmemory(&fext, boundarydge_Nfp, boundarydge_Ne, Nfield);
	requestmemory(&zGrad, Np, K, 2);

	netCDF::NcFile dataFile("init_fphys.nc", netCDF::NcFile::read);
	netCDF::NcVar fphys_v = dataFile.getVar("fphys");
	fphys_v.getVar(fphys);
	netCDF::NcVar zGrad_v = dataFile.getVar("zGrad");
	zGrad_v.getVar(zGrad);


	double *ind, *temp_ftoe1, *bot;
	requestmemory(&ind, boundarydge_Nfp, boundarydge_Ne);
	requestmemory(&temp_ftoe1, boundarydge_Ne, boundarydge_Nfp);
	requestmemory(&bot, Np, K);


	/////////////////////////////////////////////////
	for (int i = 0; i < *boundarydge_Nfp; i++)
	{
		cblas_dcopy(*boundarydge_Ne, meshunion->boundarydge_p->FToE, 2, temp_ftoe1 + i, *boundarydge_Nfp);
	}
	cblas_daxpy((*boundarydge_Ne)*(*boundarydge_Nfp), *Np, temp_ftoe1, 1, ind, 1);
	cblas_daxpy((*boundarydge_Ne)*(*boundarydge_Nfp), 1, meshunion->boundarydge_p->FToN1, 1, ind, 1);
	for (int i = 0; i < (*boundarydge_Nfp)*(*boundarydge_Ne); i++)
	{
		ind[i] = ind[i] - (*Np) - 1;
	}

	double *fext_4 = fext + 3 * (*boundarydge_Ne)*(*boundarydge_Nfp);
	double *fphys_4 = fphys + 3 * (*Np)*(*K);
	cblas_dcopy((*Np)*(*K), fphys_4, 1, bot, 1);
	for (int i = 0; i < (*boundarydge_Nfp)*(*boundarydge_Ne); i++)
	{
		fext_4[i] = bot[(int)ind[i]];
	}

	freememory(&ind);
	freememory(&temp_ftoe1);
	freememory(&bot);
	////////////////////////////////////////////////

	typedef enum {
		NdgEdgeInner = 0,
		NdgEdgeGaussEdge = 1,
		NdgEdgeSlipWall = 2,
		NdgEdgeNonSlipWall = 3,
		NdgEdgeZeroGrad = 4,
		NdgEdgeClamped = 5,
		NdgEdgeClampedDepth = 6,
		NdgEdgeClampedVel = 7,
		NdgEdgeFlather = 8,
		NdgEdgeNonLinearFlather = 9,
		NdgEdgeNonLinearFlatherFlow = 10,
		NdgEdgeNonReflectingFlux = 11
	} NdgEdgeType;

	signed char *ftype = meshunion->boundarydge_p->ftype;
	for (int i = 0; i < *boundarydge_Ne; i++)
	{
		NdgEdgeType type = (NdgEdgeType)ftype[i];
		if (ftype[i] == NdgEdgeClampedDepth)
		{
			obeindex.push_back(i);
		}
	}

	ifstream data("TideElevation.txt");//read tidal data
	if (!data.is_open())
	{
		std::cout << "Error File Path !!!" << std::endl;
		system("pause");
	}
	double point_tidal;
	while (data >> point_tidal)
		tidal.push_back(point_tidal);


	data.close();


}


NdgPhysMat::~NdgPhysMat()
{
	freememory(&fphys);
	freememory(&fext);
	freememory(&zGrad);

	freememory(&InnerEdgefm2d);
	freememory(&InnerEdgefp2d);
	freememory(&BoundaryEdgefm2d);
	freememory(&BoundaryEdgefp2d);

	std::cout << "析构NdgPhyMat" << std::endl;
}


void NdgPhysMat::matSolver()
{
	matEvaluateSSPRK22();
}


void NdgPhysMat::matEvaluateSSPRK22()
{
	clock_t begintime, endtime;
	begintime = clock();

	const int num = (*K)*(*Np)*Nvar;

	double time = startTime;
	double ftime = finalTime;

	double *fphys0;
	requestmemory(&fphys0, Np, K, Nvar);

	abstractoutputfile.ncFile_create(Np, K, Nvar);

	while (time < ftime)
	{
		double dt = UpdateTimeInterval(fphys)*0.4;

		cout << dt << endl;

		if (time + dt > ftime) {
			dt = ftime - time;
		}

		cblas_dcopy(num, fphys, 1, fphys0, 1);

		for (int intRK = 0; intRK < 2; intRK++) {

			double tloc = time + dt;
			UpdateExternalField(tloc, fphys);

			requestmemory(&frhs, Np, K, Nvar);
			EvaluateRHS(fphys, frhs, time);

			cblas_daxpy(num, dt, frhs, 1, fphys, 1);

			EvaluateLimiter(fphys);

			EvaluatePostFunc(fphys);//Update status

			freememory(&frhs);
		}

		cblas_dscal(num, 0.5, fphys, 1);
		cblas_daxpy(num, 0.5, fphys0, 1, fphys, 1);

		time = time + dt;
		UpdateOutputResult(time, fphys, Nvar);


		double timeRatio = time / ftime;
		std::cout << "____________________finished____________________: " << timeRatio << std::endl;
	}

	endtime = clock();
	std::cout << "\n\nRunning Time : " << endtime - begintime << " ms\n" << endl;

	freememory(&fphys0);
}


void NdgPhysMat::EvaluateRHS(double *fphys, double *frhs, double time)
{
	ndgquadfreestrongformadvsolver2d.evaluateAdvectionRHS(fphys, frhs, fext, InnerEdgefm2d, InnerEdgefp2d, BoundaryEdgefm2d, BoundaryEdgefp2d);
	//ndgwavecurrentvissolver2d.evaluateViscosityRHS(fphys, time);

	ndgswehorizsmagrinskydiffsolver.EvaluateDiffRHS(fphys, frhs, InnerEdgefm2d, InnerEdgefp2d, BoundaryEdgefm2d, BoundaryEdgefp2d);
	//evaluateViscosityRHS(fphys, time);

	sweabstract2d.EvaluateSourceTerm(fphys, frhs, zGrad, time);
};

//void NdgPhysMat::UpdateOutputResult(double time, double *fphys) {};
void NdgPhysMat::UpdateExternalField(double tloc, double *fphys)
{
	const int benfp = *meshunion->boundarydge_p->Nfp;
	const int bene = *meshunion->boundarydge_p->Ne;
	const int obnum = benfp * obeindex.size();

	const double delta = tidalinterval;

	int s1 = ceil(tloc / delta);//double s1 = floor(tloc / delta) + 1;
   //const int s2 = s1 + 1;
	const double alpha1 = (delta*s1 - tloc) / delta;
	double alpha2 = (tloc - delta * (s1 - 1)) / delta;

	std::vector<double> fnT;

	for (int i = 0; i < obnum; i++) {
		double temp = tidal[(s1 - 1)*obnum + i] * alpha1 + tidal[s1*obnum + i] * alpha2;
		fnT.push_back(temp);
	}

	double *fext_4 = fext + 3 * benfp * bene;
	for (int i = 0; i < obeindex.size(); i++) {
		for (int j = 0; j < benfp; j++)
		{
			fext[obeindex[i] * benfp + j] = max(fnT[i*benfp + j] - fext_4[obeindex[i] * benfp + j], 0);
		}
	}
}


void NdgPhysMat::UpdateOutputResult(double& time, double *fphys, int Nvar)
{
	abstractoutputfile.outputIntervalResult(time, fphys, Nvar, Np, K);
};

void NdgPhysMat::EvaluateLimiter(double *fphys)
{
	sweabstract2d.sweelevationlimiter2d.apply(fphys);
};

//double EvaluateTimeStep(double t, double  TimeInterval)
//{
//	double step;
//	if (t < 1e-3)
//		step = 1;
//	else
//	{
//		step = ceil(t / TimeInterval);
//	}
//	return step;
//}
//
//void NdgPhysMat::evaluateViscosityRHS(double *fphys, double time)
//{
//	int *const TNfp = meshunion->cell_p->TNfp;
//	int *const K = meshunion->K;
//	int *const Np = meshunion->cell_p->Np;
//	signed char *status = meshunion->status;
//	int *LAV = meshunion->LAV;
//
//	double Cs = 0.5;//之后移到控制字典
//	double lambda = 0.6;//之后移到控制字典
//	double time1;/////////////////////
//	double TimeInterval;/////////////////////////
//	double *Wave;/////////////////
//	double *LIFT;//////////////////
//
//	double *Wave_1 = Wave;
//	double *Wave_2 = Wave + (*Np)*(*K);
//
//
//
//	double *fm, *fp, *u, *v;
//	requestmemory(&fm, TNfp, K, Nvar);
//	requestmemory(&fp, TNfp, K, Nvar);
//	EvaluateSurfaceValueOld(fm, fp, fphys, fext);
//
//	requestmemory(&u, Np, K);
//	requestmemory(&v, Np, K);
//
//	double *fphys_2 = fphys + (*Np)*(*K);
//	double *fphys_3 = fphys + 2 * (*Np)*(*K);
//	double *fm_2 = fm + (*TNfp)*(*K);
//	double *fm_3 = fm + 2 * (*TNfp)*(*K);
//	double *fp_2 = fm + (*TNfp)*(*K);
//	double *fp_3 = fm + 2 * (*TNfp)*(*K);
//
//	for (int i = 0; i < *K; i++)
//	{
//		if (status[i] == (signed char)enumSWERegion::Wet)
//		{
//			for (int j = 0; j < *Np; j++)
//			{
//				u[i*(*Np) + j] = fphys_2[i*(*Np) + j] / fphys[i*(*Np) + j];
//				v[i*(*Np) + j] = fphys_3[i*(*Np) + j] / fphys[i*(*Np) + j];
//			}
//		}
//		else
//		{
//			for (int k = 0; k < *TNfp; k++)
//			{
//				fm_2[i*(*TNfp) + k] = -fp_2[i*(*TNfp) + k];
//				fm_3[i*(*TNfp) + k] = -fp_3[i*(*TNfp) + k];
//			}
//		}
//	}
//
//	double *qxU, *qxV, *qyU, *qyV;
//
//	requestmemory(&qxU, Np, K);
//	requestmemory(&qxV, Np, K);
//	requestmemory(&qyU, Np, K);
//	requestmemory(&qyV, Np, K);
//
//	EvaluateDerivativeX(qxU, fm, fp, hmin, u, 2);
//	EvaluateDerivativeX(qxV, fm, fp, hmin, v, 3);
//	EvaluateDerivativeY(qyU, fm, fp, hmin, u, 2);
//	EvaluateDerivativeY(qyV, fm, fp, hmin, v, 3);
//
//	double *SmagorinskyFact, *S;
//	requestmemory(&SmagorinskyFact, LAV, Np);
//	requestmemory(&S, LAV, Np);
//
//
//	for (int i = 0; i < (*LAV); i++)
//	{
//		for (int j = 0; j < *Np; j++)
//		{
//			S[i*(*Np) + j] = LAV[i];
//		}
//	}
//
//	for (int i = 0; i < (*LAV)*(*Np); i++)
//	{
//		SmagorinskyFact[i] = Cs * Cs*S[i] * pow((pow(qxU[i], 2) + 0.5*pow((qxV[i] + qyU[i]), 2) + pow(qyV[i], 2)), 0.5);
//	}
//
//	double a = time1;
//	double t1 = EvaluateTimeStep(a, TimeInterval);
//	double t2 = EvaluateTimeStep(time, TimeInterval);
//
//	double *WavrVisPara;
//	requestmemory(&WavrVisPara, Np, K);
//
//	if (t1 == t2 && time > 1e-6)
//		double Time1 = time;//Time1移动到成员变量中
//	else
//	{
//		double t = EvaluateTimeStep(time, TimeInterval);
//
//		double *H = Wave_1;//////////////////
//		double *Hrms;
//		requestmemory(&Hrms, Np, K);
//		cblas_dcopy((*Np)*(*K), H, 1, Hrms, 1);
//		cblas_dscal((*Np)*(*K), 1 / sqrt(2), Hrms, 1);
//		double *T = Wave_2;/////////////////////
//
//		double *h = fphys;
//
//		double *miu0, *L, *ub;
//		requestmemory(&miu0, Np, K);
//		requestmemory(&L, Np, K);
//		requestmemory(&ub, Np, K);
//
//		const double pi = 3.141592654;
//
//		for (size_t i = 0; i < (*Np)*(*K); i++)
//		{
//			if (Hrms[i] > 0 && h[i] > hmin)
//			{
//				//evaluate wave length
//				miu0[0] = (pow(2 * pi, 2)*h[i]) / (gra*T[i] * T[i]);
//				L[i] = (2 * pi*h[i] * sqrt(tanh(miu0[i]))) / (miu0[i] * (1 + miu0[i] * exp(-1.835 - 1.225*pow(miu0[i], 1.35))));
//
//				//evaluate WaveVisPara
//				ub[i] = gra * Hrms[i] * T[i] / (2 * L[i] * cosh(2 * pi*h[i] / L[i]));
//				WavrVisPara[i] = lambda * ub[i] * Hrms[i];
//			}
//			else
//			{
//				WavrVisPara[i] = 0;
//			}
//		}
//
//		double Time1 = time;////将来将Time1移动到成员变量中
//
//		freememory(&miu0);
//		freememory(&L);
//		freememory(&ub);
//		freememory(&Hrms);
//
//	}
//	double *VisPara;
//	requestmemory(&VisPara, Np, K);
//
//	for (size_t i = 0; i < (*Np)*(*K); i++)
//	{
//		VisPara[i] = sqrt(pow(SmagorinskyFact[i], 2) + pow(WavrVisPara[i], 2));
//	}
//
//	double *varflux;
//	requestmemory(&varflux, TNfp, K, Nvar);
//
//	for (size_t i = 0; i < (*TNfp)*(*K)*Nvar; i++)
//	{
//		varflux[i] = (fm[i] + fp[i]) / 2;
//	}
//
//	double *varflux_2 = varflux + (*TNfp)*(*K);
//	double *varflux_3 = varflux + 2 * (*TNfp)*(*K);
//
//	double *ry = meshunion->ry;
//	double *sy = meshunion->sy;
//	double *Dr = meshunion->cell_p->Dr;
//	double *Ds = meshunion->cell_p->Ds;
//	double *Js = meshunion->Js;  //
//	double *J = meshunion->J;
//	double *ny = meshunion->ny;
//	double *rx = meshunion->rx;
//	double *sx = meshunion->sx;
//	double *nx = meshunion->nx;
//
//	//Calculate  Volume 
//	double *hux, *huy, *hvx, *hvy;
//	requestmemory(&hux, Np, K);
//	requestmemory(&huy, Np, K);
//	requestmemory(&hvx, Np, K);
//	requestmemory(&hvy, Np, K);
//
//	/*double *dr_fphys2, *ds_fphys2, *dr_fphys3, *ds_fphys3;*/
//
//	Calculate_Volume(hux, rx, sx, Dr, fphys_2, Ds, LIFT, nx, Js, fm_2, varflux_2, J);
//	Calculate_Volume(huy, ry, sy, Dr, fphys_2, Ds, LIFT, ny, Js, fm_2, varflux_2, J);
//	Calculate_Volume(hvx, rx, sx, Dr, fphys_3, Ds, LIFT, nx, Js, fm_3, varflux_3, J);
//	Calculate_Volume(hvy, ry, sy, Dr, fphys_3, Ds, LIFT, ny, Js, fm_3, varflux_3, J);
//
//	double *Num_hux, *Num_huy, *Num_hvx, *Num_hvy;
//	requestmemory(&Num_hux, TNfp, K);
//	requestmemory(&Num_huy, TNfp, K);
//	requestmemory(&Num_hvx, TNfp, K);
//	requestmemory(&Num_hvy, TNfp, K);
//
//	int *eidM = meshunion->eidM;///////////////////
//	int *eidP = meshunion->eidP;///////////////////
//
//	const int num = (*TNfp)*(*K);
//
//	for (size_t i = 0; i < num; i++)
//	{
//		Num_hux[i] = (hux[eidM[i] - 1] + hux[eidP[i] - 1]) / 2.0;
//	}
//
//	for (size_t i = 0; i < num; i++)
//	{
//		Num_huy[i] = (huy[eidM[i] - 1] + huy[eidP[i] - 1]) / 2.0;
//	}
//
//	for (size_t i = 0; i < num; i++)
//	{
//		Num_hvx[i] = (hvx[eidM[i] - 1] + hvx[eidP[i] - 1]) / 2.0;
//	}
//
//	for (size_t i = 0; i < num; i++)
//	{
//		Num_hvy[i] = (hvy[eidM[i] - 1] + hvy[eidP[i] - 1]) / 2.0;
//	}
//
//	//calculate rhs
//	double *dqFlux1, *dqFlux2, *temp1, *temp2;
//	requestmemory(&dqFlux1, TNfp, K);
//	requestmemory(&dqFlux2, TNfp, K);
//	requestmemory(&temp1, TNfp, K);
//	requestmemory(&temp2, TNfp, K);
//
//	for (size_t i = 0; i < num; i++)
//	{
//		temp1[i] = nx[i] * (Num_hux[i] - hux[eidM[i] - 1]);
//	}
//	for (size_t i = 0; i < num; i++)
//	{
//		temp2[i] = ny[i] * (Num_huy[i] - huy[eidM[i] - 1]);
//	}
//	dotmul(num, temp1, temp2, dqFlux1);
//
//
//	for (size_t i = 0; i < num; i++)
//	{
//		temp1[i] = nx[i] * (Num_hvx[i] - hvx[eidM[i] - 1]);
//	}
//	for (size_t i = 0; i < num; i++)
//	{
//		temp2[i] = ny[i] * (Num_hvy[i] - hvy[eidM[i] - 1]);
//	}
//	dotmul(num, temp1, temp2, dqFlux2);
//
//	double *VisTerm_1, *VisTerm_2;
//	requestmemory(&VisTerm_1, Np, K);
//	requestmemory(&VisTerm_2, Np, K);
//	double *rx_dr_huvx, *sx_ds_huvx, *ry_dr_huvy, *sy_ds_huvy, *lift_dq1, *lift_dq2;
//	requestmemory(&rx_dr_huvx, Np, K);
//	requestmemory(&sx_ds_huvx, Np, K);
//	requestmemory(&ry_dr_huvy, Np, K);
//	requestmemory(&sy_ds_huvy, Np, K);
//	requestmemory(&lift_dq1, TNfp, K);
//	requestmemory(&lift_dq2, Np, K);
//
//	multiply(Dr, hux, rx_dr_huvx, *Np, *K, *Np);
//	dotmul((*Np)*(*K), rx, rx_dr_huvx, rx_dr_huvx);
//
//	multiply(Ds, hux, sx_ds_huvx, *Np, *K, *Np);
//	dotmul((*Np)*(*K), sx, sx_ds_huvx, sx_ds_huvx);
//
//	multiply(Dr, huy, ry_dr_huvy, *Np, *K, *Np);
//	dotmul((*Np)*(*K), ry, ry_dr_huvy, ry_dr_huvy);
//
//	multiply(Ds, huy, sy_ds_huvy, *Np, *K, *Np);
//	dotmul((*Np)*(*K), sy, sy_ds_huvy, sy_ds_huvy);
//
//	dotmul((*TNfp)*(*K), Js, dqFlux1, lift_dq1);
//	multiply(LIFT, lift_dq1, lift_dq2, *Np, *K, *TNfp);
//	dotdiv((*Np)*(*K), lift_dq2, J, lift_dq2);
//
//	for (size_t i = 0; i < (*Np)*(*K); i++)
//	{
//		VisTerm_1[i] = rx_dr_huvx[i] + sx_ds_huvx[i] + ry_dr_huvy[i] + sy_ds_huvy[i] + lift_dq2[i];
//	}
//
//	multiply(Dr, hvx, rx_dr_huvx, *Np, *K, *Np);
//	dotmul((*Np)*(*K), rx, rx_dr_huvx, rx_dr_huvx);
//
//	multiply(Ds, hvx, sx_ds_huvx, *Np, *K, *Np);
//	dotmul((*Np)*(*K), sx, sx_ds_huvx, sx_ds_huvx);
//
//	multiply(Dr, hvy, ry_dr_huvy, *Np, *K, *Np);
//	dotmul((*Np)*(*K), ry, ry_dr_huvy, ry_dr_huvy);
//
//	multiply(Ds, hvy, sy_ds_huvy, *Np, *K, *Np);
//	dotmul((*Np)*(*K), sy, sy_ds_huvy, sy_ds_huvy);
//
//	dotmul((*TNfp)*(*K), Js, dqFlux2, lift_dq1);
//	multiply(LIFT, lift_dq1, lift_dq2, *Np, *K, *TNfp);
//	dotdiv((*Np)*(*K), lift_dq2, J, lift_dq2);
//
//	for (size_t i = 0; i < (*Np)*(*K); i++)
//	{
//		VisTerm_2[i] = rx_dr_huvx[i] + sx_ds_huvx[i] + ry_dr_huvy[i] + sy_ds_huvy[i] + lift_dq2[i];
//	}
//
//
//	double *frhs_2 = frhs + (*Np)*(*K);
//	double *frhs_3 = frhs + 2 * (*Np)*(*K);
//
//	for (int i = 0; i < *K; i++)
//	{
//		if (status[i] == (signed char)enumSWERegion::Wet)
//		{
//			for (int j = 0; j < *Np; j++)
//			{
//				frhs_2[i*(*Np) + j] = frhs_2[i*(*Np) + j] + VisPara[i*(*Np) + j] * VisTerm_1[i*(*Np) + j];
//				frhs_3[i*(*Np) + j] = frhs_3[i*(*Np) + j] + VisPara[i*(*Np) + j] * VisTerm_2[i*(*Np) + j];
//			}
//		}
//
//	}
//
//
//	freememory(&rx_dr_huvx);
//	freememory(&sx_ds_huvx);
//	freememory(&ry_dr_huvy);
//	freememory(&sy_ds_huvy);
//	freememory(&lift_dq1);
//	freememory(&lift_dq2);
//
//	freememory(&temp1);
//	freememory(&temp2);
//	freememory(&dqFlux1);
//	freememory(&dqFlux2);
//
//	freememory(&Num_hux);
//	freememory(&Num_huy);
//	freememory(&Num_hvx);
//	freememory(&Num_hvy);
//
//	freememory(&hux);
//	freememory(&huy);
//	freememory(&hvx);
//	freememory(&hvy);
//
//	freememory(&varflux);
//	freememory(&VisPara);
//
//	freememory(&SmagorinskyFact);
//	freememory(&S);
//	freememory(&qxU);
//	freememory(&qxV);
//	freememory(&qyU);
//	freememory(&qyV);
//	freememory(&u);
//	freememory(&v);
//	freememory(&fm);
//	freememory(&fp);
//	freememory(&WavrVisPara);
//
//};
//
//
////signed char*status = meshunion->status;
//	//std::ofstream out("D:\\Desktop\\input.txt");
//	//if (!out)
//	//{
//	//	std::cerr << "open error!" << std::endl;
//	//}
//
//	////for (int j = 0; j < 433; j++) {
//
//	////	out << j + 1;
//	//for (size_t i = 0; i < (*K); i++)
//	//{
//	//	out << i << "    " << (int)status[i] << "\n";
//	//}
//	//out << "\n";
//
//	////}
//	//cout << "************************************************\n";
//
//
//int NdgPhysMat::Calculate_Volume(double *huvxy, double *r, double *s, double *dr, double *fphys_n, double *ds, double *lift, double *n, double *js, double *fm_n, double *varflux_n, double *j) {
//	int *Np = meshunion->cell_p->Np;
//	int *K = meshunion->K;
//	int *TNfp = meshunion->cell_p->TNfp;
//
//	const int num1 = (*Np)*(*K);
//	const int num2 = (*TNfp)*(*K);
//
//
//	double  *dr_fphys, *ds_fphys, *lift_fphys1, *lift_fphys2, *lift_fphys3;
//	requestmemory(&dr_fphys, Np, K);
//	requestmemory(&ds_fphys, Np, K);
//	requestmemory(&lift_fphys1, TNfp, K);
//	requestmemory(&lift_fphys2, TNfp, K);
//	requestmemory(&lift_fphys3, Np, K);
//
//
//	multiply(dr, fphys_n, dr_fphys, *Np, *K, *Np);
//	dotmul(num1, r, dr_fphys, dr_fphys);
//
//	multiply(ds, fphys_n, ds_fphys, *Np, *K, *Np);
//	dotmul(num1, s, ds_fphys, ds_fphys);
//
//	dotmul(num2, n, js, lift_fphys1);
//
//	for (size_t i = 0; i < num2; i++)
//	{
//		lift_fphys2[i] = varflux_n[i] - fm_n[i];
//	}
//
//	dotmul(num2, lift_fphys1, lift_fphys2, lift_fphys2);
//	multiply(lift, lift_fphys2, lift_fphys3, *Np, *K, *TNfp);
//	dotdiv(num1, lift_fphys3, j, lift_fphys3);
//
//	dotmul(num1, dr_fphys, ds_fphys, huvxy);
//	dotmul(num1, lift_fphys3, huvxy, huvxy);
//
//	freememory(&dr_fphys);
//	freememory(&ds_fphys);
//	freememory(&lift_fphys1);
//	freememory(&lift_fphys2);
//	freememory(&lift_fphys3);
//	return 0;
//}