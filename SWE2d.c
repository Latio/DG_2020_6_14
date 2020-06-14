#include "SWE2d.h"
//#include "mex.h"

/** convert mex variable to PhysVolField structure */
PhysField convertMexToPhysFieldp(const double *mxfield, int *Np_, int *K_, int Nfield_) {

	PhysField field;
	field.Np = *Np_;
	field.K = *K_;
	field.Nfield = Nfield_;
	const int Ntmp = field.Np * field.K;

	field.h = mxfield;
	field.hu = field.h + Ntmp;
	field.hv = field.hu + Ntmp;
	field.z = field.hv + Ntmp;
	return field;
}
PhysField convertMexToPhysField(const double *mxfield, int Np_, int K_, int Nfield_) {

	PhysField field;
	field.Np = Np_;
	field.K = K_;
	field.Nfield = Nfield_;
	const int Ntmp = field.Np * field.K;

	field.h = mxfield;
	field.hu = field.h + Ntmp;
	field.hv = field.hu + Ntmp;
	field.z = field.hv + Ntmp;
	return field;
}

/** Evaluate the flow rate depending on the depth threshold */
void evaluateFlowRateByDeptheThreshold(
	const double hcrit,  ///< depth threshold
	const double h,      ///< depth
	const double hu,     ///< water flux
	const double hv,     ///< water flux
	double *u,           ///< result velocity
	double *v            ///< velocity
) {
	if (h > hcrit) {
		//     const double sqrt2 = 1.414213562373095;
		//     double h4 = pow(h, 4);
		//     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
		//     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
		*u = hu / h;
		*v = hv / h;
	}
	else {
		*u = 0.0;
		*v = 0.0;
	}

	return;
}

///* Evaluate whether the cell adjacent to a given cell is dry*/
//void evaluateWetDryInterface(
//	signed char *status,
//	const double *FToE,
//	double *DryFaceFlag
//) {
//	const int *dims = mxGetDimensions(FToE);
//	double *FaceToElement = mxGetPr(FToE);
//	const int Ne = dims[1];
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int i = 0; i < Ne; i++) {
//		int Local_Element = (int)FaceToElement[2 * i];
//		int Adjacent_Element = (int)FaceToElement[2 * i + 1];
//		NdgRegionType Local_type = (NdgRegionType)status[Local_Element - 1];
//		NdgRegionType Adjacent_type = (NdgRegionType)status[Adjacent_Element - 1];
//		if (Local_type != NdgRegionWet || Adjacent_type != NdgRegionWet)
//			DryFaceFlag[i] = 1;
//	}
//}