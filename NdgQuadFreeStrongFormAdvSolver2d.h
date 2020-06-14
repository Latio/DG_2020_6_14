#pragma once
#include"MeshUnion.h"
#include"SWEPreBlanaced2d.h"



class NdgQuadFreeStrongFormAdvSolver2d
{
public:
	NdgQuadFreeStrongFormAdvSolver2d();
	~NdgQuadFreeStrongFormAdvSolver2d();

	void evaluateAdvectionRHS(double *fphys, double *frhs, double *fext, double *InnerEdgefm2d, double *InnerEdgefp2d, double *BoundaryEdgefm2d, double *BoundaryEdgefp2d);

	SWEAbstract2d sweabstract2d;
	SWEPreBlanaced2d swepreblanaced2d;
	/*InnerEdge inneredge;*/
};

