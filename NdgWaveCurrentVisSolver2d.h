#pragma once
#include<iostream>
#include"MeshUnion.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class NdgWaveCurrentVisSolver2d
{
public:
	void evaluateViscosityRHS(double *fphys, double *fext, double time, double hmin);
};

