#pragma once
class NonhydrostaticstandingWave2d
{
public:
	NonhydrostaticstandingWave2d();
	~NonhydrostaticstandingWave2d();
	double *dt;
	static double d;
	double *fexact;
	static double A;
	static int Lambda;

};

