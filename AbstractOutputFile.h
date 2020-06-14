#pragma once
#include<netcdfcpp.h>
#include<vector>
class AbstractOutputFile
{
public:
	AbstractOutputFile(std::string NCfile_name, double timeInerval, int StepPerFile);
	~AbstractOutputFile();
	void ncFile_create(int *Np, int *K, int Nvar);
	void outputIntervalResult(double& time, double *field, int Nvar, int *Np, int *K);
	void outputResult(double time, double *field, int Nvar, int *Np, int *K);

	netCDF::NcFile resultFile;
	netCDF::NcVar output_time;
	netCDF::NcVar output_fphys;
	std::string NCfile_name;

protected:
	double timePrevious;
	double timeInerval;
	int outputStep;
	int StepPerFile;


};


