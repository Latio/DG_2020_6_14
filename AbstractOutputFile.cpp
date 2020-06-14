#include "AbstractOutputFile.h"
//using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


AbstractOutputFile::AbstractOutputFile(std::string NCfile_name, double timeInerval, int StepPerFile) :timePrevious(0), outputStep(0), resultFile(NCfile_name, NcFile::replace)
{
	this->timeInerval = timeInerval;
	this->StepPerFile = StepPerFile;
}


AbstractOutputFile::~AbstractOutputFile()
{
	resultFile.close();
}

void AbstractOutputFile::outputIntervalResult(double &time, double *field, int Nvar, int *Np, int *K)
{
	if ((time - timePrevious) > timeInerval)
	{
		outputResult(time, field, Nvar, Np, K);
		timePrevious = time;
	}

};

void AbstractOutputFile::outputResult(double time, double *field, int Nvar, int *Np, int *K)
{
	int start_time[1] = { outputStep };
	int count_time[1] = { 1 };
	std::vector<size_t> startInd_time(start_time, start_time + 1);
	//startInd_time.push_back(outputStep);
	std::vector<size_t> countInd_time(count_time, count_time + 1);
	//countInd_time.push_back(1);
	output_time.putVar(startInd_time, countInd_time, &time);

	int start[4] = { outputStep,0,0,0 };
	int count[4] = { 1,Nvar,*K,*Np };
	std::vector<size_t> startInd_fphys(start, start + sizeof(start) / sizeof(start[0]));
	std::vector<size_t> countInd_fphys(count, count + sizeof(count) / sizeof(count[0]));
	output_fphys.putVar(startInd_fphys, countInd_fphys, field);

	if (outputStep == StepPerFile+1)
	{
		//resultFile.close();
	}
	else
	{
		outputStep++;
	}

};

void AbstractOutputFile::ncFile_create(int *Np, int *K, int Nvar)
{
	NcDim dimNp = resultFile.addDim("Np", *Np);
	NcDim dimK = resultFile.addDim("K", *K);
	NcDim dimNfield = resultFile.addDim("Nvar", Nvar);
	NcDim dimtime = resultFile.addDim("Nt", NC_UNLIMITED);

	std::vector<NcDim> dims_time;
	dims_time.push_back(dimtime);

	std::vector<NcDim> dims_fphys;
	dims_fphys.push_back(dimtime);
	dims_fphys.push_back(dimNfield);
	dims_fphys.push_back(dimK);
	dims_fphys.push_back(dimNp);

	output_time = resultFile.addVar("time", NC_DOUBLE, dims_time);
	output_fphys = resultFile.addVar("fphys", NC_DOUBLE, dims_fphys);

	resultFile.enddef();
};