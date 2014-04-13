//============================================================================
// Name        : testSettingsReader.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include <levmar.h>

//#include <omp.h>

#include "SettingsReader.h"
#include "IntensityCalculator.h"

using namespace std;

void modelFunction(double *x, double *f, int n, int p, void *data);
void printInfo(const double * info);

int main()
{
	SettingsReader settings;

	if(!settings.readSettings())
	{
		std::cerr << "Exit.";
		return -1;
	}

	double x[settings.getNbFitParameters()];
	double xlb[settings.getNbFitParameters()];
	double xub[settings.getNbFitParameters()];
	double f[settings.getNbDataPoints()];

	settings.initializeFitParameters(x);
	settings.initializeFitParametersBounds(xlb, xub);
	modelFunction(x, f,
					settings.getNbDataPoints(),
					settings.getNbFitParameters(),
					&settings);
	settings.saveFitData(f, "i");

	if (settings.getNbIterations() > 0)
	{
		settings.initializeFitData(f);

		double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
		opts[0] = LM_INIT_MU;
		opts[1] = 1E-15;
		opts[2] = 1E-15;
		opts[3] = 1E-4;
		opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used
		double covar[settings.getNbFitParameters()
				* settings.getNbFitParameters()];

		dlevmar_bc_dif(modelFunction, x, f, settings.getNbFitParameters(),
				settings.getNbDataPoints(), xlb, xub,
				NULL, //no scaling
				settings.getNbIterations(), opts, info, NULL, covar,
				(void *) &settings);
		printInfo(info);

		modelFunction(x, f, settings.getNbDataPoints(),
				settings.getNbFitParameters(), &settings);
		settings.saveFitParameters(x, covar);
		settings.saveFitData(f, "f");
	}
	std::cout << "End.";
	return 0;
}

void modelFunction(double *x, double *f, int n, int p, void *data)
{
	//#pragma omp parallel
	{
		SettingsReader * settings = (SettingsReader *) data;

		settings->resetFitParameters(x);

		IntensityCalculator calculators[settings->getNbReflections()];
		for(size_t irefl=0; irefl < settings->getNbReflections(); irefl++)
		{
			calculators[irefl].set(settings->getNbLayers(),
							settings->getQxlab(irefl),
							settings->getQzlab(irefl),
							settings->getPoissonRatio());
			calculators[irefl].initialize(settings->getCalculatorParameters(irefl));
		}

		//#pragma omp for
		for(size_t ipt=0;ipt<settings->getNbDataPoints(); ipt++)
		{
			f[ipt]=calculators[settings->getDataPoint(ipt).reflIndex].getIntensity(
				settings->getDataPoint(ipt).qx,
				settings->getDataPoint(ipt).qz)+
				settings->getBackground();
		};
	}
}

void printInfo(const double * info)
{
	switch(int(info[6]))
	{
	case 1:
		std::cout<<"\tSmall gradient J^T f";
		break;
	case 2:
		std::cout<<"Small Dp";
		break;
	case 3:
		std::cout<<"Max nb iterations";
		break;
	case 4:
		std::cout<<"Singular matrix";
		break;
	case 5:
		std::cout<<"No further error reduction is possible";
		break;
	case 6:
		std::cout<<"Small ||f||^2";
		break;
	default:
		std::cerr<<"Invalid parameters";
		return;
		break;
	}
	std::cout<<"In "<<info[5]<<" iterations ||f||^2 reduced from "<<sqrt(info[0])<<" to "<<sqrt(info[1]);
	std::cout<<"Number of function evaluations:"<<info[7];
	std::cout<<"Number of Jacobian evaluations:"<<info[8];
}
