/*
 * IntensityCalculator.h
 *
 *  Created on: 3 nov. 2011
 *      Author: kopp
 */

#ifndef INTENSITYCALCULATOR_H_
#define INTENSITYCALCULATOR_H_

#include <cmath>
#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <Log.h>
#include "RelaxedInterface.h"

class IntensityCalculator
{
public:
	enum {idxD, idxRHO0, idxK, idxSIGN, idxG1, idxG2, idxDQCZ, idxSQF, idxDQX, idxDQZ, idxNB};//indices of parameters in parameter array 
	IntensityCalculator();
	void set(size_t nbLayers, double Q0X, double Q0Z, double n);
	virtual ~IntensityCalculator();
	void initialize(const double * params);
	double getIntensity(double qx, double qz);
protected:
	friend double integrand(const double z, void * myself);
	double getIntensity(size_t ilayer);
	double Q0x, Q0z;	//reflection
	double Qx, Qz;		//reciprocal point
	double qx, qz;		//deviation from reflection

	// System parameters:
	double * d;					//coordinates of interfaces
		
	RelaxedInterface * I1; //Interfaces with 90deg dislocations
	RelaxedInterface * I2; //Interfaces with 60deg dislocations

	double * rho1, * rho2;		//dislocation densities
	double * rho0;				//hipothetical densities of 90-deg dislocations which would define the observed peaks' shifts
	double * k;					//the portion of 60-deg dislocation
	double * g1, * g2;			//correlation parameters
	double * grho1, * grho2; 	//parameters defining the g*rho ratio
	
	double * b1x, * b1y, * b1z;  //burgers vectors
	double * b2x, * b2y, * b2z;
	
	double * Qcdx, * Qcdz;		//dislocation peak shifts
	double * Qccz;				//concentration peak shifts
	double * dQx, * dQz;		//peak shifts due to some other factors
	double * Qcx, * Qcz;		//total peak shifts
	
	double * sqf;				//square modulus of structure factors
	
	double nu;				//Poisson ratio
	size_t nbLayers;
protected:
	gsl_integration_workspace * w;
};

#endif /* INTENSITYCALCULATOR_H_ */
