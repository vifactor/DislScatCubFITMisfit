/*
 * IntensityCalculator.cpp
 *
 *  Created on: 3 nov. 2011
 *      Author: kopp
 */

#include "IntensityCalculator.h"


IntensityCalculator::IntensityCalculator()
{
	nbLayers=0;
	Q0x=0; Q0z=0;
	nu=1.0/3;
	// System parameters:
	d=NULL;					//coordinates of interfaces
	
	I1=NULL; //Interfaces with 90deg dislocations
	I2=NULL; //Interfaces with 60deg dislocations
	
	rho1=NULL, rho2=NULL;		//dislocation densities
	rho0=NULL;				//hipothetical densities of 90-deg dislocations which would define the observed peaks' shifts
	k=NULL;					//the portion of 60-deg dislocation
	g1=NULL, g2=NULL;			//correlation parameters
	grho1=NULL, grho2=NULL; 	//parameters defining the g*rho ratio
	
	b1x=NULL; b1y=NULL, b1z=NULL;  //burgers vectors
	b2x=NULL; b2y=NULL, b2z=NULL;
	
	Qcdx=NULL; Qcdz=NULL;		//dislocation peak shifts
	Qccz=NULL;
	dQx=NULL; dQz=NULL;
	Qcx=NULL; Qcz=NULL;
	
	sqf=NULL;
	
	w = NULL;
	gsl_set_error_handler_off();
}

void IntensityCalculator::set(size_t nb, double qx, double qz, double n)
{
	nbLayers=nb;
	Q0x=qx; Q0z=qz;
	nu=n;
	// System parameters:
	d=new double[nbLayers];					//coordinates of interfaces
	
	I1=new RelaxedInterface[nbLayers]; //Interfaces with 90deg dislocations
	I2=new RelaxedInterface[nbLayers]; //Interfaces with 60deg dislocations
	
	rho1=new double[nbLayers]; rho2=new double[nbLayers];		//dislocation densities
	rho0=new double[nbLayers];				//hipothetical densities of 90-deg dislocations which would define the observed peaks' shifts
	k=new double[nbLayers];					//the portion of 60-deg dislocation
	g1=new double[nbLayers]; g2=new double[nbLayers];			//correlation parameters
	grho1=new double[nbLayers]; grho2=new double[nbLayers]; 	//parameters defining the g*rho ratio
	
	b1x=new double[nbLayers]; b1y=new double[nbLayers], b1z=new double[nbLayers];  //burgers vectors
	b2x=new double[nbLayers]; b2y=new double[nbLayers], b2z=new double[nbLayers];
	//the burgers vectors differs only in signs of bx coordinates, zhich will be set later
	for(size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		b1x[ilayer] = 0.707; b1y[ilayer] = 0.0; b1z[ilayer] = 0.0;
		b2x[ilayer] = 0.3536; b2y[ilayer] = 0.3536; b2z[ilayer] = 0.5;
	}
	
	Qcdx=new double[nbLayers]; Qcdz=new double[nbLayers];		//dislocation peak shifts
	Qccz=new double[nbLayers];
	dQx=new double[nbLayers]; dQz=new double[nbLayers];
	Qcx=new double[nbLayers]; Qcz=new double[nbLayers];
	
	sqf=new double[nbLayers];
	
	w = gsl_integration_workspace_alloc(5000);
	gsl_set_error_handler_off();
}

IntensityCalculator::~IntensityCalculator()
{
	if(d)
		delete[] d;					//coordinates of interfaces
	
	if(I1)
		delete[] I1;
	if(I2)
		delete[] I2;
		
	if(rho1)
		delete[] rho1;
	if(rho2)
		delete[] rho2;
	if(rho0)
		delete[] rho0;
	if(k)
		delete[] k;
	
	if(g1)
		delete[] g1;
	if(g2)
		delete[] g2;
	if(grho1)
		delete[] grho1;
	if(grho2)
		delete[] grho2;
	
	if(b1x)
		delete[] b1x;
	if(b1y)
		delete[] b1y;
	if(b1z)
		delete[] b1z;
	if(b2x)
		delete[] b2x;
	if(b2y)
		delete[] b2y;
	if(b2z)
		delete[] b2z;
		
	if(Qcdx)
		delete[] Qcdx;
	if(Qcdz)
		delete[] Qcdz;
	if(Qccz)
		delete[] Qccz;
	if(dQx)
		delete[] dQx;
	if(dQz)
		delete[] dQz;
	if(Qcx)
		delete[] Qcx; 
	if(Qcz)
		delete[] Qcz;
	
	if(sqf) 
		delete[] sqf;
	
	if(w)
		gsl_integration_workspace_free(w);
}

double integrand(const double z, void * params)
{
	IntensityCalculator * calc=reinterpret_cast<IntensityCalculator *>(params);
	double qx(0), qz(0);
	double QxQz=calc->Qx*calc->Qz, QzQz=calc->Qz*calc->Qz, QxQx=calc->Qx*calc->Qx;
	double wxx=0, wxz=0, wzz=0;
	for(size_t ilayer=0;ilayer<calc->nbLayers;ilayer++)
	{
		calc->I1[ilayer].preinitialize(z);
		calc->I2[ilayer].preinitialize(z);

		wxx+=calc->grho1[ilayer]*(QxQx*calc->I1[ilayer].Wxxxx(z)+QzQz*calc->I1[ilayer].Wzxzx(z));
		wxx+=calc->grho2[ilayer]*(QxQx*calc->I2[ilayer].Wxxxx(z)+QzQz*calc->I2[ilayer].Wzxzx(z));

		wxz+=calc->grho1[ilayer]*QxQz*(calc->I1[ilayer].Wxxzz(z)+calc->I1[ilayer].Wxzzx(z));
		wxz+=calc->grho2[ilayer]*QxQz*(calc->I2[ilayer].Wxxzz(z)+calc->I2[ilayer].Wxzzx(z));

		wzz+=calc->grho1[ilayer]*(QxQx*calc->I1[ilayer].Wxzxz(z)+2*QzQz*calc->I1[ilayer].Wzzzz(z)+QxQx*calc->I1[ilayer].Wyzyz(z));
		wzz+=calc->grho2[ilayer]*(QxQx*calc->I2[ilayer].Wxzxz(z)+2*QzQz*calc->I2[ilayer].Wzzzz(z)+QxQx*calc->I2[ilayer].Wyzyz(z));
	}

	for(size_t ilayer=0;ilayer<calc->nbLayers;ilayer++)
	{
		if(z<=calc->d[ilayer])
		{
			qx=calc->qx-calc->Qcx[ilayer];
			qz=calc->qz-calc->Qcz[ilayer];
			break;
		}
	}
	double detW=wxx*wzz-wxz*wxz;

	double w1xx=wzz/detW;
	double w1xz=-wxz/detW;
	double w1zz=wxx/detW;

	double factor=1.0/4.0*(w1xx*qx*qx+2.0*w1xz*qx*qz+w1zz*qz*qz);
	double res=M_PI/sqrt(detW)*exp(-factor);

	//LOG(logDEBUG)<<"q"<<qx<<"\t"<<qz<<std::endl;
	//LOG(logINFO)<<"Qc"<<Qcx<<"\t"<<Qcz<<std::endl;
	return res;
}

double IntensityCalculator::getIntensity(double qx, double qz)
{
	//std::cout<<"q"<<qx<<"\t"<<qz<<std::endl;
	this->qx=qx;
	this->qz=qz;

	double result=0;

	for(size_t ilayer=0;ilayer<nbLayers;ilayer++)
	{
		result += sqf[ilayer]*getIntensity(ilayer);
	}

	return result;
}

double IntensityCalculator::getIntensity(size_t i)
{
	double result(0), error(0);
	std::size_t neval;
	/*int status;*/

	gsl_function F;
	F.function = integrand;
	F.params = this;

	//status=gsl_integration_qag(&F, 0, sysParameters.d1, 0, 1e-6, 5000, GSL_INTEG_GAUSS15, w, &result, &error);
	if(i==0)
	{
		/*status=*/gsl_integration_qng(&F, 0, d[0], 0, 1e-6, &result, &error, &neval);
	}
	else
	{
		/*status=*/gsl_integration_qng(&F, d[i-1], d[i], 0, 1e-6, &result, &error, &neval);
	}
	/*if(status)
	{
		LOG(logWARNING)<<"error (0-d1): "<< gsl_strerror (status)<<"\tresult="<<result;
	}*/

	return result;
}

void IntensityCalculator::initialize(const double * params)
{
	size_t shift;
	//enum {idxD, idxRHO0, idxK, idxSIGN, idxG1, idxG2, idxDQCZ, idxSQF, idxDQX, idxDQZ, idxNB};
	for(size_t ilayer=0;ilayer<nbLayers;ilayer++)
	{
		shift=idxNB*ilayer;
		d[ilayer]=params[shift+idxD];
		rho0[ilayer]=params[shift+idxRHO0];
		k[ilayer]=params[shift+idxK];
		b1x[ilayer]*=params[shift+idxSIGN];
		b2x[ilayer]*=params[shift+idxSIGN];
		g1[ilayer]=params[shift+idxG1];
		g2[ilayer]=params[shift+idxG2];
		Qccz[ilayer]=params[shift+idxDQCZ];
		dQx[ilayer]=params[shift+idxDQX];
		dQz[ilayer]=params[shift+idxDQZ];
		sqf[ilayer]=params[shift+idxSQF];
		
		rho1[ilayer]=2.0*(1.0-k[ilayer])/(2.0-k[ilayer])*rho0[ilayer];
		rho2[ilayer]=2.0*k[ilayer]/(2.0-k[ilayer])*rho0[ilayer];
		
		grho1[ilayer]=g1[ilayer]*rho1[ilayer]/2.0;
		grho2[ilayer]=g2[ilayer]*rho2[ilayer]/2.0;
		
		
		I1[ilayer].setBurgersComponents(b1x[ilayer],
								b1y[ilayer],
								b1z[ilayer]);
		I2[ilayer].setBurgersComponents(b2x[ilayer],
								b2y[ilayer],
								b2z[ilayer]);
		I1[ilayer].setDepth(d[ilayer]);
		I2[ilayer].setDepth(d[ilayer]);

		I1[ilayer].setNu(nu);
		I2[ilayer].setNu(nu);
	}
	for (size_t ilayer = 0; ilayer < nbLayers; ilayer++)
	{
		Qcdx[ilayer]=0;
		Qcdz[ilayer]=0;
		for (size_t jlayer = 0; jlayer < nbLayers; jlayer++)
		{
			Qcdx[ilayer] += -Q0x * (rho1[jlayer] * I1[jlayer].wxx(d[ilayer])
								 +  rho2[jlayer] * I2[jlayer].wxx(d[ilayer]));
			Qcdz[ilayer] += -2.0 * Q0z * (rho1[jlayer] * I1[jlayer].wzz(d[ilayer])
									   +  rho2[jlayer] * I2[jlayer].wzz(d[ilayer]));
		}
		//TODO temporarily the peak center is set to (0, 0) point (only FWHM is interesting)
		Qcx[ilayer] = 0;//* M_PI * Qcdx[ilayer] + dQx[ilayer];
		Qcz[ilayer] = 0;//* M_PI * Qcdz[ilayer] + dQz[ilayer] + Qccz[ilayer];

		//LOG(logINFO)<<"Peak positions:\t"<<ilayer<<"\t"<<Qcx[ilayer] / (2 * M_PI)+Q0x <<"\t"<<Qcz[ilayer] / (2 * M_PI)+Q0z;
		//LOG(logINFO)<<"Peak positions:\t"<<ilayer<<"\t"<<Qcdx[ilayer] * (2 * M_PI) <<"\t"<< Qcdz[ilayer] * (2 * M_PI);
	}
	//we do not take into account here that in fact Q=q+Q0
	Qx=Q0x * 2 * M_PI;
	Qz=Q0z * 2 * M_PI;
}
