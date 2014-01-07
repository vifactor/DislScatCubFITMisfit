/*
 * RelaxedInterface.h
 *
 *  Created on: 29 ���. 2011
 *  modified and checked on: 23 nov. 2011
 *      Author: kopp
 */

#ifndef RELAXEDINTERFACE_H_
#define RELAXEDINTERFACE_H_

#include <cmath>

class RelaxedInterface {
public:
	RelaxedInterface();
	virtual ~RelaxedInterface();

	double wxx(double z) const;
	double wzz(double z) const;

	double getDepth() const {return Depth;}

	double Wxxxx(double z) const;
	double Wxzzx(double z) const;
	double Wxzxz(double z) const;
	double Wxxzz(double z) const;
	double Wzxzx(double z) const;
	double Wzzzz(double z) const;

	double Wyzyz(double z) const;

	void preinitialize(double z);

	void setDepth(double d);
	void setNu(double nu)
	{
		Alpha=1.0/(2.0*(1.0-nu));
		Alpha2=Alpha*Alpha;
		AlphaX=(1.0 - 2.0*Alpha)*(1.0 - 2.0*Alpha);
	}
	void setBurgersComponents(double x, double y, double z)
	{
		bx=x; by=y; bz=z;
		bx2=bx*bx;by2=by*by;bz2=bz*bz;
		//LOG(logINFO)<<"burgerComponents:\n"<<bx<<"\t"<<by<<"\t"<<bz;
	}
	double eps;//how close should we approach the interface
protected:
	double Depth;
	double Alpha;
	double bx, by, bz;

	//acceleration variables
	double z2, z3, z4, z5, z6, z7, z8;//powers of z
	double d2, d3, d4, d5, d6, d7, d8;//powers of d
	double z1d7, z2d6, z3d5, z4d4, z5d3, z6d2, z7d1;
	double zpd2, zpd5;//powers of (z+d)
	double Alpha2; //powers of alpha
	double AlphaX;//(1 - 2*Alpha)^2
	double bx2, by2, bz2;
	double denominatorZ3, denominatorD3;
};

#endif /* RELAXEDINTERFACE_H_ */
