/*
 * RelaxedInterface.cpp
 *
 *  Created on: 29 ���. 2011
 *      Author: kopp
 */

#include "RelaxedInterface.h"

RelaxedInterface::RelaxedInterface()
{
	eps=0;
}

RelaxedInterface::~RelaxedInterface()
{
}

double RelaxedInterface::wxx(double z) const
{
	if((z-Depth)<=0)
		return -bx;
	else
		return 0;
}

double RelaxedInterface::wzz(double z) const
{
	if((z-Depth)<=0)
		return bx*(-1 + 2*Alpha);
	else
		return 0;
}

double RelaxedInterface::Wxxxx(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		return denominatorD3*
			(
			bx2*(2.0*d8 + 6.0*z1d7 + 2.0*z3d5*(5.0 - 4.0*Alpha) -
				7.0*z7d1*Alpha2 - 2.0*z8*Alpha2 +
		     2.0*z2d6*(5.0 + 2.0*Alpha) - z6d2*Alpha*(2.0 + 5.0*Alpha) +
		     z5*d3*Alpha*(-8.0 + 9.0*Alpha) +
		     z4d4*(4.0 - 18.0*Alpha + 21.0*Alpha2)) +
		  bz2*(2.0*d8 - 2.0*z1d7 + 2.0*z3d5*(1.0 - 4.0*Alpha) -
		     11.0*z7d1*Alpha2 - 4.0*z8*Alpha2 +
		     2.0*z2d6*(-5.0 + 6.0*Alpha) +
		     z6d2*(4.0 - 10.0*Alpha + 5.0*Alpha2) +
		     5.0*z4d4*(4.0 - 10.0*Alpha + 9.0*Alpha2) +
		     z5d3*(16.0 - 40.0*Alpha + 45.0*Alpha2))
		     );
	}
	//LOG(logINFO)<<"z>0 "<<z-Depth;
	return -denominatorZ3*
			 (
			    bx2*(-2.0*d8*Alpha*(1.0-2.0*Alpha) +
			       z3d5*(-2.0 + 22.0*Alpha - 21.0*Alpha2) +
			       z4d4*(14.0 - 2.0*Alpha - 15.0*Alpha2) +
			       z2d6*(-6.0 + 12.0*Alpha + Alpha2) +
			       4.0*z5d3*(5.0 - 9.0*Alpha + 3.0*Alpha2) +
			       8.0*z6d2*(1.0 - 3.0*Alpha + 3.0*Alpha2) +
			       z1d7*(-2.0 - 2.0*Alpha + 11.0*Alpha2)) +
			    bz2*(2.0*d8*(1.0-2.0*Alpha)*Alpha +
			       z1d7*(2.0 + 2.0*Alpha - 11.0*Alpha2) +
			       5.0*z2d6*(2.0 - 4.0*Alpha + Alpha2) +
			       3.0*z4d4*(2.0 - 10.0*Alpha + 15.0*Alpha2) +
			       z3d5*(14.0 - 50.0*Alpha + 45.0*Alpha2))
			  );
}

double RelaxedInterface::Wxzzx(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		return denominatorD3*
				(
				-bx2*(2.0*d8 + 6.0*z1d7 + z5d3*Alpha2 -
				     7.0*z6d2*Alpha2 - 7.0*z7d1*Alpha2 - 2.0*z8*Alpha2 +
				     2.0*z2d6*(5.0 - 4.0*Alpha2) + 2.0*z3d5*(5.0 + 2.0*Alpha2) +
				     z4d4*(4.0 + 3.0*Alpha2)) +
				  bz2*(-10.0*d8 - 30.0*z1d7 + 11.0*z7d1*Alpha2 + 4.0*z8*Alpha2 +
				     z5d3*(16.0 - 21.0*Alpha2) + z6d2*(4.0 + Alpha2) -
				     5.0*z4d4*(-4.0 + 3.0*Alpha2) + 6.0*z2d6*(-5.0 + 4.0*Alpha2) +
				     2.0*z3d5*(-1.0 + 6.0*Alpha2))
				 );
	}
	return denominatorZ3*
			(
			   bz2*(2.0*d8*Alpha2 +
			      z3d5*(14.0 - 9.0*Alpha2) + 3.0*z4d4*(2.0 - 7.0*Alpha2) +
			      5.0*z2d6*(2.0 + Alpha2) +z1d7*(2.0 + 7.0*Alpha2)) +
			   bx2*(-2.0*d8*Alpha2 + z3d5*(-2.0 + Alpha2) -
			      8.0*z6d2*(-1.0 + Alpha2) + 4.0*z5d3*(5.0 + Alpha2) +
			      z4d4*(14.0 + 3.0*Alpha2) - z1d7*(2.0 + 7.0*Alpha2) -
			      z2d6*(6.0 + 7.0*Alpha2))
			 );
}

double RelaxedInterface::Wxzxz(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		return denominatorD3*
				(
				bx2*(2.0*d8 + z1d7*(6.0 - 4.0*Alpha) -
				     z5d3*(-8.0 + Alpha)*Alpha + 7.0*z7d1*Alpha2 +
				     2.0*z8*Alpha2 + z6d2*Alpha*(2.0 + 7.0*Alpha) +
				     z4d4*(4.0 + 10.0*Alpha - 3.0*Alpha2) +
				     2.0*z3d5*(5.0 + 6.0*Alpha - 2.0*Alpha2) +
				     2.0*z2d6*(5.0 + 2.0*Alpha + 4.0*Alpha2)) +
				bz2*(10.0*d8 + 10.0*z1d7*(3.0 - 2.0*Alpha) + 11.0*z7d1*Alpha2 +
				     4.0*z8*Alpha2 + z5d3*(-16.0 + 8.0*Alpha - 21.0*Alpha2) +
				     z6d2*(-4.0 + 2.0*Alpha + Alpha2) -
				     5.0*z4d4*(4.0 - 2.0*Alpha + 3.0*Alpha2) +
				     2.0*z3d5*(1.0 - 2.0*Alpha + 6.0*Alpha2) +
				     2.0*z2d6*(15.0 - 14.0*Alpha + 12.0*Alpha2))
				 );
	}

	return -denominatorZ3*
			  (
			    bx2*(2.0*d8*(-1.0 + Alpha)*Alpha +
			       z4d4*(14.0 + 26.0*Alpha - 3.0*Alpha2) +
			       4.0*z5d3*(5.0 + Alpha - Alpha2) -
			       z3d5*(2.0 - 18.0*Alpha + Alpha2) +
			       8.0*z6d2*(1.0 - Alpha + Alpha2) +
			       z2d6*(-6.0 + 7.0*Alpha2) +
			       z1d7*(-2.0 - 6.0*Alpha + 7.0*Alpha2)) +
			    bz2*(-2.0*d8*(-1.0 + Alpha)*Alpha +
			       z1d7*(2.0 + 6.0*Alpha - 7.0*Alpha2) -
			       5.0*z2d6*(-2.0 + Alpha2) +
			       3.0*z4d4*(2.0 - 6.0*Alpha + 7.0*Alpha2) +
			       z3d5*(14.0 - 22.0*Alpha + 9.0*Alpha2))
			   );
}

double RelaxedInterface::Wxxzz(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		//TODO CHECK
		return denominatorD3*
				(
				 bx2*(d8*(2.0 - 4.0*Alpha) + 6.0*z1d7*(1.0 - 2.0*Alpha) -
				     z5d3*Alpha2 + 7.0*z6d2*Alpha2 + 7.0*z7d1*Alpha2 +
				     2.0*z8*Alpha2 + z4d4*(4.0 - 8.0*Alpha - 3.0*Alpha2) -
				     2.0*z2d6*(-5.0 + 10.0*Alpha + 2.0*Alpha2) +
				     2.0*z3d5*(5.0 - 10.0*Alpha + 4.0*Alpha2)) +
				 bz2*(d8*(2.0 - 4.0*Alpha) + 11.0*z7d1*Alpha2 +
				     4.0*z8*Alpha2 + 2.0*z1d7*(-1.0 + 2.0*Alpha) +
				     z5d3*(16.0 - 32.0*Alpha - 5.0*Alpha2) +
				     5.0*z4d4*(4.0 - 8.0*Alpha + Alpha2) +
				     2.0*z3d5*(1.0 - 2.0*Alpha + 4.0*Alpha2) +
				     z6d2*(4.0 - 8.0*Alpha + 5.0*Alpha2) -
				     2.0*z2d6*(5.0 - 10.0*Alpha + 6.0*Alpha2))
				  );
	}

	return denominatorZ3*
			(
			  -bz2*(2.0*d8*Alpha2 +
			      5.0*z2d6*(2.0 - 4.0*Alpha + 3.0*Alpha2) +
			      z3d5*(14.0 - 28.0*Alpha + 5.0*Alpha2) -
			      3.0*z4d4*(-2.0 + 4.0*Alpha + 5.0*Alpha2) +
			      z1d7*(2.0 - 4.0*Alpha + 9.0*Alpha2)) +
			   bx2*(2.0*d8*Alpha2 +
				  8.0*z6d2*(-1.0 + 2.0*Alpha) +
			      z4d4*(-14.0 + 28.0*Alpha - 17.0*Alpha2) +
			      z3d5*(2.0 - 4.0*Alpha + Alpha2) -
			      4.0*z5d3*(5.0 - 10.0*Alpha + 6.0*Alpha2) +
			      z1d7*(2.0 - 4.0*Alpha + 9.0*Alpha2) +
			      z2d6*(6.0 - 12.0*Alpha + 13.0*Alpha2))
			 );
}

double RelaxedInterface::Wzxzx(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
	return denominatorD3*
			(
			 bx2*(2.0*d8 + 7.0*z7d1*Alpha2 + 2.0*z8*Alpha2 -
			     z5d3*Alpha*(8.0 + Alpha) + 2.0*z1d7*(3.0 + 2.0*Alpha) +
			     z6d2*Alpha*(-2.0 + 7.0*Alpha) +
			     z4d4*(4.0 - 10.0*Alpha - 3.0*Alpha2) -
			     2.0*z3d5*(-5.0 + 6.0*Alpha + 2.0*Alpha2) +
			     2.0*z2d6*(5.0 - 2.0*Alpha + 4.0*Alpha2)) +
			  bz2*(10.0*d8 + 11.0*z7d1*Alpha2 + 4.0*z8*Alpha2 +
			     10.0*z1d7*(3.0 + 2.0*Alpha) +
			     z6d2*(-4.0 - 2.0*Alpha + Alpha2) -
			     5.0*z4d4*(4.0 + 2.0*Alpha + 3.0*Alpha2) +
			     2.0*z3d5*(1.0 + 2.0*Alpha + 6.0*Alpha2) +
			     2.0*z2d6*(15.0 + 14.0*Alpha + 12.0*Alpha2) -
			     z5d3*(16.0 + 8.0*Alpha + 21.0*Alpha2))
			  );
	}

	return -denominatorZ3*
			 (
			    bx2*(2.0*d8*Alpha*(1.0 + Alpha) +
			       z4d4*(14.0 - 26.0*Alpha - 3.0*Alpha2) -
			       4.0*z5d3*(-5.0 + Alpha + Alpha2) +
			       8.0*z6d2*(1.0 + Alpha + Alpha2) -
			       z3d5*(2.0 + 18.0*Alpha + Alpha2) +
			       z2d6*(-6.0 + 7.0*Alpha2) +
			       z1d7*(-2.0 + 6.0*Alpha + 7.0*Alpha2)) +
			    bz2*(-2.0*d8*Alpha*(1.0 + Alpha) +
			       z1d7*(2 - 6.0*Alpha - 7.0*Alpha2) -
			       5.0*z2d6*(-2.0 + Alpha2) +
			       3.0*z4d4*(2.0 + 6.0*Alpha + 7.0*Alpha2) +
			       z3d5*(14.0 + 22.0*Alpha + 9.0*Alpha2))
			   );
}

double RelaxedInterface::Wzzzz(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		return denominatorD3*
			(
			 bx2*(2.0*d8*AlphaX +
			     6.0*z1d7*AlphaX + z6d2*(2.0 - 9.0*Alpha)*Alpha +
			     z5d3*(8.0 - 7.0*Alpha)*Alpha - 7.0*z7d1*Alpha2 -
			     2.0*z8*Alpha2 + z4d4*(4.0 + 2.0*Alpha + Alpha2) +
			     2.0*z3d5*(5.0 - 16.0*Alpha + 12.0*Alpha2) +
			     2.0*z2d6*(5.0 - 22.0*Alpha + 24.0*Alpha2)) +
			  bz2*(2.0*d8*AlphaX - 2.0*z1d7*AlphaX -
			     11.0*z7d1*Alpha2 - 4.0*z8*Alpha2 +
			     2.0*z3d5*(1.0 - 4.0*Alpha2) +
			     z6d2*(4.0 - 6.0*Alpha + Alpha2) +
			     5.0*z4d4*(4.0 - 6.0*Alpha + 5.0*Alpha2) -
			     2.0*z2d6*(5.0 - 14.0*Alpha + 8.0*Alpha2) +
			     z5d3*(16.0 - 24.0*Alpha + 29.0*Alpha2))
			  );
	}
	return -denominatorZ3*
			 (
			    bz2*(-2.0*d8*Alpha +
			       z1d7*(2.0 - 10.0*Alpha + Alpha2) +
			       z3d5*(14.0 - 6.0*Alpha + Alpha2) +
			       5.0*z2d6*(2.0 - 4.0*Alpha + Alpha2) +
			       3.0*z4d4*(2.0 + 2.0*Alpha + 3.0*Alpha2)) +
			    bx2*(2.0*d8*Alpha -
			    	z1d7*(2.0 - 10.0*Alpha + Alpha2) +
			       8.0*z6d2*(1.0 - Alpha + Alpha2) +
			       z2d6*(-6.0 + 12.0*Alpha + Alpha2) +
			       4.0*z5d3*(5.0 - 11.0*Alpha + 5.0*Alpha2) +
			       z3d5*(-2.0 - 14.0*Alpha + 15.0*Alpha2) +
			       z4d4*(14.0 - 54.0*Alpha + 37.0*Alpha2))
			   );
}

double RelaxedInterface::Wyzyz(double z) const
{
	if(z==Depth)
	{
		z=Depth+eps;
	}
	if((z-Depth)<0)
	{
		return by2*z2/(4.0*M_PI*(d3 - Depth*z2));
	}

	return by2*d2/(4.0*M_PI*(z3 - d2*z));
}

void RelaxedInterface::preinitialize(double z)
{
	z2=pow(z, 2.0);
	z3=pow(z, 3.0);
	z4=pow(z, 4.0);
	z5=pow(z, 5.0);
	z6=pow(z, 6.0);
	z7=pow(z, 7.0);
	z8=pow(z, 8.0);

	zpd2=pow(z+Depth, 2.0);
	zpd5=pow(z+Depth, 5.0);

	z1d7=z*d7;
	z2d6=z2*d6;
	z3d5=z3*d5;
	z4d4=z4*d4;
	z5d3=z5*d3;
	z6d2=z6*d2;
	z7d1=z7*Depth;

	denominatorD3=1.0/(8.*M_PI*(Depth- z)*d3*zpd5);
	denominatorZ3=1.0/(8.*M_PI*(Depth- z)*z3*zpd5);
}

void RelaxedInterface::setDepth(double d)
{
	Depth=d;
	d2=pow(d, 2.0);
	d3=pow(d, 3.0);
	d4=pow(d, 4.0);
	d5=pow(d, 5.0);
	d6=pow(d, 6.0);
	d7=pow(d, 7.0);
	d8=pow(d, 8.0);
}
