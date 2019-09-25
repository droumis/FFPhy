// $Id: numUtils.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __NUMUTILS_H__
#define __NUMUTILS_H__

#include <math.h>

// TODO add:
// inline double getRes(double xMin, double xMax, int nPoints)
// inline double interpx(int i, int ni, double min, double max)
// gets the smallest multiple of xRes greater, less respectively than x
inline double getCeilingMultiple(double x, double xInc) 
{
	int q = (int)(x/xInc);
	double _x = q*xInc;
	if (fabs(x-_x) < 0.01f)
		return x;
	return _x+xInc;
}
inline double getFloorMultiple(double x, double xInc)
{ 
	int q = (int)(x/xInc);
	double _x = q*xInc;
	if (fabs(x-_x) < 0.01f)
		return x;
	return _x-xInc;
}

inline double normalize(double x, double xMin, double xMax)
{ return (x-xMin)/(xMax-xMin); }

inline double rerange(double x, double xMin, double xMax, double rangeMin, double rangeMax)
{ return rangeMin + (rangeMax-rangeMin)*normalize(x, xMin, xMax); }

#endif
