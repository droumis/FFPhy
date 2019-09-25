/* $Id: numerics.h,v 1.4 2007/01/10 17:04:38 chengs Exp $
   authors     : Sen Cheng
   created     : 2005/02/22
 */

#ifndef AF_numerics_H
#define AF_numerics_H

#include "hippoIO.h"

#include <stddef.h>
#include <assert.h>
#include <math.h>

namespace AFilter {

const double tiny= 1e-20;

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

inline double sq(const double x) { return x*x; }

#define MOD(x,y) ((x)%(y) >= 0 ? (x)%(y) : (x)%(y) + (y))

inline int mod(int x, int y) { int m= x%y; return m >= 0 ? m : m+y; }


inline double angdist(double x, double y)
{
    double u= fmod(y-x,2*M_PI);
    if(u > M_PI) u-= 2*M_PI;
    if(u < -M_PI) u+= 2*M_PI;
    return u;
}

inline void circstats(int n, double x[], double p[], double &m, double &r)
{
    double u= 0, v= 0, norm= 0;
    for (int j= 0; j < n; j++) {
        u+= p[j]*cos(x[j]);
        v+= p[j]*sin(x[j]);     
        norm+= p[j];
    }
    u/= norm; v/= norm;
    m= atan2(v,u);
    if(m<0) m+= 2*M_PI;
    r= 1-sqrt((u*u)+(v*v));
}


template <typename T>
struct lazyVar
{
	private:
		T _v;
	public:
		bool valid;
		inline lazyVar() {valid=false;}
		/* handy operators */
		inline operator T() { assert(valid); return _v; } 
		inline operator T() const { assert(valid); return _v; } 
		inline T& operator=(const T& v) { valid=true; return _v=v; }
		inline T& operator+=(const T& v) { assert(valid); return _v+=v; }
		inline T& operator-=(const T& v) { assert(valid); return _v-=v; }
};

template <typename T> struct vec
{
    T *x;
    bool valid;
    bool alloc;
    size_t n;

    inline vec() { alloc= valid= false; n= 0; }
    ~vec() { dealloc(); }

    void allocate (size_t _n) { 
        if(alloc && _n==n) return;
        dealloc();
        x= new T[_n]; 
        n=_n; alloc= valid= true; 
    }
    inline void dealloc() { 
        if(alloc && x) delete[] x; alloc= valid= false; n= 0; 
    }
    size_t getN() const { return n; }
    inline void copy(const vec<T>& v) { 
        assert(valid && v.valid && n== v.n);
        for(size_t i=0; i<n; ++i) x[i]= v[i];
    }

    inline T& operator[] (size_t i) { assert(valid && i >= 0 && i < n); return x[i]; } 
    inline T& operator[] (size_t i) const { assert(valid && i >= 0 && i < n); return x[i]; } ;
    inline bool operator== (T* _x) { assert(valid); return x== _x; }
    inline T* operator= (T* _x) { dealloc(); valid=true; return x= _x; }
    inline vec<T>& operator= (const vec<T>& v) { 
        if(v.valid)
            if(v.alloc) {
                allocate(v.n);
                copy(v);
            } else (*this)= v.x;
        else dealloc();
        return *this;
    }

}; // struct vec


} // namespace AFilter

#endif   // AF_numerics_H
