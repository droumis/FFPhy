/* $Id: CSplinesFct.h,v 1.4 2006/07/04 15:30:19 chengs Exp $
   File name   : model/CSplinesFct.h
   authors     : Sen Cheng
   created     : 2004/05/12
  
   Data structures and functions to manipulate splines used in adaptive
   filtering. 
 */

#ifndef CSPLINESFCT_H
#define CSPLINESFCT_H

#include <iostream>
#include "../aux/hippoIO.h"

using namespace std;

namespace AFilter {

/* ************************************************************
                          class  CSFct
   ************************************************************ */

class CSStorage;

class CSFct {
protected:
    double *param;

    bool allocated;
    int nParam;

protected:
    virtual void allocParam(int n);
    virtual void dealloc();
//    double* getParamPointer() const { return param; }
public:
    CSFct();
    virtual ~CSFct();

    inline void getParam(double p[]) { memcpy(p, param, nParam*sizeof(double)); }
    inline void setParam(double p[]) { memcpy(param, p, nParam*sizeof(double)); }
    inline void setParam(double p) { for(int i=0; i<nParam; i++) param[i]= p; }
    inline void getParamPartial(double z[], int n, int a[]) {
        for (int i= 0; i < n; i++) {z[i]= param[a[i]];
//            hippo_Assert(a[i]<nParam, "index out of bounds");
        }
    }
    inline void setParamPartial(double z[], int n, int a[]) {
        for (int i= 0; i < n; i++) param[a[i]]= z[i];
    }
    inline void addParamPartial(double z[], int n, int a[]) {
        for (int i= 0; i < n; i++) param[a[i]]+= z[i];
    }
    inline void subtractParamPartial(double z[], int n, int a[]) {
        for (int i= 0; i < n; i++) param[a[i]]-= z[i];
    }


    friend class CSStorage;
};

/* ************************************************************
                          class  CSplinesFct
   ************************************************************ */


/** 1-dim cubic cardinal splines function.
 * 
 * @author Sen Cheng
 * @version 1.1
 */
class CSplinesFct : public CSFct {
protected:
    int      ncpx;	      // the numbers of spatial control points
    double   *cpx;	   // spatial control points in x
public:
    CSplinesFct();
    CSplinesFct(int _ncpx, double* _cpx);
    virtual ~CSplinesFct();
    
    virtual void init(int _ncpx, double* _cpx);

    virtual double getMinX() const { return cpx[1]; }
    virtual double getMaxX() const { return cpx[ncpx-2]; }

    virtual void getPartialIndices(int a[], int csegx);

    virtual int findSegment(double x) const;

    virtual double eval(double x, int csegx= -1)const ;
    virtual double evalPartial(double x, double partialx[], int csegx=-1)const ;
    virtual void release();

    virtual double integral() const;
};

/* ************************************************************
                          class  CSplinesFct2d
   ************************************************************ */

class CSplinesFct2d : public CSFct {
protected:
    double            *cpx;	   // spatial control points in x
    double            *cpy;	   // spatial control points in y
    
    int            ncpx;	      // the numbers of spatial control points
    int            ncpy;	      // number of control points for theta phasse
protected:
    virtual int getLinear(int ix, int iy) const { return ix*ncpy+iy; }
public:
    CSplinesFct2d();
    CSplinesFct2d(int _ncpx, int _ncpy, double* _cpx, double* _cpy);
    virtual ~CSplinesFct2d();
    
    virtual void init(int _ncpx, int _ncpy, double* _cpx, double* _cpy);

    virtual void getPartialIndices(int a[], int csegx, int csegy);

    virtual double getMinX() const { return cpx[1]; }
    virtual double getMaxX() const { return cpx[ncpx-2]; }
    virtual double getMinY() const { return cpy[1]; }
    virtual double getMaxY() const { return cpy[ncpy-2]; }

    virtual void findSegments(double x, double y, int &csegx, int &csegy) const;

    // return value of function at (x,y)
    virtual double eval(double x, double y, int csegx=-1, int csegy=-1)const ;
    
    // calc. value and partial derivatives of function at (x,y)
    virtual double evalPartial(double x, double y, double partialx[], 
            double partialy[], int csegx=-1, int csegy=-1)const ;

    virtual void release();

//    friend ostream& operator << (ostream &, CSplinesFct2d &);
};

/* ************************************************************
                          class  CSplinesFctLinCirc
   ************************************************************ */

class CSplinesFctLinCirc : public CSplinesFct2d {
protected:
    virtual int getLinear(int ix, int iy) const { return ix*ncpy+mod(iy,ncpy-1); }
public:
    
    CSplinesFctLinCirc() {};
//    CSplinesFctLinCirc(int _ncpx, int _ncpy, double* _cpx, double* _cpy);
    virtual ~CSplinesFctLinCirc() {};
    
    virtual double getMinY() const { return cpy[0]; }
    virtual double getMaxY() const { return cpy[ncpy-1]; }

    virtual void findSegments(double x, double y, int &csegx, int &csegy) const;

}; // class CSplinesFctLinCirc

}; // namespace AFilter

#endif // CSLPLINESFCT_H
