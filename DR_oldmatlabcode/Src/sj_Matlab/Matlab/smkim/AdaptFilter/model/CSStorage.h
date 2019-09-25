#ifndef CSStorage_H
#define CSStorage_H

/* $Id: CSStorage.h,v 1.6 2008/06/27 00:21:28 chengs Exp $

   Sen Cheng, Tue Jun  6 14:06:42 PDT 2006
   definition of class CSStorage
*/

#include "../aux/defs.h"
#include "../aux/numerics.h"
#include "CSplinesFct.h"
#include "CSFunction.h"

namespace AFilter {

class CSFunction;
class CSCompress;
class VOp;


/* ************************************************************
                          class CSStorage
   ************************************************************ */
    
class CSStorage {
protected: // data structures
    int startindex, endindex, T;                  // 
    CSFunction *super;
    int nFct;
    int nUpdate, *nPerFct, nTotal;

    // pointers to current values of the spline control points
    double **param;

    mutable lazyVar<double> min, max;
protected: // methods
public:
    CSStorage() : nFct(0), nUpdate(0), nTotal(0), param(0) {};
    virtual ~CSStorage() { if(param) delete[] param; };
    virtual CSStorage* alloc_new() const { hippo_Virtual; }
    virtual CSStorage* alloc_copy() const { hippo_Virtual; }

    virtual const char* getName() const =0;

    virtual void init(int start, int end, CSFunction *f);

    virtual void calcMinMax(const VOp *op) const { hippo_Virtual; }
    
    virtual double getMin(const VOp *op) const { calcMinMax(op);return min; } ;
    virtual double getMax(const VOp *op) const { calcMinMax(op);return max; } ;

    // @@temp. these functions are not necessary
    void setInitial(double x){ hippo_Virtual; };
//    void setInitialEst(int t) const{ hippo_Virtual; };
//    void setInitialToCurrent() { hippo_Virtual; };

    virtual void updateEst(int tStart, int tEnd, double diff[]) {
        hippo_Virtual; 
    }; 

    virtual void allocForConverge(){ hippo_Virtual; };
    virtual void copyForConverge(){ hippo_Virtual; };
    virtual void averageWithCopy() { hippo_Virtual; };
    virtual void maxDiff(double &abs, double &rel) const { hippo_Virtual; };; 

    virtual void getAllParam(int t, double z[]) const { hippo_Virtual; };
    virtual void setAllParam(int t, double z[]) { hippo_Virtual; };
    virtual void setAllParam(int t, double z) { hippo_Virtual; };
    virtual void moveTo(int t) const { hippo_Virtual; }

    virtual void setMexInput(const mxArray *in){ hippo_Virtual; };
    virtual mxArray* getMexOutput() const{ hippo_Virtual; return 0;};

    friend class CSCompress;
}; // class CSStorage

}; //namespace AFilter

#endif   // CSStorage_H
