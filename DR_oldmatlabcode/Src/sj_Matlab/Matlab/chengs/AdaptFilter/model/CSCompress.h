#ifndef CSCompress_H
#define CSCompress_H

/* $Id: CSCompress.h,v 1.4 2008/08/24 19:47:03 chengs Exp $

   Sen Cheng, Tue Jun  6 14:06:42 PDT 2006
   definition of class CSCompress
*/

#include "../aux/defs.h"
#include "CSStorage.h"


/* ************************************************************
                          class CSCompress
   ************************************************************ */
namespace AFilter {
    
class CSCompress : public CSStorage {
protected: // data structures

    int outputInterval;

    int nInterval;
    double *tParam;
    bool tParamAllocated;

    mutable int tCurr;
    mutable lazyVar<double> convertMin, convertFac, convertMax;

private: // methods
    void calcConversion() const ;
    inline int getTimestep(int t) const { 
        return (int) floor((t-startindex)/outputInterval); 
    }
    inline void updateParam(int index) const {
        double *tmp= tParam + index*nTotal;
        for (traj_type i= 0; i < nFct; i++) {
            memcpy(param[i], tmp, nPerFct[i]*sizeof(double));
            tmp+=nPerFct[i];
        }
    }
protected:
public:
    CSCompress();
    CSCompress(const CSStorage &);
    virtual ~CSCompress();

    virtual const char* getName() const {return "CSCompress"; };
    void init(int _startindex, int _endindex, CSFunction *f);

    virtual void moveTo(int t) const;
    virtual void getAllParam(int t, double z[]) const;
    virtual void setAllParam(int t, double z[]);
    virtual void setAllParam(int t, double z);

    virtual void calcMinMax(const VOp *op) const;
//    virtual double getMin(const VOp *op) const;
//    virtual double getMax(const VOp *op) const;


    // @@ don't need right now, but maybe later
//    void updateDynamics(int tStart, int tEnd, double diff[]);
    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput() const;

}; // class CSCompress

}; //namespace AFilter

#endif   // CSCompress_H
