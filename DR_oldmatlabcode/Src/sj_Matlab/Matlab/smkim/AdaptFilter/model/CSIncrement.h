#ifndef CSIncrement_H
#define CSIncrement_H

/* $Id: CSIncrement.h,v 1.5 2008/06/27 00:21:28 chengs Exp $

   Sen Cheng, Tue Jun  6 14:06:42 PDT 2006
   definition of class CSIncrement
*/

#include "../aux/defs.h"
#include "CSStorage.h"

namespace AFilter {

/* ************************************************************
                          class CSIncrement
   ************************************************************ */
    
class CSIncrement : public CSStorage {
protected: // data structures
    mutable int tIni, tCurr;
    mutable traj_type traj;
    double *initial, *updates;

    struct {
        bool initial, updates;

        bool inds, trajs;
    } allocated;

    mutable int nt;
    traj_type *trajs;   // new
    int *inds;         // new

    CSIncrement *convCopy;
private: // methods
protected: // methods
    void setInitialFromCurrent();
    void setToInitial();
    void resetInitialTo(int t);
public:
    CSIncrement();
    virtual ~CSIncrement();
    virtual CSIncrement* alloc_new() const { return new CSIncrement; }
    virtual CSIncrement* alloc_copy() const;

    virtual const char* getName() const {return "CSIncrement"; };
    virtual void init(int _startindex, int _endindex, CSFunction *);

    // @@temp. these functions are not necessary
//    void setInitial(double x);
//    void setInitialEst(int t) const;
//

    virtual void calcMinMax(const VOp *op) const;
    virtual void moveTo(int t) const;

    virtual void getAllParam(int t, double z[]) const;
    virtual void setAllParam(int t, double z[]);
    virtual void setAllParam(int t, double z);

//    virtual double getMin(const VOp *op) const;
//    virtual double getMax(const VOp *op) const;

    virtual void allocForConverge() { convCopy= alloc_copy(); }
    virtual void copyForConverge();
    virtual void averageWithCopy();
    virtual void maxDiff(double &abs, double &rel) const;

    void updateEst(int tStart, int tEnd, double diff[]);


    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput() const;

}; // class CSIncrement

}; //namespace AFilter

#endif   // CSIncrement_H
