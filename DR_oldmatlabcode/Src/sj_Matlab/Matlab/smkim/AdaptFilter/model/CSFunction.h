/* $Id: CSFunction.h,v 1.5 2008/08/24 19:47:03 chengs Exp $
   Sen Cheng,  Tue Jun  6 15:24:02 PDT 2006
 */

#ifndef CSFunction_H
#define CSFunction_H

#include <vector>
#include "../aux/defs.h"
#include "../aux/numerics.h"
#include "../aux/hippoIO.h"
#include "VDynamics.h"
#include "CSplinesFct.h"
#include "Operators.h"
//#include "CSStorage.h"

namespace AFilter {

class CSStorage;

/* ************************************************************
                          class  CSFunction
   ************************************************************ */


/**
 */
class CSFunction : public virtual VDynamics {

protected:
    int nUpdate, *nPerFct, nTotal;

    // current time and function
    vector<CSFct *> fct;
    VOp *op;
    bool allocOp;

    CSStorage *store;
    bool inputCompressed;
    bool outputCompress;
    int outputInterval;

    bool rectify;

    // set all control points to this value initially
    lazyVar<double> initAll;

    double conv, convper;

    // variable for mapping control points across different trajectories
    typedef struct {
        traj_type t1, t2;
        double min, max;
    } mapItem;
    bool mapping;
    vector<mapItem> map;

protected:
    void copyInit(CSFunction *v) const;
public:
    
    CSFunction();
    virtual ~CSFunction();

    virtual CSFunction* alloc_new() const =0;

    void init(int _startindex, int _endindex, int _T, int _nid, traj_type *_id);
    
    inline int* getNPerFct() const { return nPerFct; }
    inline int getNUpdate() const { return nUpdate; }
    inline int getNTotal() const { return nTotal; }

    //@@ for comp with old code, remove
    virtual int getNParam() const { return nUpdate; };
    virtual int getNAllParam() const { return  nTotal; }; 

    virtual bool getParam(int t, double x[], int a[]) const;

    virtual void getAllParam(int t, double z[]) const;
    virtual void setAllParam(int t, double z[]);
    virtual void setAllParam(int t, double z);
    virtual void setValues(int t, double z);

    virtual double getMinZ() const;
    virtual double getMaxZ() const;

    virtual int getOutputInterval() const { return outputInterval; }

    virtual void getActiveIndices(int t, int &ntraj, traj_type trajs[], 
            int ind[]) const =0;

    virtual void updateDynamics(int tStart, int tEnd, double diff[]); 
    
    virtual void moveTo(int timeindex) const;

    virtual const vector<CSFct *>* getFunctionObjs() const { return &fct;};

    virtual void allocForConverge();
    virtual void copyForConverge();
    virtual bool hasConverged();
    virtual void averageWithCopy();

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
//    virtual void addMexOutput(mxArray* );

}; //class CSFunction {

}; // namespace AFilter {

#endif // CSFunction_H
