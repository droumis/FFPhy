#ifndef CSPLINES_H
#define CSPLINES_H

/* $Id: CSplines.h,v 1.8 2008/09/01 18:24:22 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description

   need to call init before use
*/

#include "VDynamics1d.h"
#include "../aux/defs.h"
//#include "CSplinesFct.h"
//#include "CSFunction.h"
#include "CSStorage.h"


/* ************************************************************
                          class CSplines
   ************************************************************ */
namespace AFilter {


/** 1-d cardinal splines model
 * @author: Sen Cheng
 * @date:   May 2004
 */
class CSplines : public VDynamics1d, public CSFunction {
private:
    mutable lazyVar<double> xmin, xmax;
protected:
    /* at each timestep, store only 4 control points that change
       0:xval1   1:xval2   2:xval3   3:xval4 
    */

    vector<CSplinesFct *> fctCurr;
    int            *ncpx;	      // the numbers of spatial control points
    double         **cpx;  	   // the list of temporal control points

    int *csegx;

    struct {
        bool controls;
    } allocated;

protected:
//    void alloc(int NTraj, int NSteps);
    virtual void dealloc();
    virtual void initObj();

//    virtual void updateFromTheta(int time, traj_type traj, int tsegx) const;
    virtual void zero_out();
//    virtual void addToParam(traj_type traj, double p[], int segx= -1) const;
//    virtual void setParam_p(traj_type traj, double p[], int segx= -1) const;
//    virtual void getParam_p(traj_type traj, double p[], int segx= -1) const;
    void copyInit(CSplines *v) const;
public:
    CSplines();
    CSplines(CSplines &);
    virtual ~CSplines();

    virtual const char* getName() const {return "CSplines"; };
    virtual CSplines* alloc_new() const { return new CSplines; }
    virtual CSplines* link_copy();

    virtual void init(int startindex, int endindex, int _T, double *x, int nid, traj_type *id);

    virtual double getMinX() const;
    virtual double getMaxX() const;

//    virtual void getTheta(int t, traj_type traj, double out[]) const;
//    virtual traj_type getParam(int t, double [], int a[]) const;
//    virtual void getAllParam(int t, double z[]) const;

    virtual void getActiveIndices(int t, int &ntraj, traj_type trajs[], 
            int ind[]) const;


//    virtual CSFct* getFunctionObjs() const { return fctCurr; }

//    virtual void setInitial(double);
//    virtual void setInitialToCurrent();
//    virtual void setInitialEst(int t= -1) const;

//    virtual void moveTo(int timeindex) const;
//    virtual double eval(int t, double x, traj_type id)const ;

    virtual double eval(double x, traj_type id)const ;

    virtual void  eval(double x, traj_type id, Grad &g)const ;

    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(double x, traj_type id, double partial[])const ;

//    virtual void updateDynamics(int tStart, int tEnd, double diff[]);
    virtual double getIntegral(int timestep, traj_type traj) const;
    //
//    virtual bool hasConverged(VDynamics *);
//    virtual CSplines* allocForConverge();
//    virtual void copyForConverge(VDynamics *);


    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();


//    friend class CardSplines;
};

} // namespace AFilters

#endif   // CSPLINES_H
