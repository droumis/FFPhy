/* $Id: AdaptModel.h,v 1.13 2008/09/01 18:24:22 chengs Exp $
   File name   : model/AdaptModel.h
   authors     : Sen Cheng
   created     : Sun 17 Oct 2004 10:51:14 AM PDT
  
 */
#ifndef ADAPTMODEL_H
#define ADAPTMODEL_H

#include "../aux/defs.h"
#include "../aux/TData.h"
#include "../aux/hippoIO.h"
#include "../aux/Grad.h" 

namespace AFilter {

/** Abstract class for models to be used in adaptive filtering algorithms.
 * 
 * @author Sen Cheng
 */
class AdaptModel 
{
private:
protected:
    int startindex;
    int endindex;
    TData *data;

    bool convCopy;
    bool real; // =true for parameters, =false for variances
private:
    virtual void calcRescaled (double *&rescaled, int &n) const;
protected:

    AdaptModel(TData *_data=0) : startindex(0), endindex(0), data(_data), convCopy(false), real(true) { }
    virtual AdaptModel* alloc() const =0;

public:
    virtual ~AdaptModel() {};

    virtual const char* getName() const =0;

    virtual const TData* getData() const { return data; }

    virtual void init() =0;
    virtual void setUnreal() { real= false; };
    virtual void moveTo(int timeindex) const =0;
//    virtual void update(int t) {hippo_Virtual; };

    /* */
    virtual double eval(int timeindex) const =0;

//    virtual double evalAtTime(double time) const { 
//        return eval(data->timeToIndex(time));
//    }

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void eval(int timeindex, Grad &g) const =0;

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void evalNoMove(int timeindex, Grad &g) const =0;

    /** calculate f(x) and df/dz, where z are the parameters
    */
    virtual double evalGrad(int timeindex, double partial[]) const =0;

    /// For filtering algorithm. Evaluate without moving.
    virtual double evalNoMove(int timeindex, double partial[]) const =0;

//    /** calculate log f(x) and d(log f)/dz, where z are the parameters
//    */
//    virtual double evalLogGrad(int timeindex, double partial[]) const
//    {hippo_Virtual; return 0; };

    virtual int getStartindex() const { return startindex;};
    virtual int getEndindex() const { return endindex;};

    virtual int getNParam() const {hippo_Virtual; return 0;};
    virtual int getNAllParam() const {hippo_Virtual; return 0;};

    virtual bool getParam(int t, double z[], int a[]) const =0;
    virtual void getAllParam(int t, double x[]) const =0;

    virtual void setParam(int t, double x[]) { hippo_Virtual; }
    virtual void setAllParam(int t, double x) { hippo_Virtual; }
    virtual void setAllParam(int t, double x[]) { hippo_Virtual; }

    virtual void setInitialEst(int t) =0;
 
    virtual void updateModel(int tStart, int tEnd, double diff[]) =0;

    //     virtual int getNChangeParam() { hippo_Virtual; return 0; };
    //     virtual void copyTo(double *) { hippo_Virtual; }
    //     virtual void getMaxDiff(double *, double &maxdiff, double &reldiff, double cut) { 
    // 	hippo_Virtual; }

    virtual void allocForConvergence() { hippo_Virtual; }
    virtual void copyForConverge() { hippo_Virtual; }
    virtual bool hasConverged() { hippo_Virtual; return false; }
    virtual void averageWithCopy() { hippo_Virtual; }

    virtual AdaptModel* link_copy() { hippo_Virtual; return 0; };
//    virtual AdaptModel* link_copy() = 0;

    // functions for M-step of EM algorithm
    virtual void reestimate(AdaptModel *var) {hippo_Virtual; };
    virtual double evalCore(int t) const { hippo_Virtual; return 0;}
    virtual double evalCoreVar(int t) const { hippo_Virtual; return 0;}

    virtual void setMexInput(const mxArray *in) =0;
    virtual mxArray* getMexOutput()= 0;

/**** for generating spike trains ******/

    /* */
    virtual double evalAtIsi(int timeindex, double isi) {hippo_Virtual; return 0;};

}; // class AdaptModel

}; // namespace AFilter

#endif   // ADAPTMODEL_H
