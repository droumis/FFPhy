#ifndef VDYNAMICS1D_H
#define VDYNAMICS1D_H

/* $Id: VDynamics1d.h,v 1.6 2008/09/01 18:24:22 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VDynamics.h"
#include "../aux/hippoIO.h"
#include "../aux/Grad.h"

namespace AFilter {

class VDynamics1d : public virtual VDynamics
{
private:
protected:
    double *x;
    mutable lazyVar<double> xmin, xmax;
public:
    VDynamics1d() {};
    VDynamics1d(VDynamics1d &) {};
    virtual ~VDynamics1d() {};
    
    virtual int getDim() const { return 1; }

    virtual void init(int _startindex, int _endindex, int _T, 
            double *_x, int _nid, traj_type *_id) 
    {
        VDynamics::init(_startindex, _endindex, _T, _nid, _id);
        x= _x;
    };

    virtual void setInitialFrom(VDynamics1d *) {hippo_Virtual; };


    virtual double getMinX() const { return xmin; };
    virtual double getMaxX() const { return xmax; };

    virtual void setMinX(double _xmin) { xmin= _xmin; };
    virtual void setMaxX(double _xmax) { xmax= _xmax; };

    virtual double evalAtTime(int t) const { 
        moveTo(t);
        return eval(x[tCurr], id[tCurr]); 
    }

    virtual void evalAtTime(int t, Grad &results) const {
        moveTo(t);
        eval(x[tCurr], id[tCurr], results);
    }

    virtual double evalGradAtTime(int t, double partial[]) const {
        moveTo(t);
        return evalGrad(x[tCurr], id[tCurr], partial);
    }

    /** calculate f(x) */
    virtual double evalAtTime(int t, double _x, traj_type id) const {
        moveTo(t);
        return eval(_x, id);
    }

    virtual double evalNoMove(int t, double partial[]) const {
        return evalGrad(x[t], id[t], partial);
    }

    virtual void evalNoMove(int t, Grad &results) const {
        eval(x[t], id[t], results);
    }

    /** calculate f(x) */
    virtual double eval(double x, traj_type id) const =0;

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void eval(double x, traj_type traj, Grad &g) const =0;

    /** calculate f(x) and df/dz, where z are the parameters */
    virtual double evalGrad(double x, traj_type traj, double partial[]) const =0;

//    /** calculate natural log f(x) */
//    virtual double evalLog(double x, traj_type id) const =0;

//    /** calculate log f(x) and d(log f)/dz, where z are the parameters, log is
//     * natural */
//    virtual double evalLogGrad(double x, traj_type traj, double partial[]) const
//    {hippo_Virtual; return 0;};

    virtual VDynamics1d* link_copy() { hippo_Virtual; return 0; };

//    virtual void setMexInput(const mxArray *in) { VDynamics::setMexInput(in); };
//    virtual mxArray* getMexOutput() { return VDynamics::getMexOutput(); }
}; // class VDynamics1d

}; // namespace AFilter 
#endif   // VDYNAMICS1D_H
