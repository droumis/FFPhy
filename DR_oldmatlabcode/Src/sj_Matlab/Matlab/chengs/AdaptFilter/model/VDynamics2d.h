/* $Id: VDynamics2d.h,v 1.7 2008/09/01 21:07:44 chengs Exp $
   File name   : model/VDynamics2d.h
   authors     : Sen Cheng
   created     : Sun 17 Oct 2004 10:49:28 AM PDT
  
 */

#ifndef VDYNAMICS2D_H
#define VDYNAMICS2D_H

#include "../aux/defs.h"
#include "../aux/hippoIO.h"
#include "../aux/Grad.h"
#include "VDynamics.h"

namespace AFilter {

/** Abstract class for 2-dim dynamics.
 * 
 * @author Sen Cheng
 */
class VDynamics2d : public virtual VDynamics
{
private:
protected:
    double *x, *y;
    mutable lazyVar<double> xmin, xmax, ymin, ymax;
public:
    VDynamics2d() {};
    VDynamics2d(VDynamics2d &) {};
    virtual ~VDynamics2d() {};

    virtual int getDim() const { return 2; }

    virtual void init(int _startindex, int _endindex, int _T, 
            double *_x, double *_y, int _nid, traj_type *_id) 
    {
        VDynamics::init(_startindex, _endindex, _T, _nid, _id);
        x= _x; y= _y;
    };

    virtual void setInitialFrom(VDynamics2d *) {hippo_Virtual; };


    virtual double getMinX() const { return xmin; };
    virtual double getMaxX() const { return xmax; };
    virtual double getMinY() const { return ymin; };
    virtual double getMaxY() const { return ymax; };
    //
    // set the boundaries of the domain
    virtual void setMinX(double _xmin) { xmin= _xmin; };
    virtual void setMaxX(double _xmax) { xmax= _xmax; };
    virtual void setMinY(double _ymin) { ymin= _ymin; };
    virtual void setMaxY(double _ymax) { ymax= _ymax; };

    virtual double evalAtTime(int t) const { 
        moveTo(t);
        return eval(x[tCurr], y[tCurr], id[tCurr]); 
    }

    virtual void evalAtTime(int t, Grad &results) const {
        moveTo(t);
        eval(x[tCurr], y[tCurr], id[tCurr], results);
    }

    virtual double evalGradAtTime(int t, double partial[]) const {
        moveTo(t);
        return evalGrad(x[tCurr], y[tCurr], id[tCurr], partial);
    }

    /** calculate f(x) */
    virtual double evalAtTime(int t, double _x, double _y, traj_type _id) const
    {
        moveTo(t);
        return eval(_x, _y, _id);
    }

    virtual double evalNoMove(int t, double partial[]) const {
        return evalGrad(x[t], y[t], id[t], partial);
    }

    virtual void evalNoMove(int t, Grad &results) const {
        eval(x[t], y[t], id[t], results);
    }

    /** calculate f(x) */
    virtual double eval(double x, double y, traj_type id) const =0;

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void eval(double x, double y, traj_type traj, Grad &g) const { hippo_Virtual;}

    /** calculate natural log f(x) */
    virtual double evalLog(double x, double y, traj_type id) const
    { hippo_Virtual; return 0;}

    /** calculate f(x) and df/dz, where z are the parameters */
    virtual double evalGrad(double x,double y, 
            traj_type traj, double partial[]) const =0;

    /** calculate log f(x) and d(log f)/dz, where z are the parameters, log is
     * natural */
    virtual double evalLogGrad(double x,double y, 
            traj_type traj, double partial[]) const
    {hippo_Virtual; return 0;};

    virtual VDynamics2d* link_copy() { hippo_Virtual; return 0; };

//    virtual void setMexInput(const mxArray *in) { VDynamics::setMexInput(in); };
//    virtual mxArray* getMexOutput() { return VDynamics::getMexOutput(); }

}; // class VDynamics2d

}; // namespace AFilter 

#endif   // VDYNAMICS2D_H
