#ifndef CONST_H
#define CONST_H

/* $Id: Const.h,v 1.12 2008/09/01 21:07:44 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../model/VDynamics1d.h"
#include "../model/VDynamics2d.h"
#include "../model/AdaptModel.h"
//#include "../aux/numerics.h"
#include "Operators.h"

/* ************************************************************
                          class Const
   ************************************************************ */

namespace AFilter {


class Const : public VDynamics1d, public VDynamics2d, public AdaptModel
{
private:
    double constant;
    double xsav;
    double conv, convper;
private:
    void zero_out();
protected:
    AdaptModel* alloc() const { return new Const; }
public:
    Const(Const &);
    Const(TData *data= 0);
    virtual ~Const();
    
    const char* getName() const { return "Const"; }

    virtual int getDim() const { return 0; }

    void init();
    void init(int _startindex, int _endindex, int _T, double *_x, 
                int _nFct, traj_type *_id);
    void init(int _startindex, int _endindex, int _T, double *_x, double *y, 
                int _nFct, traj_type *_id);

    virtual void moveTo(int t) const { tCurr= t; };

/*** begin VDynamics ***/

    virtual double evalAtTime(int t) const  { 
        moveTo(t);
        return id[t]<0? 1: constant; 
    }

    virtual void evalAtTime(int t, Grad &results) const  {
        moveTo(t);
        results.f= (id[t]<0) ? 1 : constant;
    }

    virtual double evalGradAtTime(int t, double partial[]) const {
        moveTo(t);
        return id[t]<0? 1: constant; 
    }

    virtual double evalNoMove(int t, double partial[]) const {
        return id[t]<0? 1: constant; 
    }

    virtual void evalNoMove(int t, Grad &results) const {
        results.f= (id[t]<0) ? 1 : constant;
    }

/*** end VDynamics ***/

/*** begin VDynamics1d ***/

    /** calculate f(x) */
    virtual double eval(int t, double _x, traj_type traj) const { 
        return traj<0? 1: constant; 
    }

    /** calculate f(x) */
    virtual double eval(double _x, traj_type traj) const { 
        return traj<0? 1: constant; 
    }

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    inline void eval(double x, traj_type traj, Grad &g) const {
        g.f= (traj<0) ? 1 : constant;
    }

    /** calculate f(x) and df/dz, where z are the parameters */
    inline double evalGrad(double x, traj_type traj, double partial[]) const {
        return (traj<0) ? 1 : constant;
    }
/*** end VDynamics1d ***/

/*** begin VDynamics2d ***/

    /** calculate f(x) */
    virtual double eval(int t, double _x, double _y, traj_type traj) const { 
        return traj<0? 1: constant; 
    }

    /** calculate f(x) */
    virtual double eval(double _x, double _y, traj_type traj) const { 
        return traj<0? 1: constant; 
    }

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    inline void eval(double _x, double _y, traj_type traj, Grad &g) const {
        g.f= (traj<0) ? 1 : constant;
    }

    /** calculate f(x) and df/dz, where z are the parameters */
    inline double evalGrad(double _x, double _y, traj_type traj, double partial[]) const {
        return (traj<0) ? 1 : constant;
    }
/*** end VDynamics2d ***/

    double eval(int t) const { return constant; }
    /** calculate f(x) */

    virtual double evalCore(int t) const { 
        return id[t]<0? 1: constant; 
    }
    virtual double evalCoreVar(int t) const { 
        return id[t]<0? 1: constant; 
    }

    virtual void eval(int t, Grad &g) const { g.f= constant; }

    virtual double evalGrad(int t, double p[]) const {
        return id[t]<0? 1: constant; 
    }

    /// For filtering algorithm. Evaluate without moving.
//    virtual double evalGradNoMove(int t, double p[]) const 
//    {
//        return id[t]<0? 1: constant; 
//    }

//    virtual double evalGrad(int t, double x, traj_type tr, double p[]) const {
//        return evalGrad(t,p);
//    }
    virtual double evalGrad(int t, double x, double y,
			      traj_type tr, double p[]) const {
        return evalGrad(t,p);
    }

    virtual double getMinZ() const { return constant; }
    virtual double getMaconstant() const { return constant; }

    virtual int getNParam() const { return 0; };
    virtual int getNAllParam() const { return 0; };

    virtual bool getParam(int t, double p[], int a[]) const {
        return true; 
    };
    virtual void getAllParam(int t, double p[]) const { };

    virtual void setParam(int t, double p[]) { }
    virtual void setAllParam(int t, double p) { constant= p; hippo_Print(constant); }
    virtual void setAllParam(int t, double p[]) { }
    virtual void setInitialEst(int t) { };
    virtual void setValues(int t, double z) { constant= z; } 
    
    virtual void updateDynamics(int tStart, int tEnd, double diff[])  {
    }

    virtual void updateModel(int tStart, int tEnd, double diff[]) { };
    
    virtual void allocForConverge() { }
    virtual void copyForConverge() { }
    virtual bool hasConverged() { return true; }
    virtual Const* link_copy() { 
        Const *c= new Const;
        VDynamics::copyInit(c);
        c->constant= constant;
        c->conv= conv; c->convper= convper;
        return c; 
    };

    virtual void reestimate(AdaptModel *);
    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
    
    virtual double evalAtIsi(int t, double isi) { return eval(t); }
};

} // namespace AFilter

#endif   // CONST_H
