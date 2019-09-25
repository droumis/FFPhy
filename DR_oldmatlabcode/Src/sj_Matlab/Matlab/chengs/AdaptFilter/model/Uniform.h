#ifndef UNIFORM_H
#define UNIFORM_H

/* $Id: Uniform.h,v 1.5 2008/08/24 19:47:04 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../model/VDynamics1d.h"
#include "../model/VDynamics2d.h"

/* ************************************************************
                          class Uniform
   ************************************************************ */

namespace AFilter {

/** NOT implemented correctly. Please fix before using.
 */
class Uniform : public VDynamics1d, public VDynamics2d
{
private:
    int T;
    double *a;
public:
    Uniform();
    Uniform(Uniform &);
    virtual ~Uniform();

    const char* getName() const { return "Uniform"; }
    virtual void init(int _startindex, int _endindex, int _T, 
            double *x, int nid, traj_type *id) 
    {
        VDynamics1d::init(_startindex, _endindex, _T, _x, nid, _id) ;
    };
    virtual void init(int _startindex, int _endindex, int _T, 
            double *x, double *y, int nid, traj_type *id) 
    {
        VDynamics2d::init(_startindex, _endindex, _T, _x, _y, nid, _id) ;
    };
    
    
    virtual double eval(int t, double x, traj_type traj) const { return a[t];}
    virtual double eval(int t, double x, double y, traj_type traj) const
    { return a[t];}
    
    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(int t, double x,
			      traj_type traj, double partial[]) const
    { partial[0]=1/a[t]; return a[t];}

    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(int t, double x, double y,
			      traj_type traj, double partial[]) const
    { partial[0]=1/a[t]; return a[t];}
    

    virtual int getNParam() { return 1; };
    virtual void setInitialToCurrent(int _tIni) { }
    virtual void setInitialEst(int t) const { };
    
    virtual void updateDynamics(int tStart, int tEnd, double diff[]){hippo_Empty; };
    
    void copyTo(double p[]) {hippo_Empty;}
    virtual void getMaxDiff(double p[], double &maxdiff, double &reldiff, 
			    double cut) {hippo_Empty;}

    virtual Uniform* link_copy() { hippo_Empty; return 0; };
    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
    
};

} // namespace AFilter

#endif   // UNIFORM_H
