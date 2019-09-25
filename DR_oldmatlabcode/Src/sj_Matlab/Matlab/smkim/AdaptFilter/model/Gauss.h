#ifndef GAUSSMODEL_H
#define GAUSSMODEL_H

/* $Id: Gauss.h,v 1.9 2008/09/01 21:07:44 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   Gaussian model for spatial firing rate.
   lambda(x)= exp(alpha - 1/2 (x-mu)' * sigma^-1 * (x-mu))

   lambda(t)= m;

   same for all trajectories @@ -> change

   ordering of param[6]:  
   0: alpha,   1: mu_x, 2: mu_y, 3: sigma_x, 4: sigma_y, 5: r
*/

#include "../model/VDynamics2d.h"

/* ************************************************************
                          class Gauss
   ************************************************************ */

namespace AFilter {

class Gauss : public VDynamics2d
{
private:
    double   *alpha;           //  log(max rate)[time]
    double   *mx;              //     mx[time]
    double   *my;              //     my[time]
    double   *Sx;              // std x[time]
    double   *Sy;              // std y[time]
    double   *r;               //      r[time]
    
    bool alloc;
public:
    Gauss();
    Gauss(Gauss &);
    virtual ~Gauss();
    
    
    virtual const char* getName() const {return "Gauss"; };
    
    virtual void init(int startindex, int endindex, int _T, double *x, double *y, int nid, traj_type *id);
    
    virtual void moveTo(int t) const { tCurr= t; };

    virtual double eval(int t, double x, double y, traj_type traj) const ;

    /** calculate f(x) */
    virtual double eval(double x, double y, traj_type id) const
    { return eval(tCurr, x, y, id);}
    
    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(double x, double y, traj_type traj, 
			      double partial[]) const ;
    
    virtual int getNParam() const { return 6; };
    virtual int getNAllParam() const { return 6; };

    virtual void setInitialToCurrent(int _tIni) {hippo_Empty;}
    virtual void setInitialEst(int t)const ;
    
    virtual void updateDynamics(int tStart, int tEnd, double diff[]);
    
    void copyTo(double p[]) {hippo_Empty;};
    virtual void getMaxDiff(double p[], double &maxdiff, double &reldiff, 
			    double cut) {hippo_Empty;};
    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

private:
    void dealloc();


};
} // namespace AFilter

#endif   // GAUSSMODEL_H
