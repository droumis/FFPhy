#ifndef GAUSSMIXMODEL_H
#define GAUSSMIXMODEL_H

/* $Id: GaussMix.h,v 1.5 2008/08/24 19:47:04 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   Gaussian model for spatial firing rate.
   lambda(x)= exp(alpha - 1/2 (x-mu)' * sigma^-1 * (x-mu))

   lambda(t)= m;

   same for all trajectories @@

   ultimately, split model into spatial and temporal -> PosPhase_Isi
*/

#include "../model/AdaptModel.h"
#include "Gauss.h"

/* ************************************************************
                          class GaussMix
   ************************************************************ */

namespace AFilter {

class GaussMix : public AdaptModel
{
private:
    int nMix;                  // number of mixture components
    Gauss *comp;          // Gaussian components
    
    double mult;

    TData *data;
    int T;                     // number of timesteps
    
    bool alloc;
    static const int nCompParam= 6; // number of parameters per component
public:
    GaussMix(TData *_data=0);
    ~GaussMix();
    
    virtual void init(double *x, double *y, traj_type *id);

    virtual double eval(int startindex, int endindex, int timeindex, 
            double x, double y, double isi, traj_type traj) const;
    
    virtual double eval(int startindex, int endindex, int timeindex, 
            double x, double y, traj_type traj) const;
    
    /* partial[6]:  0: alpha,   1: mu1,     2: mu2, 
                    3: sigma_x, 4: sigma_y, 5: r
		    ...
    */
//     virtual double evalPartial(int timeindex, double x, double y, 
// 			       traj_type traj, double partial[])const ;
//     virtual double evalLogPartial(int timeindex, double x, double y, 
// 			       traj_type traj, double partial[]) const;
    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(int timeindex, double x, double y, 
			       traj_type traj, double partial[]) const;
    
//    virtual TFilter* allocFilter(char []);

    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
      
    virtual int getNTimesteps() { return T; }
    virtual int getNParam() { return nCompParam*nMix; }
    virtual void setInitialEst()const ;
    virtual void updateDynamics(int tStart, int tEnd, double diff[]);

private:
    void dealloc();
};
} // namespace AFilter

#endif   // GAUSSMIXMODEL_H
