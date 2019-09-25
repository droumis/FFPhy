/* $Id: CSplinesNorm.h,v 1.5 2006/07/04 15:30:19 chengs Exp $
   File name   : /bach/AdaptFilter/model/CSplinesNorm.h
   authors     : Sen Cheng
   created     : Tue 02 Nov 2004 12:06:39 PM PST
  
 */
#ifndef CSplinesNorm_H
#define CSplinesNorm_H


#include "../aux/defs.h"
#include "CSplines.h"
#include "CSplinesFct.h"


/* ************************************************************
                          class CSplinesNorm
   ************************************************************ */

namespace AFilter {

/** normalized 1-d cardinal splines
 *  authors: Sen Cheng
 */
class CSplinesNorm : public CSplines {
protected:
    
    double **dNdTheta;      // gradient of norm[traj][csegx]

    mutable struct { // save when norm was calculated last
        double factor;  // 1/area of function
        int t;          // timestep 
        traj_type traj;   // trajectory
    } norm_last; 
protected:
    virtual void freeDNdTheta();
    virtual void allocDNdTheta();

    void zero_out();
    virtual void initObj();

    virtual void calc_norm(int t, traj_type traj) const;

//    virtual void getParam_p(traj_type traj, double p[], int segx= -1) const;
//    virtual void setParam_p(traj_type traj, double p[], int segx= -1) const;
public:
    CSplinesNorm();
    CSplinesNorm(CSplinesNorm &);
    virtual ~CSplinesNorm();

    virtual const char* getName() const {return "CSplinesNorm"; };

    bool  getParam(int t, double z[], int a[]) const;

    void getAllParam(int t, double z[]) const;
    // parameters are normalized before returned (so it could be used in regular
    // CardSplines
//    virtual void getTheta(int t, traj_type traj, double out[]) const;

    virtual double eval(double x, traj_type id)const ;
    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(double x, traj_type id, double partial[])const ;

    virtual double getIntegral(int timestep, traj_type traj) const;


};

} // namespace AFilters

#endif   // CSplinesNorm_H
