#ifndef ASCENTFILTER_H
#define ASCENTFILTER_H

/* $Id: AscentFilter.h,v 1.7 2008/10/23 21:02:33 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VFilter.h"
#include "../model/IsiModel.h"
#include "../aux/defs.h"

/* ************************************************************
                          class AscentFilter
   ************************************************************ */

namespace AFilter {

class AscentFilter : public VFilter
{
protected:

    int niter;
    int ran_iterations;
    /* multipliers for the learning rate  */
    double forwardMult, backMult; 

    bool alternatePass;

    int nparam;
    double maxGradient;
    double *eps;
    double *conv, *convper;    // absolute and relative convergence criteria
    
    TData *data;
    IsiModel *model;

    IsiModel* prevmodel;
public:
    AscentFilter(TData *data=0, AdaptModel *model=0);
    virtual ~AscentFilter();
    
    virtual const char* getName() const { return "AscentFilter"; };
    virtual void init();
    virtual void runFilter();
    virtual void close();

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

protected:
    // fixing= {-1, 0, 1} for fixing no component, the first (spatial) component
    // and the second (isi) component
    void runPass(int dir, int fixing);
};

} // namespace AFilter

#endif   // ASCENTFILTER_H