/* $Id: MutualInfo.h,v 1.2 2008/08/24 19:47:06 chengs Exp $
   
   definition of class AFilter::MutualInfo
*/

#ifndef MUTUALINFO_H
#define MUTUALINFO_H

#include "VStat2d.h"

/* ************************************************************
                        class MutualInfo
   ************************************************************ */

namespace AFilter {
    
/** Compute the mutual information of 2-d probability distribution.
 *
 * Inputs: <br>
 *  opts.name= 'MutualInfo';<br>
 *  opts.sizex= ;   % size of bin in x-variable (cm) <br>
 *  opts.sizey= ;   % size of bin in y-variable (cm) <br>
 *
 * Needs from StatLimits:
 *  linpos, phase, traj

    @author   Sen Cheng
*/
class MutualInfo : public VStat2d
{
private:
    /// number of bins
    int nbinx, nbiny;
    
    // firing probabilty
    double **p;
    
    // marginal distributions
    double *px, *py;

private:
    void allocate();
    void dealloc();
    void zeroOut();
public:
    MutualInfo(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~MutualInfo();
    
    virtual const char* getName() const { return "MutualInfo"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // MUTUALINFO_H
