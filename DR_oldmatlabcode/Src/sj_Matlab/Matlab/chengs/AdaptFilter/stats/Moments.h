/* $Id: Moments.h,v 1.2 2008/08/24 19:47:06 chengs Exp $
   
   definition of class AFilter::Moments
*/

#ifndef Moments_H
#define Moments_H

#include "VStat1d.h"

/* ************************************************************
                        class Moments
   ************************************************************ */

namespace AFilter {
    
/** Compute moments of a 1-d probability distribution.
 *
 * Inputs (in Matlab syntax): <br>
 *  opts.name= 'Moments';<br>
 *  opts.sizex= ;   % size of bin in x-variable (cm) <br>
 *
 * Outputs: <br>
 *     [area, mean (or center-of-mass), std, skewness]
 *
 * Needs from StatLimits: <br>
 *    linpos, traj

    @author   Sen Cheng
*/
class Moments : public VStat1d
{
private:
    // size of bins
    double sx;

    /// number of bins
    int nbinx;
    
    // firing probabilty
    double *p;
    
    // variables
    double *x;

private:
    void allocate();
    void dealloc();
    void zeroOut();
public:
    Moments(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~Moments();
    
    virtual const char* getName() const { return "Moments"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // Moments_H
