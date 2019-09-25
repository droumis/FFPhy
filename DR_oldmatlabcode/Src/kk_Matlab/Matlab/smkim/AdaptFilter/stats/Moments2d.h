/* $Id: Moments2d.h,v 1.3 2008/08/24 19:47:06 chengs Exp $
   
   definition of class AFilter::Moments2d
*/

#ifndef Moments2d_H
#define Moments2d_H

#include "VStat2d.h"

/* ************************************************************
                        class Moments2d
   ************************************************************ */

namespace AFilter {
    
/** Compute moments of a 2-d probability distribution.
 *
 * Inputs (in Matlab syntax): <br>
 *  opts.name= 'Moments2d';<br>
 *  opts.sizex= ;   % size of bin in x-variable (cm) <br>
 *  opts.sizey= ;   % size of bin in y-variable (cm) <br>
 *  opts.circy= 0;  % treat y as circular variable (cm) <br>
 *
 * Outputs
 *  [meanx, meany, varx, vary, medianx, mediany, mutual info]
 *
 * Needs from StatLimits:
 *  linpos, phase, traj

    @author   Sen Cheng
*/
class Moments2d : public VStat2d
{
private:
    // size of bins
    double sx, sy;

    /// number of bins
    int nbinx, nbiny;
    
    // firing probabilty
    double **p;
    
    // marginal distributions
    double *px, *py;

    // variables
    double *x, *y;

    // whether y is circular variable
    bool circy;

private:
    void allocate();
    void dealloc();
    void zeroOut();
public:
    Moments2d(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~Moments2d();
    
    virtual const char* getName() const { return "Moments2d"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // Moments2d_H
