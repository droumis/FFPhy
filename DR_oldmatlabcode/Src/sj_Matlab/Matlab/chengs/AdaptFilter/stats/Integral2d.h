/* $Id: Integral2d.h,v 1.1.1.1 2006/04/24 14:21:47 chengs Exp $
   
   definition of class AFilter::Integral2d
*/

#ifndef Integral2d_H
#define Integral2d_H

#include "VStat2d.h"

/* ************************************************************
                        class Integral2d
   ************************************************************ */

namespace AFilter {
    
/** Integrate 2-d probability distribution.
 *
 * Matlab syntax: 
 *  opts.name= 'Integral2d';<br>
 *  opts.dx= *;<br>
 *  opts.dy= *;<br>

    @author   Sen Cheng
*/
class Integral2d : public VStat2d
{
private:
    /// size of integration bins
    double dx, dy;
    
private:
    void zeroOut();
public:
    Integral2d(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~Integral2d();
    
    virtual const char* getName() const { return "Integral2d"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // Integral2d_H
