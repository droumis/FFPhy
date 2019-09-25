/* $Id: CircMeasures.h,v 1.2 2007/01/12 04:56:29 chengs Exp $
   
   definition of class AFilter::CircMeasures
*/

#ifndef CircMeasures_H
#define CircMeasures_H

#include "LinCorr.h"

/* ************************************************************
                        class CircMeasures
   ************************************************************ */

namespace AFilter {
    
/** Compute circular measures of 2-d probability distribution.
 *
 * Matlab syntax: 
 *  opts.name= 'CircMeasures';<br>
 *  opts.sizex= *;<br>
 *  opts.sizey= *;<br>
 *
 * Outputs
 *  [mean change, mean dispersion, mean angle, change weighted by min, change w.
 *  by mean]

    @author   Sen Cheng
*/
class CircMeasures : public LinCorr
{
public:

public:
    CircMeasures(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~CircMeasures();
    
    virtual const char* getName() const { return "CircMeasures"; }
    
    virtual void run();
    
};
    
} // namespace AFilter

#endif   // CircMeasures_H
