/* $Id: LinCorrOffset.h,v 1.3 2008/10/23 21:24:10 chengs Exp $
   
   definition of class AFilter::LinCorrOffset
*/

#ifndef LinCorrOffset_H
#define LinCorrOffset_H

#include "LinCorrShift.h"

/* ************************************************************
                        class LinCorrOffset
   ************************************************************ */

namespace AFilter {
    
/** Compute the correlation coefficient of linear-circular probability
 * distribution by shifting circular variable such that the *absolute*
 * value of the correlation is maximized.
 *
 * Matlab syntax: 
 *  opts.name= 'LinCorrOffset';<br>
 *  opts.nbinx= *;<br>
 *  opts.nbiny= *;<br>

    @author   Sen Cheng
*/
class LinCorrOffset : public LinCorrShift
{
protected:
public:
    virtual void setOffset();
public:
    LinCorrOffset(const TData *d, const AdaptModel *m, const StatLimits *l)
    : LinCorrShift(d,m,l) {};
    
    virtual const char* getName() const { return "LinCorrOffset"; }
};
    
} // namespace AFilter

#endif   // LinCorrOffset_H
