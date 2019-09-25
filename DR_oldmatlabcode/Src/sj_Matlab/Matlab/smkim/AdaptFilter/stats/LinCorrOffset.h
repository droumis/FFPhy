/* $Id: LinCorrOffset.h,v 1.1 2007/01/12 04:57:17 chengs Exp $
   
   definition of class AFilter::LinCorrOffset
*/

#ifndef LinCorrOffset_H
#define LinCorrOffset_H

#include "LinCorrShift.h"

/* ************************************************************
                        class LinCorrOffset
   ************************************************************ */

namespace AFilter {
    
/** Compute the linear correlation coefficient of 2-d probability distribution.
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
