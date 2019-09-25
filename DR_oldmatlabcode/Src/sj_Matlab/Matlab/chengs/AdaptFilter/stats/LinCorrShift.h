/* $Id: LinCorrShift.h,v 1.4 2008/10/23 21:24:10 chengs Exp $
   
   definition of class AFilter::LinCorrShift
*/

#ifndef LinCorrShift_H
#define LinCorrShift_H

#include "LinCorr.h"

/* ************************************************************
                        class LinCorrShift
   ************************************************************ */

namespace AFilter {
    
/** Compute the correlation coefficient of linear-circular probability
 * distribution by shifting circular variable such that the correlation is
 * maximized.
 *
 * Matlab syntax: 
 *  opts.name= 'LinCorrShift';<br>
 *  opts.nbinx= *;<br>
 *  opts.nbiny= *;<br>

    @author   Sen Cheng
*/
class LinCorrShift : public LinCorr
{
protected:
    double *yreal;
    double offset;
    traj_type traj;
    int ana;
public:
    void allocate();
    void dealloc();
    void zeroOut();

    virtual void setOffset();
public:
    LinCorrShift(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~LinCorrShift();
    
    virtual const char* getName() const { return "LinCorrShift"; }
    
    virtual void run();
};
    
} // namespace AFilter

#endif   // LinCorrShift_H