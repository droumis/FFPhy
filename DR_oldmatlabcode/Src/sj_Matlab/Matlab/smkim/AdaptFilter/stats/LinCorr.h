/* $Id: LinCorr.h,v 1.3 2007/01/12 04:56:29 chengs Exp $
   
   definition of class AFilter::LinCorr
*/

#ifndef LinCorr_H
#define LinCorr_H

#include "VStat2d.h"

/* ************************************************************
                        class LinCorr
   ************************************************************ */

namespace AFilter {
    
/** Compute the linear correlation coefficient of 2-d probability distribution.
 *
 * Matlab syntax: 
 *  opts.name= 'LinCorr';<br>
 *  opts.sizex= *;<br>
 *  opts.sizey= *;<br>

    @author   Sen Cheng
*/
class LinCorr : public VStat2d
{
protected:
    // size of bins
    double sx, sy;

    // number of bins
    int nbinx, nbiny;
    
    // firing probabilty
    double **p;
    
    // marginal distributions
    double *px, *py;

    double *x, *y;
    double minx, maxx;
    double miny, maxy;
public:
    virtual void allocate();
    virtual void dealloc();
    virtual void zeroOut();

    virtual void setXY();
    virtual void normNmarg();
    virtual double calcCorrCoef();

public:
    LinCorr(const TData *d, const AdaptModel *m, const StatLimits *l);
    virtual ~LinCorr();
    
    virtual const char* getName() const { return "LinCorr"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // LinCorr_H
