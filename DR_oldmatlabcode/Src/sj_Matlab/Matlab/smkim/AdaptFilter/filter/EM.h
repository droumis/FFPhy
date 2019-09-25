#ifndef EM_H
#define EM_H

/* $Id: EM.h,v 1.4 2008/08/24 19:47:00 chengs Exp $
   
   Sen Cheng, 2006/06/01
   
   program description
*/

#include "KalmanSmooth.h"

/* ************************************************************
                          class EM
   ************************************************************ */

namespace AFilter {

/** EM algorithm for estimating dynamics parameters.
 */
class EM : public VFilter
{
private:
    /// max number of EM iterations
    int niter; 

    TData *data;
    AdaptModel *model;

    /// Kalman Smoother object
    KalmanSmooth* kalman;
public:
    EM(TData *data=0, AdaptModel *model=0);
    virtual ~EM();

    virtual const char* getName() const { return "EM"; }
    
    virtual void init();
    virtual void runFilter();
    virtual void close();

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

private:
};

} // namespace AFilter

#endif   // EM_H
