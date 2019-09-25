#ifndef KalmanSmooth_H
#define KalmanSmooth_H

/* $Id: KalmanSmooth.h,v 1.5 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2006/06/01
   
   program description
*/

#include "VFilter.h"

/* ************************************************************
                          class KalmanSmooth
   ************************************************************ */

namespace AFilter {

class TData;
class AdaptModel;

/**  Kalman-smoothing-like algorithm to estimate states of a dynamical system
 * with Gaussian state-transitions and point process outputs.
 */
class KalmanSmooth : public VFilter
{
protected:
    TData *data;

    int nState, nTotal;

    double x1;
    double Q1, Qt;
    double *Q;

     /** sum E[ {x(t-1)-x(t)} {x(t-1)-x(t)}' ] */
    double *scov;

    int *nOcc; // "occupancy"
public: //@@
    AdaptModel *model;
    AdaptModel* var; // for variance of state estimate
public:
    KalmanSmooth(TData *data=0, AdaptModel *model=0);
    virtual ~KalmanSmooth();
    
    virtual const char* getName() const { return "KalmanSmooth"; };
    virtual void init();
    virtual void runFilter();
    virtual void close();

    virtual void runPass(int fixing);
    virtual int getNTotal() const { return nTotal; }
    virtual int getNState() const { return nState; }
     /** sum E[ {x(t-1)-x(t)} {x(t-1)-x(t)}' ] */
    virtual double* getSCov() const { return scov; }
    virtual double* getQ() const { return Q; }
    virtual int* getOcc() const { return nOcc; }

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

private:
};

} // namespace AFilter

#endif   // KalmanSmooth_H
