/* $Id: RescaledSpikes.h,v 1.1.1.1 2006/04/24 14:21:47 chengs Exp $
   File name   : RescaledSpikes.h
   authors     : Sen Cheng
   created     : Mon 08 Nov 2004 12:18:47 PM PST
  
 */
#ifndef RescaledSpikes_H
#define RescaledSpikes_H

#include "VStatAM.h"

/* ************************************************************
                        class RescaledSpikes
   ************************************************************ */

namespace AFilter {
    
/** Compute rescaled spike times z(n)= exp(-int ( lambda(t) dt)).
 * The integration runs from t(n-1) to t(n), i.e., from the time of the previous
 * spike to the time of the current spike.
 * By the time rescaling theorem P(z) should be uniform on [0,1]
 * Note, large interspike intervals t(n)-t(n-1) are mapped to small z(n), and
 * vice versa.
 *
 * opt.name=    'RescaledSpikes'; <br>
 * opt.intStep= 0.001;  % step size for integration [0.001] <br>
 *   
 @author   Sen Cheng
*/
class RescaledSpikes : public VStatAM
{
private:
    /** Width of integration step, default= 0.001 sec. */
    double intStep;

    /** Number of  spikes. */
    int nSpikes;

    /** Rescaled spike times. */
    double *rescaled;
private:
public:
    RescaledSpikes(const TData*, const AdaptModel*, const StatLimits*);
    virtual ~RescaledSpikes();
    
    virtual const char* getName() const { return "RescaledSpikes"; }
    
    virtual void run();
    
    virtual void setMexInput(const mxArray *in);
};
    
} // namespace AFilter

#endif   // RescaledSpikes_H
