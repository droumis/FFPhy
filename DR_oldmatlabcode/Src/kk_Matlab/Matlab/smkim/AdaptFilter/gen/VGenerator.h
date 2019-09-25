#ifndef TGENERATOR_H
#define TGENERATOR_H

/* $Id: VGenerator.h,v 1.1 2008/08/24 19:47:02 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/defs.h"

/* ************************************************************
                          class VGenerator
   ************************************************************ */

namespace AFilter {

class AdaptModel;
class TData;

/** Abstract base class for spike generation algorithms.
 */
class VGenerator
{
private:
    int nspikes;
    double *spiketimes;
    int *ispikes;

    int nalloc;

protected:
    VGenerator();

    virtual void init(TData *, AdaptModel *);
    virtual void appendSpike(double spiketime, int tindex);

public:
    virtual ~VGenerator();

    inline int getNSpikes() const { return nspikes; }
    /// main function: simulates spike trains
    virtual void runGenerator(TData *, AdaptModel *) =0;

    virtual void setMexInput(const mxArray *in) =0;

    virtual mxArray* getMexOutput();


};

}; // namespace AFilter

#endif   // TGENERATOR_H
