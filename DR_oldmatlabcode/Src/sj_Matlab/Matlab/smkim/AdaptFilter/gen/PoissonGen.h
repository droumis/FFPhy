#ifndef POISSONGEN_H
#define POISSONGEN_H

/* $Id: PoissonGen.h,v 1.2 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   Generate data without history dependence

*/

#include "VGenerator.h"

/* ************************************************************
                          class PoissonGen
   ************************************************************ */

namespace AFilter {

class PoissonGen : public VGenerator
{
private:
public:
    PoissonGen();
    virtual ~PoissonGen();
    
    virtual void runGenerator(TData *, AdaptModel *);

    virtual void setMexInput(const mxArray *in) {};
};

}; // namespace AFilter

#endif   // POISSONGEN_H
