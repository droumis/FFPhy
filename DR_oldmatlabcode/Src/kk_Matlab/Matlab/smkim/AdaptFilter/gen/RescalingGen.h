#ifndef RESCALINGGEN_H
#define RESCALINGGEN_H

/* $Id: RescalingGen.h,v 1.2 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   Generate data using the time rescaling theorem 
     [cf. Brown EN, et al., Neural Comput. (2002) 14(2):325-46]
   simulate the data by generating a draw from a unit exponential distribution
   and integrating lambda until the integral equals or exceeds the draw

*/

#include "VGenerator.h"

/* ************************************************************
                          class RescalingGen
   ************************************************************ */

namespace AFilter {

class RescalingGen : public VGenerator
{
private:

public:
    RescalingGen();
    virtual ~RescalingGen();
    
    virtual void runGenerator(TData *, AdaptModel *);

    virtual void setMexInput(const mxArray *in) {};
    
};

}; // namespace AFilter

#endif   // RESCALINGGEN_H
