#ifndef TFILTER_H
#define TFILTER_H

/* $Id: VFilter.h,v 1.1 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/defs.h"

/* ************************************************************
                          class VFilter
   ************************************************************ */

namespace AFilter {

class VFilter 
{
private:
public:
    VFilter() {};
    virtual ~VFilter() {};

    virtual const char* getName() const =0;
    virtual void init() =0;
    virtual void runFilter() =0;

    virtual void setMexInput(const mxArray *in=0) =0;
    virtual mxArray* getMexOutput() =0;
}; // class VFilter

}; // namespace AFilter

#endif   // TFILTER_H
