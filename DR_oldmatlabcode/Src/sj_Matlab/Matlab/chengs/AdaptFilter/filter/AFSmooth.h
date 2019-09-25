#ifndef AFSmooth_H
#define AFSmooth_H

/* $Id: AFSmooth.h,v 1.2 2008/08/24 19:47:00 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VFilter.h"
#include "../aux/defs.h"
#include "AscentFilter.h"

/* ************************************************************
                          class AFSmooth
   ************************************************************ */

namespace AFilter {

class AFSmooth : public AscentFilter
{
protected:
public:
    AFSmooth(TData *data=0, AdaptModel *model=0);
    virtual ~AFSmooth();
    
    virtual const char* getName() const { return "AFSmooth"; };
    virtual void runFilter();

};

} // namespace AFilter

#endif   // AFSmooth_H
