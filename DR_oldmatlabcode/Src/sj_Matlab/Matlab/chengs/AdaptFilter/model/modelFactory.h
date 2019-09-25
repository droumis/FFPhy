#ifndef MODELFACTORY_H
#define MODELFACTORY_H

/* $Id: modelFactory.h,v 1.3 2008/08/24 19:47:05 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/


#include "../aux/defs.h"
#include "AdaptModel.h"
#include "VDynamics1d.h"
#include "VDynamics2d.h"

/* ************************************************************
   ************************************************************ */

namespace AFilter {

    AdaptModel* allocModel(const mxArray *opts, TData *data= 0, bool warn= true);

    VDynamics* allocDynamics(const mxArray *opts, TData *data= 0, bool warn= true);
    VDynamics1d* allocDynamics1d(const mxArray *opts, TData *data= 0, bool warn= true);

    VDynamics2d* allocDynamics2d(const mxArray *opts, TData *data= 0, bool warn= true);

} // namespace AFilter

#endif   // MODELFACTORY_H
