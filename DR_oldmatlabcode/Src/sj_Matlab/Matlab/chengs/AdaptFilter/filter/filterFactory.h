#ifndef FILTERFACTORY_H
#define FILTERFACTORY_H

/* $Id: filterFactory.h,v 1.2 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/


#include "../aux/defs.h"
#include "../model/AdaptModel.h"
#include "../gen/VGenerator.h"

#include "VFilter.h"

/* ************************************************************
   ************************************************************ */

namespace AFilter {

    VFilter* allocFilter(const mxArray *opts, TData *data, AdaptModel* allocModel);

} // namespace AFilter

#endif   // FILTERFACTORY_H
