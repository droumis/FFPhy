#ifndef GENFACTORY_H
#define GENFACTORY_H

/* $Id: genFactory.h,v 1.2 2008/08/24 19:47:02 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/


#include "../aux/defs.h"
#include "../model/AdaptModel.h"

#include "VGenerator.h"

/* ************************************************************
   ************************************************************ */

namespace AFilter {
    VGenerator* allocGenerator(const mxArray *opts, TData *data, AdaptModel* model);

} // namespace AFilter

#endif   // GENFACTORY_H
