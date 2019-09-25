/* $Id: AFSmooth.cc,v 1.2 2008/08/24 19:47:00 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "../model/IsiModel.h"
#include "AFSmooth.h"

/* ************************************************************
                          class AFSmooth
   ************************************************************ */

using namespace AFilter;

AFSmooth::AFSmooth(TData *_data, AdaptModel *_model)
: AscentFilter(_data, _model)
{
}

AFSmooth::~AFSmooth()
{
}

void AFSmooth::runFilter()
{
    // run AF until convergence
    AscentFilter::runFilter();

    // save last forward pass
    model->copyForConverge();

    // run backward pass
    if(alternatePass) {
        runPass(-1, 1); // temporal parameters fixed 
        runPass(-1, 0); // spatial parameters fixed 
    } else {
        runPass(-1, -1); // none fixed
    }

    // average forward and backward passes
    model->averageWithCopy();
}

