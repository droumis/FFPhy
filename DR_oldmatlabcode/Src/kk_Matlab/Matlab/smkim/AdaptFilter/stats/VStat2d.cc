/* $Id: VStat2d.cc,v 1.4 2008/08/24 19:47:06 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VStat2d.h"
#include "../model/IsiModel.h"

/* ************************************************************
                          class VStat2d
   ************************************************************ */

using namespace AFilter;

VStat2d::VStat2d(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat(d,m,l)
{
}

VStat2d::VStat2d(VStat2d &)
{
    hippo_Empty;
}

VStat2d::~VStat2d()
{
}

void VStat2d::init_obj()
{
    const char *name= model->getName();
    if (strcmp(name,"PosPhase_Isi")==0) {
        dyn2= dynamic_cast<const IsiModel*>(model)->get2dSpatialModel();
    } else {
        printf("stat= '%s', model= '%s'\n", this->getName(), name);
        hippo_ErrorMsg("The statistic above cannot be used with the model.");
        dyn2= 0;
    }
}

