/* $Id: VStat1d.cc,v 1.4 2008/08/24 19:47:06 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VStat1d.h"
#include "../model/IsiModel.h"

/* ************************************************************
                          class VStat1d
   ************************************************************ */

using namespace AFilter;

VStat1d::VStat1d(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat(d,m, l)
{
}

VStat1d::VStat1d(VStat1d &)
{
    hippo_Empty;
}

VStat1d::~VStat1d()
{
}

void VStat1d::init_obj()
{
    const char *name= model->getName();
    if (strcmp(name,"Pos_Isi")==0) {
        dyn1= dynamic_cast<const IsiModel*>(model)->get1dSpatialModel();
    } else {
        printf("stat= '%s', model= '%s'\n", this->getName(), name);
        hippo_ErrorMsg("The statistic above cannot be used with the model.");
        dyn1= 0;
    }
//    hippo_Mark;
//    hippo_Print(getName());
}

