/* $Id: filterFactory.cc,v 1.10 2008/09/01 21:07:44 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"

#include "filterFactory.h"
#include "AscentFilter.h"
//#include "AFSmooth.h"
//#include "EM.h"
//#include "KalmanSmooth.h"

/* ************************************************************
   ************************************************************ */

using namespace AFilter;

VFilter* AFilter::allocFilter(const mxArray *opts, TData *data, AdaptModel* model)
{
    VFilter *filter;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "filter.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     


    if (strcmp(name,"AscentFilter")==0) {
        filter= new AFilter::AscentFilter(data, model);
//    } else if (strcmp(name,"AFSmooth")==0) {
//        filter= new AFilter::AFSmooth(data, model);
//    } else if (strcmp(name,"KalmanSmooth")==0) {
//        filter= new AFilter::KalmanSmooth(data, model);
//    } else if (strcmp(name,"EM")==0) {
//        filter= new AFilter::EM(data, model);
    } else {
        filter= 0;
        MX_Error("unknown adaptive filter algorithm: %s ", name);
    }
    filter->setMexInput(opts);
    filter->init();

    mxFree(name);
    return filter;
}
