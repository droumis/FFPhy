/* $Id: statsFactory.cc,v 1.7 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "statsFactory.h"
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"

//#include "../model/AdaptModel.h"
//#include "../model/VDynamics1d.h"
//#include "../model/VDynamics2d.h"

#include "MutualInfo.h"
#include "RescaledSpikes.h"
#include "LinCorr.h"
#include "LinCorrOffset.h"
#include "LinCorrShift.h"
#include "Integral2d.h"
#include "Moments.h"
#include "Moments2d.h"
#include "CircMeasures.h"


/* ************************************************************
   ************************************************************ */

using namespace AFilter;

VStat* AFilter::allocStat(const mxArray *opts, TData * data, AdaptModel *model,
        StatLimits *lim, bool warn)
{
    if(!opts || mxIsEmpty(opts)) return 0;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "field 'name' was not defined in stats");
    MX_StringAllocAssign(mxTmp,name);     

    VStat *stat= 0;
    if (strcmp(name,"RescaledSpikes")==0) {
        stat= new RescaledSpikes(data,model,lim);
    } else if (strcmp(name,"Moments")==0) {
        stat= new Moments(data,model,lim);
    } else if (strcmp(name,"Moments2d")==0) {
        stat= new Moments2d(data,model,lim);
    } else if (strcmp(name,"MutualInfo")==0) {
        stat= new MutualInfo(data,model,lim);
    } else if (strcmp(name,"LinCorr")==0) {
        stat= new LinCorr(data,model,lim);
    } else if (strcmp(name,"LinCorrOffset")==0) {
        stat= new LinCorrOffset(data,model,lim);
    } else if (strcmp(name,"LinCorrShift")==0) {
        stat= new LinCorrShift(data,model,lim);
    } else if (strcmp(name,"Integral2d")==0) {
        stat= new Integral2d(data,model,lim);
    } else if (strcmp(name,"CircMeasures")==0) {
        stat= new CircMeasures(data,model,lim);
    } else {
        if (warn) MX_Warn("Skipping unknown statistic '%s'", name);
    }

    if(stat) {
        stat->init_obj();
        stat->setMexInput(opts);
    }
    mxFree(name);
    return stat;
}

