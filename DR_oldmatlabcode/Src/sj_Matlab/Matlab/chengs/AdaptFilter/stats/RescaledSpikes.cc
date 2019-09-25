/* $Id: RescaledSpikes.cc,v 1.3 2006/12/06 20:04:16 chengs Exp $
   
   implementation of class AFilter::RescaledSpikes
*/

#include <math.h>

#include "RescaledSpikes.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class RescaledSpikes
   ************************************************************ */

using namespace AFilter;

RescaledSpikes::RescaledSpikes(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStatAM(d,m,l)
{
    intStep= 0;
    nSpikes= 0;
    rescaled= 0;
}

RescaledSpikes::~RescaledSpikes()
{
    if(rescaled) delete[] rescaled;
}

void RescaledSpikes::run()
{
    intStep= data->getTime(1)-data->getTime(0);
    // skipped spikes that occur at times earlier or later than the model is defined
    int i= data->getStartspike(), endspike= data->getEndspike();
    double tStart= data->getTime(model->getStartindex());
    double tEnd= data->getTime(model->getEndindex());
    while(i <= endspike && data->getSpikeTime(i) < tStart) i++;
    while(i <= endspike && data->getSpikeTime(endspike) > tEnd) endspike--;

    // allocated space for rescaled spikes times
    if(rescaled) delete[] rescaled;
    rescaled= new double[endspike-i+1];

    // rescale valid spike times
    nSpikes= 0;
    double tLast= tStart, tNext;  // time of last and next spike
    for( ; i <= endspike; i++) {
        tNext= data->getSpikeTime(i);
        int jEnd= data->timeToIndex(tNext);
        if(data->fieldID[jEnd]<0) continue; //@@ reject invalid spikes

        double lambdaTInt= 0;
        int jStart= data->timeToIndex(tLast);

        // integrate the interspike interval 
        // first integration step covers only part of intStep
        model->moveTo(jStart);
        lambdaTInt+= model->eval(jStart)* (tLast-data->getTime(jStart));
        for(int j= jStart+1; j < jEnd; j++) {
            if(data->fieldID[j]<0) continue; 
            model->moveTo(j);
            lambdaTInt += model->eval(j)* intStep; 
        }
        // last integration step covers only part of intStep, add rest
        model->moveTo(jEnd);
        lambdaTInt += model->eval(jEnd)*(tNext-data->getTime(jEnd));

        rescaled[nSpikes]= exp(-lambdaTInt);

        //@@debug
//        if(rescaled[nSpikes] > .9) {
//            printf("dt= %f, x= %f\n", tSpike-tStart, data->posarray[data->timeToIndex(t)]);
//        }

        nSpikes++;
        tLast= tNext;
    }

    // save analyses 
    init(1, nSpikes);
    for(i=0; i < nSpikes; i++) {
//        printf("i= %d,\tz= %f\n", i, rescaled[i]);
        set(0, i, 0, rescaled[i]);
    }
}

void RescaledSpikes::setMexInput(const mxArray *in)
{
    VStat::setMexInput(in);
//    MX_FieldScalarDefault(in, "intStep", intStep, 0.001, double);
//    hippo_Print(intStep);
}
