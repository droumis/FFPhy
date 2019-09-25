/* $Id: AdaptModel.cc,v 1.8 2008/08/24 19:47:02 chengs Exp $
   authors     : Sen Cheng
   created     : 2005/02/20
 */

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "AdaptModel.h"

#include <math.h>

using namespace AFilter;

void AdaptModel::calcRescaled (double *&rescaled, int &n) const
{
    double intStep= data->getTime(1)-data->getTime(0);

    // skipped spikes that occur at times earlier or later than the model is defined
    int i= data->getStartspike(), endspike= data->getEndspike();

    while(data->getSpikeTimeIndex(i) < startindex && i <= endspike) i++;
    while(data->getSpikeTimeIndex(endspike) > endindex && i <= endspike)  endspike--;

    // allocated space for rescaled spikes times
    n= endspike-i+1;
    if(n<=0) {
        n= 0;
        rescaled= 0;
        return;
    }
    rescaled= new double[n];

    // rescale valid spike times
    n= 0;
    double tLast= data->getTime(startindex); // time of last
    double tNext;                            // time of next spike
    for( ; i <= endspike; i++) {
        tNext= data->getSpikeTime(i);
        int jEnd= data->timeToIndex(tNext);
        if(data->fieldID[jEnd]<0) continue; //@@ reject invalid spikes

        int jStart= data->timeToIndex(tLast);

        // integrate the interspike interval 
//         first integration step covers only part of intStep
        double lambdaTInt= 0;
        this->moveTo(jStart);
        lambdaTInt+= this->eval(jStart)* (tLast-data->getTime(jStart));
        for(int j= jStart+1; j < jEnd; j++) {
            if(data->fieldID[j]<0) continue; 
            this->moveTo(j);
            lambdaTInt += this->eval(j)* intStep; 
        }
        // last integration step covers only part of intStep, add rest
        this->moveTo(jEnd);
        lambdaTInt += this->eval(jEnd)*(tNext-data->getTime(jEnd));

        rescaled[n]= 1-exp(-lambdaTInt);
        n++;
        tLast= tNext;
    }
//    hippo_Print(n);
}


mxArray* AdaptModel::getMexOutput() 
{
    mxArray *out= mxCreateStructMatrix(1,1, 0,  0);

    if(real) {
        mxAddField(out, "rescaled_isi");
        // compute rescaled isi's
        double *rescaled;
        int n;
        calcRescaled(rescaled, n);
        // set rescaled output
        MX_AssignDouble(out, "rescaled_isi", n, 1, rescaled);
        if(n > 0) delete[] rescaled;
    }

    return out;
}
