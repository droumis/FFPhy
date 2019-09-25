/* $Id: VGenerator.cc,v 1.1 2008/08/24 19:47:02 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <sys/time.h>
#include <math.h>
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
#include "../aux/TData.h"
#include "../model/AdaptModel.h"
#include "VGenerator.h"

/* ************************************************************
                          class VGenerator
   ************************************************************ */

using namespace AFilter;

VGenerator::VGenerator()
    : nspikes(0), spiketimes(0), ispikes(0), nalloc(0)
{
    // seed random number generator
    struct timeval    tmptime;
    struct timezone    tmpzone;
    gettimeofday(&tmptime, &tmpzone);
    srandom(tmptime.tv_sec);
}

VGenerator::~VGenerator()
{
    if (spiketimes) free(spiketimes);
    if (ispikes) free(ispikes);
}

void VGenerator::init(TData *d, AdaptModel *model)
{
    traj_type dir= 0; //@@
    double deltaT= d->timearray[1] - d->timearray[0];
    
    // Expected number of spikes:   E(n)= \int \lambda(t) dt 
    double LambdaInt= 0;
    for (int t=0 ;t < d->ntimesteps; t++) {
        if (d->fieldID[t] != dir) continue;
        model->moveTo(t);
        LambdaInt+=  model->eval(t);
    }
    LambdaInt*= deltaT;

    hippo_Assert(LambdaInt >= 0, "lambda should be non-negative");
    hippo_Assert(LambdaInt < d->ntimesteps, "predicted too many spikes");
    
    nalloc= (int) round(LambdaInt);
    hippo_CAlloc(spiketimes, nalloc, double);
    hippo_CAlloc(ispikes, nalloc, int);

    nspikes= 0;
}

void VGenerator::appendSpike(double time, int tindex)
{
//     hippo_Print(nspikes);
//     hippo_Print(nalloc);
    hippo_Assert(nspikes <= nalloc, 
	      "should always have more memory allocated than spikes saved.");
    if (nspikes== nalloc) {
        nalloc*= 2;
        spiketimes= (double*) realloc(spiketimes, nalloc*sizeof(double));
        hippo_Assert(spiketimes, "could not allocated more needed memory");
        
        ispikes= (int*) realloc(ispikes, nalloc*sizeof(int));
        hippo_Assert(ispikes, "could not allocated more needed memory");
    }
    spiketimes[nspikes]= time;
    ispikes[nspikes]= tindex;
    nspikes++;
}    

mxArray* VGenerator::getMexOutput()
{
    mxArray *out;
    if (nspikes > 0) {
        const char *fields[]= { "spiketimes", "ispikes"};
        const int nFields= sizeof(fields)/sizeof(fields[0]);
        out= mxCreateStructMatrix(1,1, nFields,  fields);
        
        mxArray *mxtmp;
        double *ptr;


        
        MX_CreateDoubleAssign(mxtmp, 1, nspikes,ptr);
        memcpy(ptr, spiketimes, nspikes* sizeof(double));
        mxSetField(out, 0, "spiketimes", mxtmp);

        MX_CreateDoubleAssign(mxtmp, 1, nspikes,ptr);
        for(int i=0; i < nspikes; i++) ptr[i]= ispikes[i];
        mxSetField(out, 0, "ispikes", mxtmp);
    } else {
        MX_CreateDouble(out, 0, 0);
        hippo_Assert(out, "no matlab output allocated");
    }
    return out;
}

 
