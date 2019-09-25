/* $Id: RescalingGen.cc,v 1.10 2008/10/23 21:24:10 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <sys/time.h>
#include <math.h>
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
#include "../aux/TData.h"
#include "../model/AdaptModel.h"

#include "RescalingGen.h"

/* ************************************************************
                          class RescalingGen
   ************************************************************ */

using namespace AFilter;

RescalingGen::RescalingGen()
{
}

RescalingGen::~RescalingGen()
{
}

// generate spikes using time rescaling theorem
void RescalingGen::runGenerator(TData *data, AdaptModel *model)
{
    traj_type dir=0; /* assume that this is to go in the positive motion direction ??@@*/
    int t= 0; 
    double dT= data->deltaTime;
    double currenttime= data->timearray[t];
    double lastspiketime= -1;
    double LambdaT;

    double lastlambda_2= 0.5*model->eval(0);
    init(data,model);

    // Draw exponential distribution with a parameter of one, subtract from
    // LambdaTInt
    // Integrate conditional prob. in LambdaTInt
    // Every time LambdaTInt crosses 0, generate a spike.

    double LambdaTInt= log((double) random() / (double) RAND_MAX);
    
    for (int t= 0; t < data->ntimesteps; t++) {

        // trajectory is as selected
        //            if (data->fieldID[t] != dir) continue;
        if (data->fieldID[t] < 0) continue;

        currenttime= data->timearray[t];

        LambdaT= lastlambda_2;

        if(lastspiketime < 0)
            lastlambda_2= 0.5*model->eval(t);
        else
            lastlambda_2= 0.5*model->evalAtIsi(t, currenttime-lastspiketime);
        LambdaT+= lastlambda_2;

        // integrate probability
        LambdaTInt += LambdaT * dT; 

        // then generate spike
        if (LambdaTInt > 0) {
            /* put a spike at the current time */
            appendSpike(currenttime, t);
            lastspiketime= currenttime;
            /* draw from an exponential distribution with a parameter of one */ 
            LambdaTInt+= log((double) random() / (double) RAND_MAX);
        }
    }
}

void RescalingGen::setMexInput(const mxArray *in)
{
    MX_FieldScalarDefault(in, "rand_seed", rand_seed, -1, unsigned int);

    if (rand_seed<0) {
        /* initialize the random number generator */
        struct timeval    tmptime;
        struct timezone    tmpzone;
        gettimeofday(&tmptime, &tmpzone);
        srandom(tmptime.tv_sec);
    } else {
        srandom(rand_seed);
    }
}
