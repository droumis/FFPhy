/* $Id: RescalingGen.cc,v 1.8 2008/08/24 19:47:01 chengs Exp $
   
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
    /* initialize the random number generator */
    struct timeval    tmptime;
    struct timezone    tmpzone;
    gettimeofday(&tmptime, &tmpzone);
    srandom(tmptime.tv_sec);
}

RescalingGen::~RescalingGen()
{
}

// generate spikes using time rescaling theorem
void RescalingGen::runGenerator(TData *data, AdaptModel *model)
{
    traj_type dir=0; /* assume that this is to go in the positive motion direction ??@@*/
    double dT= data->deltaTime;
    double currenttime;
    double lastspiketime= -1;
    double LambdaT, LambdaTInt= 0;

    double lastlambda_2= 0.5*model->eval(0);
    init(data,model);
    // run through time steps
    int t= 0; 
    while (t < data->ntimesteps) {
        /* draw from an exponential distribution with a parameter of one */ 
        LambdaTInt+= log((double) random() / (double) RAND_MAX);

        // integrate conditional prob. until greater than -log(rand)
        while (LambdaTInt<0 && t<data->ntimesteps) {

            currenttime= data->timearray[t];

            // trajectory is as selected
            if (data->fieldID[t] != dir) continue;

            LambdaT= lastlambda_2;

            if(lastspiketime < 0)
                lastlambda_2= 0.5*model->eval(t);
            else
                lastlambda_2= 0.5*model->evalAtIsi(t, currenttime-lastspiketime);
            LambdaT+= lastlambda_2;
            
            // integrate probability
            LambdaTInt += LambdaT * dT; 

            t++;
        }

        // then generate spike
        if (LambdaTInt > 0 && t <= data->ntimesteps) {
            /* put a spike at the current time */
            appendSpike(currenttime, t-1);
            lastspiketime= currenttime;
        }
    }
}
