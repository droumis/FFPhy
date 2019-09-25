/* $Id: PoissonGen.cc,v 1.2 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
#include "../aux/TData.h"
#include "../model/AdaptModel.h"
#include "PoissonGen.h"

/* ************************************************************
                          class PoissonGen
   ************************************************************ */

using namespace AFilter;

PoissonGen::PoissonGen()
{
}

PoissonGen::~PoissonGen()
{
}

// void PoissonGen::setMexInput(const mxArray *in)
// {
// }

// mxArray* PoissonGen::getMexOutput()
// {
//     mxArray *out;
//     return out;
// }

//generate spikes using Poisson statistics
void PoissonGen::runGenerator(TData *d, AdaptModel *model)
{
    traj_type dir=0; /* assume that this is to go in the positive motion direction ??@@*/
    double dT= d->deltaTime;
    double currenttime;
    double lastspiketime= -1;
    double lambda, ran;

    init(d,model);

    // run through time steps
    for (int t=0 ;t < d->ntimesteps; t++) {

        // trajectory is as selected
        if (d->fieldID[t] != dir) continue;

        currenttime = d->timearray[t];

        if(lastspiketime < 0)
            lambda= model->eval(t);
        else
            lambda= model->evalAtIsi(t, currenttime-lastspiketime);

        ran= (double) random() / (double) RAND_MAX;
        if (ran < lambda*dT) {
            appendSpike(currenttime, t);
            lastspiketime= currenttime;
        }
    }
}
