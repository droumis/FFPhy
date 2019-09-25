/* $Id: calcStatistics.cc,v 1.5 2008/08/24 19:47:06 chengs Exp $
   calcStatistics.c 
  
   Compute relevant statistics for given model and data.
   The model can derived from classes AdaptModel, VDynamics1d or VDynamics2d.
*/

#include <string.h>
#include <math.h>
#include <iostream>

#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"

#include "../aux/TData.h"
#include "../model/modelFactory.h"
#include "../stats/statsFactory.h"

using namespace std;
using namespace AFilter;

void Usage(void) ;

/******************************************************************************
  INTERFACE FUNCTION to be called from matlab
*/

/* 
   stats= calcStatistics(data, model, statsParam, limits)
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check numbers of arguments */
    if ( (nlhs != 1) | ((nrhs != 4))) {
        mexErrMsgTxt("incorrect number of parameters, type \"help calcStatistics\" for usage\n");
    }

    cerr << "\nCalculating quantities... \n" ;

    // read inputs and initialize objects
    TData *data= new TData(prhs[0]);
    StatLimits *limits= new StatLimits(prhs[3]);

    AdaptModel  *am= allocModel(prhs[1], data);

    const mxArray *sopts= prhs[2];
    VStat  **stats;
    int nstats= 0, ntmp= 0;

    // allocate statistics objects
    hippo_Assert(!mxIsEmpty(sopts), "no statistic objects defined in input"); 
    if(mxIsCell(sopts)) { // list of stats
        nstats= mxGetNumberOfElements(sopts);
        stats= new VStat*[nstats];
        for (int i= 0; i < nstats; i++) {
            const mxArray *mxtmp = mxGetCell(sopts, i);
            stats[ntmp]= allocStat(mxtmp,data,am,limits);
            if(stats[ntmp]) ntmp++;
        }
    } else {            // just one stat
        nstats= 1;
        stats= new VStat*;
        stats[0]= allocStat(sopts,data,am,limits);
        if(stats[0]) ntmp= 1;
    }
    nstats= ntmp;

    // calculate statistics
    for (int i= 0; i < nstats; i++) {
        stats[i]->run();
    }

    // output
    plhs[0]= collectMexOutput(nstats, stats);

    // free allocated memory
    for (int i= 0; i < nstats; i++) {
        delete stats[i];
    }
    delete[] stats;
    delete am;
    delete limits;
    delete data;

    cerr << "...done with stats." << endl;
    return;
}
