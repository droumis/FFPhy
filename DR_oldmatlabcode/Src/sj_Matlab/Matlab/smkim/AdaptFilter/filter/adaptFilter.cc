/*
    adaptFilter.c 
    2 D version for spatial and temporal spline

    Here we fix the temporal piece, estimate the spatial piece, fix the spatial dynaimics, and so on...


*/

#include <string.h>
#include <math.h>
#include <iostream>

#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"

#include "../aux/TData.h"
#include "../model/modelFactory.h"
#include "../filter/filterFactory.h"

using namespace std;
using namespace AFilter;

void Usage(void) ;

/******************************************************************************
  INTERFACE FUNCTION to be called from matlab

   model= adaptFilter(data, model, filter_opts)
   [model, filter]= adaptFilter(data, model, filter_opts)
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check numbers of arguments */
    if ( nlhs==0 || nlhs>2 || nrhs!=3) {
        mexErrMsgTxt("Incorrect number of parameters. Try \"help adaptFilter\".\n");
    }

    cerr << "\nRunning adaptive filtering ...\n";

    // read inputs and initialize objects
    TData *data= new TData(prhs[0]);
    AdaptModel *model= allocModel(prhs[1], data);
    VFilter *filter= allocFilter(prhs[2], data, model);
    
    // run the adaptive filter algorithm until convergence or limit of
    // iterations is reached
    cerr << "Filter: '" << filter->getName() << endl;
    filter->runFilter();
    
    // assign outputs
    plhs[0]= model->getMexOutput();
    if(nlhs>= 2) plhs[1]= filter->getMexOutput();

    delete filter;
    delete model;
    delete data;

    cerr << "... done with filter." << endl;
    return;
}
