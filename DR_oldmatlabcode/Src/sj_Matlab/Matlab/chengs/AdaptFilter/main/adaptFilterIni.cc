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
*/

/* 
   [output]= adaptFilter(data, initModel, modelParam, filterParam)
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check numbers of arguments */
    if ( (nlhs != 1) || ((nrhs != 4))) {
	mexErrMsgTxt("incorrect number of parameters, type \"help adaptFilter\" for usage\n");
    }

    // read inputs and initialize objects
    TData *data= new TData(prhs[0]);
    
    hippo_Mark;
    AdaptModel *inimodel= allocModel(prhs[1], data);
    hippo_Mark;
    AdaptModel *model= allocModel(prhs[2], data);
    hippo_Mark;
    model->setInitialFrom(inimodel);

    hippo_Mark;
    TFilter *filter= allocFilter(prhs[3], data, model);
    
    // run the adaptive filter algorithm until convergence or limit of
    // iterations is reached
//    hippo_Mark;
    hippo_Mark;
    filter->runFilter();
    
    // assign outputs
//    filter->close();
    hippo_Mark;
    plhs[0]= model->getMexOutput();
    hippo_Mark;

    delete filter;
    delete model;
    delete inimodel;
    delete data;
    return;
}
