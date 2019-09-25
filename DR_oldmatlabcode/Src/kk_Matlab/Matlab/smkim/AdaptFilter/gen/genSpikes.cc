/* $Id: genSpikes.cc,v 1.2 2008/08/24 19:47:02 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <string.h> 

#include "../aux/TData.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"

#include "../model/modelFactory.h"
#include "../gen/genFactory.h"


/******************************************************************************
  INTERFACE FUNCTION to be called from matlab

  spiketimes= genSpikes(data, model, generator)
*/

using namespace AFilter;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /* Check numbers of arguments */
    if ( (nlhs != 1) || ((nrhs != 3))) {
	mexErrMsgTxt("incorrect number of parameters, type \"help adaptFilter\" for usage\n");
    }

    cerr << "\nGenerating spikes... \n" ;

    // read inputs and initialize objects
    // read inputs and initialize objects
    TData       *data= new TData(prhs[0]);
    AdaptModel *model= allocModel(prhs[1], data);
    VGenerator    *gen= allocGenerator(prhs[2], data, model);
    
    gen->runGenerator(data, model);
    
    plhs[0]= gen->getMexOutput();
    
    cerr << "... " << gen->getNSpikes() << " spikes";
    
    delete gen;
    delete model;
    delete data;
    cerr << " generated." << endl;
    return;
}
