/* $Id: EM.cc,v 1.6 2006/09/25 15:38:38 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "EM.h"
//#include <math.h>
#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 

/* ************************************************************
                          class EM
   ************************************************************ */

using namespace AFilter;

EM::EM(TData *_data, AdaptModel *_model)
{
    data= _data;
    model= _model;
    kalman= new KalmanSmooth(_data, _model);
}

EM::~EM()
{
    delete kalman;
}

void EM::init()
{
    hippo_Assert(data, "data object not set");
    hippo_Assert(model, "model object not set");
    kalman->init();
}

void EM::runFilter()
{
    int nTotal= kalman->getNTotal(); 

//    hippo_Print(n);
    double *scov= kalman->getSCov();
    double *Q= kalman->getQ();
    int *occ= kalman->getOcc();
    for(int iter= 0; iter <= niter; iter++) {
        cerr << "**** iteration " << iter << endl; 
        if(iter==0) kalman->runFilter();
        else  { 
            kalman->runPass(1); // temporal parameters fixed 
            kalman->runPass(0); // spatial parameters fixed 
        }
        if(iter<niter)  {
//        hippo_Print(nTotal);
        for(int i=0; i < nTotal; i++) { 
            Q[i]= occ[i]? scov[i]/ occ[i] : 0; 
        }
        cerr << "Q:  ";
        for(int i=0; i < nTotal-30; i++) cerr << Q[i] << '\t';
        cerr << "\n----\n";
//        cerr << "scov:  ";
//        for(int i=nTotal-20; i < nTotal; i++) cerr << scov[i] << '\t';
//        for(int i=0; i < nTotal; i++) cerr << scov[i] << '\t';
//        cerr << "\n----\n";
        // stop EM once convergence criterion is reached
//        kalman->model->reestimate(kalman->var);
        }
    }
}

void EM::close()
{
}

void EM::setMexInput(const mxArray *in)
{
    MX_FieldScalarDefault(in, "niter", niter, 10, int);
    kalman->setMexInput(in);
}


mxArray* EM::getMexOutput()
{
    mxArray *out= kalman->getMexOutput();
    mxAddField(out, "niter");
    mxSetField(out, 0, "niter", mxCreateScalarDouble(niter));
    return out;
}


