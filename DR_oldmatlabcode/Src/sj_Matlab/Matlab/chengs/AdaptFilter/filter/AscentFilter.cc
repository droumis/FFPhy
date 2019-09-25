/* $Id: AscentFilter.cc,v 1.10 2008/10/23 21:23:57 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "AscentFilter.h"

/* ************************************************************
                          class AscentFilter
   ************************************************************ */

using namespace AFilter;

AscentFilter::AscentFilter(TData *_data, AdaptModel *_model)
: alternatePass(false), nparam(0), eps(0), data(_data)
{
    model= dynamic_cast<IsiModel*>(_model);
    hippo_Assert(model, "AscentFilter only works with classes derived from IsiModel.");
}

AscentFilter::~AscentFilter()
{
}

void AscentFilter::init()
{
    hippo_Assert(data, "data object not set");
    hippo_Assert(model, "model object not set");
    hippo_Assert(model->getNParam()== nparam, 
            "eps vector and model should have same number of elements");

    model->allocForConvergence();

    hippo_Assert(data->getNValidSpikes(), "data contains no spikes");
}

void AscentFilter::runFilter()
{
    // run the adaptive filter algorithm until convergence or limit of
    // iterations is reached

    bool converged=false;
    if (ran_iterations==0)
        model->setInitialEst(model->getEndindex());

    cerr << "iteration " ;
    while ((ran_iterations < niter) && !converged) {
        ++ran_iterations;

        cerr << ran_iterations << " ";

        model->copyForConverge();

//    int nisi= model->getIsiModel()->getNParam();
//    int nspatial= model->getSpatialModel()->getNParam();

        /* we first run the filter backwards to get initial estimate, and then forwards */
        if(alternatePass) {
//            cerr << "\t\t - alternating passes" << endl;
//            cerr << "pass back,T\n";
//            if (nspatial>0)
                runPass(-1, 1); // backward pass with temporal parameters fixed 
//            cerr << "pass back,S\n";
//            if (nisi>0)
                runPass(-1, 0); // backward pass with spatial parameters fixed 
//            cerr << "pass forward,S\n";
//            if (nspatial>0)
                runPass(1,0);   // forward pass with spatial parameters fixed 
//            cerr << "pass forward,T\n";
//            if (nisi>0)
                runPass(1,1);   // forward pass with temporal parameters fixed 
        } else {
            cerr << "\t\t - combined passes" << endl;
            runPass(-1, -1);                // backward pass
            runPass(1, -1);                 // forward pass
        }
        converged= model->hasConverged();
    }
    cerr << endl;
}

void AscentFilter::runPass(int dir, int fixing)
{
    double            lambda;  // values of lambda
    double            dT= data->deltaTime;
    double            inno, *diff= new double[nparam], mag;
    double mult;
    traj_type         id;
    int               dNt, dt;
    int i, t, begin, end;

    hippo_Assert(dir== 1 || dir == -1, "invalid direction");

    if (dir== 1) { // run forward pass
        begin= model->getStartindex(); end= model->getEndindex(); 
        mult= forwardMult;
        dt= 0;
    } else if(dir== -1) {       // run backward pass
        begin=  model->getEndindex(); end= model->getStartindex();
        mult= backMult; 
        dt= -1;
    } else {
        hippo_Error("unknown dir= ", dir);
        return;
    }

    for(t= begin; t != end; t+= dir) { 
        model->moveTo(t);

        id = data->fieldID[t+dt];

        if(id== -1) {
            // trajectory not valid, do not update splines
            /* we should not use the current point, as the fieldID is -1, so we
             * do not record any updates */

            for(i=0; i < nparam; i++) { diff[i]= 0; }
        } else { // valid trajectory, calc. & update function

            // calculate innovations
            lambda= model->evalNoMove(t+dt, diff);
            hippo_Assert(isfinite(lambda), "lambda is infinite");
            dNt= data->dN[t+dt];
            inno = dNt - (lambda*dT);

            // calculate ascent vector and its magnitude
            mag=0;
            for(i=0; i < nparam; i++) {
                diff[i]*= mult*eps[i]*inno;
                mag+= diff[i]*diff[i];
            }
            mag= sqrt(mag);

            if (mag > maxGradient) {
                // print debugging info about gradient
#ifdef DEBUG
                cout << "t= " << t << "\t" << "large= " << mag << "\n";
                for(i=0; i < nparam; i++) {
                    cout << diff[i] << "\t";
                }
                cout << '\n';
                hippo_Print(inno);
                hippo_Print(lambda);
#endif
                // rescale gradient that are too large
                for(i=0; i < nparam; i++) diff[i]*= (0.5*maxGradient/mag);
            }

        } // if valid traj
        model->updateModel(t, t+dir, diff, fixing);

    } // for t
    delete[] diff;
}

void AscentFilter::close()
{
}
void AscentFilter::setMexInput(const mxArray* MexParam)
{
    MX_FieldScalarDefault(MexParam, "niter", niter, 20, int);
    MX_FieldScalarDefault(MexParam, "ran_iterations", ran_iterations, 0, int);
//         hippo_Print(niter);
    hippo_Assert(niter < 1000,"too many iterations requested");

    MX_FieldScalarDefault(MexParam, "forward_mult", forwardMult, 1, double);
    MX_FieldScalarDefault(MexParam, "back_mult", backMult, 1, double);
    MX_FieldScalarDefault(MexParam, "maxGradient", maxGradient, 1000, double);
    MX_FieldScalarDefault(MexParam, "alternatePass", alternatePass, true, bool);

    MX_FieldAssign(MexParam, "eps", eps, nparam);

    //    cerr << nparam << '\n';
    //    for(int i= 0; i<nparam; i++) cerr << i << ":\t" << eps[i] << '\n';
//         hippo_Print(ran_iterations);
    //     hippo_Print(eps);
}

mxArray* AscentFilter::getMexOutput()
{
    mxArray *out= mxCreateStructMatrix(1,1, 0,  0);
    mxArray *mxtmp;

    mxAddField(out, "name");
    mxtmp= mxCreateString("AscentFilter");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);

    mxAddField(out, "niter");
    mxSetField(out, 0, "niter", mxCreateScalarDouble(niter));

    mxAddField(out, "ran_iterations");
    mxSetField(out, 0, "ran_iterations", mxCreateScalarDouble(ran_iterations));

    mxAddField(out, "forwardMult");
    mxSetField(out, 0, "forwardMult", mxCreateScalarDouble(forwardMult));

    mxAddField(out, "backMult");
    mxSetField(out, 0, "backMult", mxCreateScalarDouble(backMult));

    mxAddField(out, "maxGradient");
    mxSetField(out, 0, "maxGradient", mxCreateScalarDouble(maxGradient));

    mxAddField(out, "alternatePass");
    mxSetField(out, 0, "alternatePass", mxCreateLogicalScalar(alternatePass));

    double *dblptr;
    MX_CreateDoubleAssign(mxtmp,1,nparam,dblptr)
    memcpy(dblptr, eps, nparam * sizeof(double));
    mxAddField(out, "eps");
    mxSetField(out, 0, "eps", mxtmp);

    return out;
}


