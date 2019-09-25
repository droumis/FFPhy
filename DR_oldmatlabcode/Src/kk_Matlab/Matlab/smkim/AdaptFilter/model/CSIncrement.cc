/* $Id: CSIncrement.cc,v 1.6 2008/06/27 00:21:28 chengs Exp $
   
   Sen Cheng, Wed Oct  6 08:39:54 PDT 2004
   implementation of class CSIncrement
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSIncrement.h"
#include "Operators.h"

/* ************************************************************
                          class CSIncrement
   ************************************************************ */

using namespace AFilter;

CSIncrement::CSIncrement()
{
    updates= initial= 0;
    allocated.updates= false;
    allocated.initial= false;
    convCopy= 0;
    tIni= -1;

    trajs= 0;
    inds= 0;
}

CSIncrement::~CSIncrement()
{
    if (allocated.updates && updates) free(updates);
    if (allocated.initial && initial) free(initial);

    if (allocated.trajs && trajs) free(trajs);
    if (allocated.inds && inds) free(inds);
    if (convCopy) delete convCopy;
}

CSIncrement* CSIncrement::alloc_copy() const
{
    CSIncrement *tmp= new CSIncrement;
    tmp->init(startindex, endindex, super);
    return tmp;
}

void CSIncrement::copyForConverge()
{
    hippo_Assert(convCopy, "conCopy was not allocated");
    memcpy(convCopy->updates, updates, T*nUpdate*sizeof(double));
    memcpy(convCopy->initial, initial, nTotal*sizeof(double));
    convCopy->tIni= tIni;
    convCopy->tCurr= tCurr;
}

void CSIncrement::averageWithCopy()
{
    hippo_Assert(convCopy, "conCopy was not allocated");

    // average initial
    moveTo(convCopy->tIni);
    setInitialFromCurrent();
    for(int i=0; i<nTotal; ++i)
        initial[i]= 0.5*(initial[i]+convCopy->initial[i]);
    setToInitial();

    // average updates
    for(int i=0; i<T*nUpdate; ++i)
        updates[i]= 0.5*(updates[i]+convCopy->updates[i]);
}


void CSIncrement::init(int _startindex, int _endindex, CSFunction *f)
{
    CSStorage::init(_startindex, _endindex, f);
    if(!allocated.updates && !updates) {
        hippo_CAlloc(updates, T*nUpdate, double);
        allocated.updates= true;
    }
    if(!allocated.initial && !initial) {
        hippo_CAlloc(initial,nTotal, double);
        allocated.initial= true;
        tIni= startindex;
    } else
        setToInitial();

    hippo_CAlloc(trajs, f->getNFct(), traj_type);
    allocated.trajs= true;

    hippo_CAlloc(inds, nUpdate, int);
    allocated.inds= true;
}   

void CSIncrement::calcMinMax(const VOp *op) const
{
    double tmp;
    min= DBL_MAX; max= -DBL_MAX;

    //extreme parameter value at first timestep
    moveTo(startindex);
    for(int i=0; i<nFct; i++)
        for(int j=0; j<nPerFct[i]; j++) {
            tmp= op->op(param[i][j]);
            if(tmp < min) min= tmp;
            if(tmp > max) max= tmp;
        }

    //look for extreme values after updates in every timestep
    for(int t=startindex; t<endindex; t++)  {
        moveTo(t); // updates variables: param, nt, trajs and inds
        if(nt <= 0) continue;
        for(int i= 0; i < nUpdate; i++)  {
            tmp= op->op(param[trajs[0]][inds[i]]);
            if(tmp < min) min= tmp;
            if(tmp > max) max= tmp;
        }
    }
}

void CSIncrement::setInitialFromCurrent()
{
    double *tmp= initial;
    for (int i = 0; i < nFct; i++)  {
        memcpy(tmp, param[i], nPerFct[i]*sizeof(double));
        tmp+= nPerFct[i];
    }
    tIni= tCurr;
}

void CSIncrement::setToInitial()
{
    double *tmp= initial;
    for (int i = 0; i < nFct; i++)  {
        memcpy(param[i], tmp, nPerFct[i]*sizeof(double));
        tmp+= nPerFct[i];
    }
    tCurr= tIni;
}

void CSIncrement::resetInitialTo(int t)
{
    setToInitial();
    getAllParam(t, initial);
    tIni= t;
}

void CSIncrement::getAllParam(int t, double z[]) const
{
    moveTo(t);
    for (int i = 0; i < nFct; i++) {
        memcpy(z, param[i], nPerFct[i]*sizeof(double));
        z+= nPerFct[i];
    }
}

void CSIncrement::setAllParam(int t, double z[])
{
    for (int i = 0; i < nFct; i++)  {
        memcpy(param[i], z, nPerFct[i]*sizeof(double));
        z+=nPerFct[i];
    }
    tCurr= t;
    setInitialFromCurrent();
}

void CSIncrement::setAllParam(int t, double z)
{
    for(int i = 0; i < nFct; i++) 
        for(int j = 0; j < nPerFct[i]; j++)  {
            param[i][j]= z;
        }
    tCurr= t;
    setInitialFromCurrent();

    hippo_Assert(initial[0]== z, "setAllParam was not successful");
}

void CSIncrement::moveTo(int t) const
{
    if(t== tCurr) return;

    if(t > tCurr) { // move forward
        for( ; tCurr < t; tCurr++) {
            super->getActiveIndices(tCurr, nt, trajs, inds);
            if(nt <= 0) continue;
            double *u= updates+nUpdate*tCurr;
            for(int it= 0; it < nt; it++) 
                for(int i= 0; i < nUpdate; i++)  {
                    param[trajs[it]][inds[i]]+= u[i];
                }
        }
    } else {       // move backwards
        for(tCurr-= 1; tCurr >= t; tCurr--) {
            super->getActiveIndices(tCurr, nt, trajs, inds);
            if(nt <= 0) continue;
            double *u= updates+nUpdate*tCurr;
            for(int it= 0; it < nt; it++) 
                for(int i= 0; i < nUpdate; i++)  {
                    param[trajs[it]][inds[i]]-= u[i];
                }
        }
        tCurr++;
    }

    hippo_Assert(t==tCurr, "current time messed up");
}


void CSIncrement::updateEst(int tStart, int tEnd, double diff[])
{
    hippo_Assert(tEnd >= startindex && tEnd <= endindex, "tEnd out of bounds");
    int dir= tEnd-tStart;
    hippo_Assert(dir==1 || dir==-1, "-1 or 1");

    int t= (tStart < tEnd)? tStart : tEnd;

    //@@ what if traj at tStart and tEnd are not the same? -> mapping ?
    const int kOffset= nUpdate*t;
    if(dir==1) {
        for (int k= 0; k < nUpdate; k++) updates[kOffset+k]= diff[k];
    } else {
        for (int k= 0; k < nUpdate; k++) updates[kOffset+k]= -diff[k];
    }

    if(tEnd== startindex || tEnd== endindex) resetInitialTo(tEnd);
}

inline double minVal(const double a, const double b) { return a<b? a : b; }

void CSIncrement::maxDiff(double &abs, double &rel) const
{
    hippo_Assert(convCopy, "conCopy was not allocated");
    double diff;
//    double norm;
    abs= rel= 0;
    for(int k= 0; k < nUpdate*T; k++) {
        diff= fabs(updates[k]- convCopy->updates[k]);
        if (diff > abs) abs= diff;
//        norm= minVal(fabs(updates[k]) , fabs(convCopy->updates[k]));
//        norm= fabs(updates[k]);
//        if(norm > 1e-6) {
//            cerr << diff << '\t' << norm << '\n';
//            diff/= norm;
//            if (diff > rel) rel= diff;
//        }
    }
}

void CSIncrement::setMexInput(const mxArray *in)
{
    mxArray *mxtmp;
    // get  parameter values
    if(allocated.updates && updates) {
        free(updates);
        updates= 0;
        allocated.updates= false;
    }
    int tmp;
    mxtmp= mxGetField(in, 0, "updates");
    if(mxtmp) {
        MX_FieldAssign(in, "updates", updates, tmp);
//        hippo_Assert(T*nUpdate== tmp,
//                "updates in input has incorrect number of parameters");
    } 

    if(allocated.initial && initial) {
        free(initial);
        initial= 0;
        allocated.initial= false;
    }
    mxtmp= mxGetField(in, 0, "initial");
    if(mxtmp) {
        MX_Assign(mxtmp, initial);
//        hippo_Assert(mxGetM(mxtmp)*mxGetN(mxtmp)== nTotal,
//                "initial has to have nTotal parameters");
    } 
    MX_FieldScalarDefault(in, "tIni", tIni, startindex, int);

//    hippo_Print(initial);
//    hippo_Print(tIni);
}    

mxArray* CSIncrement::getMexOutput()  const
{
    mxArray *mxtmp;
    const char *fields[]= { "updates", "initial", "tIni" };
    const int nFields= sizeof(fields)/sizeof(fields[0]);
    mxArray *out= mxCreateStructMatrix(1,1, nFields,  fields);

    double *dblPtr;
    MX_CreateDoubleAssign(mxtmp, nUpdate,T, dblPtr);
    memcpy(dblPtr, updates, T*nUpdate* sizeof(double));
    mxSetField(out, 0, "updates", mxtmp);

    //    for(int i=0; i< 16; i++) cout << updates[i] << '\t';

    MX_CreateDoubleAssign(mxtmp, nTotal,1, dblPtr);
    memcpy(dblPtr, initial, nTotal* sizeof(double));
    mxSetField(out, 0, "initial", mxtmp);

//    hippo_Print(tIni);
    mxSetField(out, 0, "tIni", mxCreateScalarDouble(tIni));


    return out;
}

