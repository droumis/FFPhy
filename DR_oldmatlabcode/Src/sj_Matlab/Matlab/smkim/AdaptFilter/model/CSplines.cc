/* $Id: CSplines.cc,v 1.12 2008/08/24 19:47:03 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <stdlib.h>
#include <string.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSplines.h"
//#include "CardSplines.h"
//#include "CSIncrement.h"
//#include "CSCompress.h"

/* ************************************************************
                          class CSplines
   ************************************************************ */

using namespace AFilter;

CSplines::CSplines()
{
    zero_out();
}

CSplines::CSplines(CSplines &)
{
    zero_out();
    hippo_Empty;
}

CSplines::~CSplines()
{
    dealloc();
}

void CSplines::zero_out()
{
    nUpdate= 4; T= 0; nFct=0; id=0; 
    ncpx= 0; cpx= 0; csegx= 0;

    allocated.controls= false;
}

void CSplines::dealloc() 
{
    if(allocated.controls) {
        if(nFct>1)  {
            if (ncpx) delete[] ncpx; 
            if (cpx) delete[] cpx; 
        } else if(nFct==1)  {
            if (ncpx) delete ncpx; 
            if (cpx) delete cpx; 
        }
        if (csegx) delete[] csegx;
        if(nPerFct) delete[] nPerFct; // is allocated in child class

    }
    for (int i = 0; i < nFct; i++) delete fctCurr[i];
//    if(fctCurr) delete[] fctCurr;
    fctCurr.erase(fctCurr.begin(), fctCurr.end());
    zero_out();
}

void CSplines::init(int _startindex, int _endindex, int _T, double *_x, 
        int _nFct, traj_type *_id) 
{
    VDynamics1d::init(_startindex, _endindex, _T, _x, _nFct, _id); 
    initObj(); // fct obj's have to be allocated before calling CSFunction::init
    if(!nPerFct) nPerFct= new int[nFct];
    for(int i=0; i<nFct; i++) nPerFct[i]= ncpx[i];

    CSFunction::init(_startindex, _endindex, _T, _nFct, _id);

    hippo_Assert(csegx==0, "csegx was previously assigned");
    csegx= new int[T];
    allocated.controls= true;

    traj_type tid;
    for(int i=startindex; i <= endindex; i++) {
        tid= id? id[i] : 0;
        if (tid < 0) {
            csegx[i]= -1;
        } else {
            // check whether control points cover all occuring values
            hippo_Assert(x[i] >= cpx[tid][1],"change control points");
            hippo_Assert(x[i] < cpx[tid][ncpx[tid]-2],"change control points");

            csegx[i]= fctCurr[tid]->findSegment(x[i]);

//            if(i==160) cerr << "t= " << i << ", csegx= "  << csegx[i] << "\n";
            hippo_Assert(csegx[i]>0 && csegx[i] < ncpx[tid]-1, "should not get here");
        }
    }

    // check wether mapping is correct
    if(mapping) {
        vector<mapItem>::iterator iter;
        for(iter= map.begin(); iter!= map.end(); ++iter) {
            tid= iter->t1;
            int i= 0;
            while( cpx[tid][i]< iter->min & i<ncpx[tid]) i++;
            while( cpx[tid][i]< iter->max & i<ncpx[tid]) {
                hippo_Assert(cpx[tid][i]== cpx[iter->t2][i], 
                        "mapped control points must be located at same values");
                i++;
            }
        }
    }
}   

double CSplines::getMinX() const 
{ 
//    hippo_Mark;
    if(!xmin.valid) {
        xmin= DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMinX() < xmin) xmin= fctCurr[traj]->getMinX();
//    hippo_Print(cpx[traj][1]);
        }
//    hippo_Print(xmin);
    }
    return xmin;
};

double CSplines::getMaxX() const 
{ 
//    hippo_Mark;
    if(!xmax.valid) {
        xmax= -DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMaxX() > xmax) xmax= fctCurr[traj]->getMaxX();
        }
    }
    return xmax;
};

//void CSplines::getParam_p(traj_type traj, double p[], int segx) const
//{
//    fctCurr[traj]->getParamPartial(p, segx);
//}

//void CSplines::setParam_p(traj_type traj, double p[], int segx) const
//{
//    fctCurr[traj]->setParamPartial(p, segx);
//}

//void CSplines::addToParam(traj_type traj, double p[], int segx) const
//{
//    fctCurr[traj]->addParamPartial(p, segx);
//}


//traj_type CSplines::getParam(int t, double z[], int a[]) const
//{
//    traj_type tid= CSFunction::getParam(t,z,a);
//    if(tid >= 0) {
//        int ind=0;
//        for(int i = 0; i < tid; i++) ind+= ncpx[i];
//        for(int i= 0; i<nUpdate; i++) a[i]+= ind;
//    }
//    return tid;
//}


void CSplines::initObj()
{
//    hippo_Mark;
//    fctCurr= new CSplinesFct*[nFct];
    // init. the current values 
    fctCurr.erase(fctCurr.begin(), fctCurr.end());
    fct.erase(fct.begin(), fct.end());
    for (int i= 0; i < nFct; i++) {
        CSplinesFct *tmp= new CSplinesFct;
        fctCurr.push_back(tmp);
        fct.push_back(tmp);
        tmp->init(ncpx[i], cpx[i]);
    }
}

void CSplines::getActiveIndices(int t, int &ntraj, traj_type trajs[], int ind[]) const
{
    ntraj= 0;
    trajs[0]= id? id[t] : 0;
    int tsegx= csegx[t];
    if (trajs[0] != -1 && tsegx!=-1) {
        ntraj++;
        fctCurr[trajs[0]]->getPartialIndices(ind, tsegx);
        if(!mapping) return;
        vector<mapItem>::const_iterator i;
        for(i= map.begin(); i!= map.end(); ++i) {
            if(i->t1== trajs[0] & x[t]>i->min & x[t]<i->max) {
                // add mapping   
            ntraj++; 
            trajs[ntraj-1]= i->t2; 
            }
        }
    }
}

//double CSplines::eval(int t, double x, traj_type traj) const
//{
//    moveTo(t);
//    return eval(x,traj);
//}

double CSplines::eval(double _x, traj_type traj) const
{
    if(traj < 0) return 1;

    hippo_Assert(traj >= 0 && traj < nFct, "traj out of bound");
    hippo_Assert(_x >= cpx[traj][1], "x coordinate smaller than control points");
    hippo_Assert(_x < cpx[traj][ncpx[traj]-2], "x coordinate larger than control points");

//    cerr << x << '\t';
    double result= op->op(fctCurr[traj]->eval(_x));
    hippo_Assert(isfinite(result), "result is infinite");
    return result;
}


void  CSplines::eval(double _x, traj_type _id, Grad &g) const 
{
    if(_id < 0) { g.markInvalid(nUpdate); return; }

    g.f= fctCurr[_id]->evalPartial(_x, g.df);
//    g.f= fctCurr[traj]->evalPartial(x, g.df, csegx[tCurr]);
    op->op(nUpdate, g);
}

double CSplines::evalGrad(double _x, traj_type traj, double diff[]) const
{
    if(traj < 0)  {
        for(int i=0; i<nUpdate; i++) diff[i]= 0;
        return 1;
    }

//    hippo_Print(t)
//    hippo_Print(csegx[t])

    hippo_Assert(traj >= 0, "trajectory cannot be negative");

//    traj_type tsegx= csegx[tCurr];

//    hippo_Assert(tsegx>= 0, "x segment is invalid. Shouldn't get here.");
//    hippo_Assert(x >= cpx[traj][tsegx]  && x <= cpx[traj][tsegx+1],
//        "inconsistency");
//    double result= fctCurr[traj]->evalPartial(x, diff, tsegx);
    double result= fctCurr[traj]->evalPartial(_x, diff);
    op->op(result, nUpdate, diff);
    return result;
}

double CSplines::getIntegral(int t, traj_type traj) const
{
    moveTo(t);
    return fctCurr[traj]->integral();
}

CSplines* CSplines::link_copy()
{
    CSplines *v= new CSplines;
    copyInit(v);
    v->store= store->alloc_new();
    v->store->init(startindex, endindex, v);
    return v;
}

void CSplines::copyInit(CSplines *v) const
{
    v->ncpx= ncpx;
    v->cpx= cpx;
    v->csegx= csegx;
    v->allocated.controls= false;

    v->initObj();
}


void CSplines::setMexInput(const mxArray *in)
{
    CSFunction::setMexInput(in);
    //    id= data->fieldID;
    int tmp;
    // cp is matlab cell array with column vectors
    mxArray *cellPtr, *mxtmp;

    cellPtr= mxGetField(in,0,"cpx");
    hippo_Assert(cellPtr, "no control points defined, errors from now on");
    if(mxIsCell(cellPtr)) {
        tmp = mxGetNumberOfElements(cellPtr);
        nFct= tmp;
        ncpx= new int[nFct];
        cpx= new double*[nFct];
        for (int i = 0; i < nFct; i++) {
            mxtmp = mxGetCell(cellPtr, i);
            ncpx[i] = mxGetN(mxtmp) * mxGetM(mxtmp);
            cpx[i]  = mxGetPr(mxtmp);
        }
    } else {
        nFct= 1;
        ncpx= new int;
        cpx= new double*;
        ncpx[0] = mxGetN(cellPtr) * mxGetM(cellPtr);
        cpx[0]  = mxGetPr(cellPtr);
    }
    allocated.controls= true;

    nTotal= 0;
    for (int i = 0; i < nFct; i++) {
        nTotal+= ncpx[i];
    }
}    

mxArray* CSplines::getMexOutput() 
{
    mxArray *out= CSFunction::getMexOutput();
    mxAddField(out, "name");
    mxAddField(out, "cpx");
//    mxAddField(out, "csegx");

    mxArray* mxtmp= mxCreateString("CSplines");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);

    double *dblPtr;
    mxArray *cellPtr;
    cellPtr= mxCreateCellMatrix(1,nFct);
    for (int i = 0; i < nFct; i++) {
        MX_CreateDoubleAssign(mxtmp, 1,ncpx[i], dblPtr);
        memcpy(dblPtr, cpx[i], ncpx[i] * sizeof(double));
        mxSetCell(cellPtr, i, mxtmp);
    }
    mxSetField(out, 0, "cpx", cellPtr);

//    MX_CreateDoubleAssign(mxtmp, 1, T, dblPtr);
//    for(int t=0; t < T; t++) dblPtr[t]=  csegx[t];
//    mxSetField(out, 0, "csegx", mxtmp);

    return out;
}
