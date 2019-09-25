/* $Id: CSFunction.cc,v 1.8 2008/09/01 18:24:22 chengs Exp $
   
   Sen Cheng, Thu Jun  8 14:00:46 PDT 2006
*/

//#include <stdlib.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSFunction.h"
#include "CSIncrement.h"
#include "CSIncRectify.h"
#include "CSCompress.h"

/* ************************************************************
                          class CSFunction
   ************************************************************ */

using namespace AFilter;

CSFunction::CSFunction() : map()
{
//    fct=  0;
    op= 0;
    allocOp= false;
    store= 0;
    inputCompressed= outputCompress= rectify= false;
    outputInterval= 0;
    nPerFct= 0;
    mapping= false;
}

CSFunction::~CSFunction()
{
    if(store) delete store;
    if(allocOp && op) delete op;
//    hippo_Print(nPerFct);
//    if(nPerFct) delete[] nPerFct; // is allocated in child class
}

void CSFunction::init(int _startindex, int _endindex, int _T,
        int _nid, traj_type *_id) 
{

    VDynamics::init(_startindex, _endindex, _T, _nid, _id);


    if(!store) store= rectify? new CSIncRectify: new CSIncrement;
    hippo_Assert(nPerFct, "nPerFct not allocated in child class!");
    hippo_Assert(fct.size(), "fct object not allocated in child class!");
    store->init(_startindex, _endindex, this);

    if(!op) {
        op= (rectify)? VOp::newOperator("Rectify") : VOp::newOperator("Id");
        allocOp= true;
    }
}   

bool CSFunction::getParam(int t, double z[], int a[]) const
{
    moveTo(t);

    traj_type *trajs= new traj_type[getNParam()], tid;
    int nt;
    getActiveIndices(t, nt, trajs, a);
    tid= trajs[0];
    if(tid >= 0)  {
        fct[tid]->getParamPartial(z, nUpdate, a);
        if(tid>0) {
            int offset= 0;
            for(int i=0; i<tid; ++i) offset+= nPerFct[i];
            for(int i= 0; i<nUpdate; ++i) a[i]+= offset;
        }
        return true;
    } else {
        for(int i= 0; i < nUpdate; i++) a[i]= -1;
        return false;
    }
}

void CSFunction::getAllParam(int t, double z[]) const { store->getAllParam(t,z); };

void CSFunction::setAllParam(int t, double z[]) { store->setAllParam(t,z); }

void CSFunction::setAllParam(int t, double z) { store->setAllParam(t,z); }

void CSFunction::setValues(int t, double z) { 
    store->setAllParam(t,op->inv(z)); 
}

double CSFunction::getMinZ() const { return store->getMin(op); } ;
double CSFunction::getMaxZ() const { return store->getMax(op); } ;

void CSFunction::updateDynamics(int tStart, int tEnd, double diff[]) {
    store->updateEst(tStart, tEnd, diff);
}

void CSFunction::moveTo(int t) const
{
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
    store->moveTo(t);
    tCurr= t;
}

void CSFunction::allocForConverge()
{
    store->allocForConverge();
}

void CSFunction::copyForConverge()
{
    store->copyForConverge();
}

bool CSFunction::hasConverged()
{
    double maxdiff, reldiff;
    store->maxDiff(maxdiff, reldiff);
    cerr <<  "  maxdiff= " << maxdiff << "\t"
        <<  "reldiff= " << reldiff << "\n";

    return maxdiff < conv && reldiff < convper;
}

void CSFunction::averageWithCopy() { store->averageWithCopy(); }

void CSFunction::copyInit(CSFunction *v) const
{
//    v->tIni= tIni;

    v->nFct= nFct;
    v->nUpdate= nUpdate;
    v->nPerFct= nPerFct;
    v->nTotal= nTotal;
    v->id= id;
    v->op= op;
    v->allocOp= false;
//        v->op= VOp::newOperator("Id"); //@@
//        v->allocOp= true;

    v->inputCompressed= inputCompressed;
    v->outputCompress= outputCompress;
    v->outputInterval= outputInterval;
    v->conv= conv;
    v->convper= convper;
    v->rectify= rectify;
}


void CSFunction::setMexInput(const mxArray *in)
{
    VDynamics::setMexInput(in);
    mxArray *mxtmp;

    mxtmp= mxGetField(in, 0, "rectify");
    if(mxtmp) {
        MX_FieldScalar(in, "rectify", rectify, bool);
    } else
        rectify= false;

    mxtmp= mxGetField(in, 0, "initAll");
    if(mxtmp) {
        MX_FieldScalar(in, "initAll", initAll, double);
    }

    MX_FieldScalarDefault(in, "conv", conv, 2, double);
    MX_FieldScalarDefault(in, "convper", convper, .1, double);

    mxtmp= mxGetField(in, 0, "operator");
    if(mxtmp) {
        op= VOp::newOperator(mxtmp);
        allocOp= true;
    } else { hippo_Warning("operator not defined"); }

    if(op) {
        if(rectify && strcmp(op->getOpName(),"Id")==0) 
            hippo_ErrorMsg("rectify is inconsistent with OpId!");
        if(strcmp(op->getOpName(),"Rectify")==0) rectify= true;
    }

    MX_FieldScalarDefault(in, "outputCompress", outputCompress, false, bool);
    MX_FieldScalarDefault(in, "inputCompressed", inputCompressed, false, bool);
    if(!outputCompress && inputCompressed) hippo_ErrorMsg("do not know how to uncompress");
    if(outputCompress || inputCompressed) 
        MX_FieldScalarDefault(in, "outputInterval", outputInterval, 1, int);


    map.clear();
    mxtmp= mxGetField(in,0,"map");
    if(mxtmp) {
        hippo_Print(mxGetM(mxtmp));
        hippo_Assert(mxGetM(mxtmp)==4, "parameter 'map' must be 4-by-x array");
        int n= mxGetN(mxtmp);
        double *d= mxGetPr(mxtmp);
        mapItem item;
        for(int in=0; in<n; in++) {
            item.t1= (traj_type) *d++; 
            item.min= *d++; 
            item.max= *d++; 
            item.t2= (traj_type) *d++; 

            map.push_back(item);
        }
        cerr << "mapping control points\n";
        vector<mapItem>::iterator iter;
        for(iter= map.begin(); iter!= map.end(); ++iter) {
            cerr << "  " << iter->min
                 << "cm, " << iter->max
                << "cm: map trajectories " << iter->t1
                 << " -> " << iter->t2 << "\n";
        }
    }
    mapping= map.empty() ? false : true;

//    hippo_Print(outputCompress);
//    hippo_Print(inputCompressed);

    if(inputCompressed)  {
        store= new CSCompress;
    //@@ need rectified compressed as well!
    } else if(rectify) {
        store= new CSIncRectify;
    } else  {
        store= new CSIncrement;
    }

    mxtmp= mxGetField(in, 0, "parameters");
    if(mxtmp) store->setMexInput(mxtmp);

}    

mxArray *CSFunction::getMexOutput() 
{
    mxArray *out= VDynamics::getMexOutput();
    mxArray *mxtmp;
    mxAddField(out, "name");
    mxSetField(out, 0, "name", mxCreateString(getName()));

    mxAddField(out, "rectify");
    mxSetField(out, 0, "rectify", mxCreateLogicalScalar(rectify));

    mxAddField(out, "operator");
    mxSetField(out, 0, "operator", op->getMexOutput());

    bool compresstmp= inputCompressed;
     
    if(inputCompressed == outputCompress) {
        mxtmp= store->getMexOutput();
        compresstmp= inputCompressed;
    } else if(!inputCompressed && outputCompress) {
        CSCompress s(*store);
        mxtmp= s.getMexOutput();
        mxAddField(out, "outputInterval");
        mxSetField(out, 0, "outputInterval", mxCreateScalarDouble(outputInterval));
        compresstmp= true;
    } else {
        mxtmp= 0;
        hippo_ErrorMsg("cannot uncompress");
    }
//    hippo_Print(compresstmp);
    mxAddField(out, "inputCompressed");
    mxSetField(out, 0, "inputCompressed", mxCreateLogicalScalar(compresstmp));

    mxAddField(out, "outputCompress");
    mxSetField(out, 0, "outputCompress", mxCreateLogicalScalar(outputCompress));

    mxAddField(out, "parameters");
    mxSetField(out, 0, "parameters", mxtmp);

    mxAddField(out, "conv");
    mxSetField(out, 0, "conv", mxCreateScalarDouble(conv));

    mxAddField(out, "convper");
    mxSetField(out, 0, "convper", mxCreateScalarDouble(convper));

    if(initAll.valid) {
        mxAddField(out, "initAll");
        mxSetField(out, 0, "initAll", mxCreateScalarDouble(initAll));
    }

    hippo_Print(store->getMin(op));
    hippo_Print(store->getMax(op));

    // save mapping of control pts.
    double *dblptr;
    if(mapping) {
        MX_CreateDoubleAssign(mxtmp,4,map.size(),dblptr)
        vector<mapItem>::iterator iter;
        for(iter= map.begin(); iter!= map.end(); ++iter) {
            *dblptr++= iter->t1;
            *dblptr++= iter->min;
            *dblptr++= iter->max;
            *dblptr++= iter->t2;
        }
        mxAddField(out, "map");
        mxSetField(out, 0, "map", mxtmp);
    }

    return out;
}
