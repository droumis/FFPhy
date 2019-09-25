/* $Id: IsiModel.cc,v 1.3 2008/10/23 21:24:10 chengs Exp $

   Sen Cheng, 2004/../..

   program description
   */

#include <stdlib.h>
#include <string.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "modelFactory.h" 
#include "IsiModel.h"


/* ************************************************************
   class IsiModel
 ************************************************************ */
using namespace AFilter;
using namespace AFilter;


IsiModel::IsiModel(TData *_data)
: AdaptModel(_data), xp(0), isi(0), nxp(0), nisi(0), totalnxp(0), isiID(0)
{
    allocated.isiID= false;
}

IsiModel::~IsiModel()
{
    if(xp) delete xp;
    if(isi) delete isi;
    if(allocated.isiID && isiID) delete[] isiID;
}

void IsiModel::init()
{

    hippo_Assert(data, "no data object defined");
    hippo_Assert(xp, "no spatial model object defined");
    hippo_Assert(isi, "no isi model object defined");

    int n= 0;
//    while(n < data->getNTimesteps() && (!mxIsFinite(data->posarray[n]) || !mxIsFinite(data->phasearray[n]) || !mxIsFinite(data->isi[n]) ))
    while(n < data->getNTimesteps() && (!mxIsFinite(data->posarray[n]) ))
        n++;
    startindex= n;
//    while(n < data->getNTimesteps() && mxIsFinite(data->posarray[n]) && 
//            mxIsFinite(data->phasearray[n]) && mxIsFinite(data->isi[n]) ) 
    while(n < data->getNTimesteps() && mxIsFinite(data->posarray[n])) 
        n++;
    endindex= n-1;
    if((double)(endindex-startindex)/(double)data->getNTimesteps() < 0.95)
        hippo_Error("too many invalid timesteps", 
                    data->getNTimesteps()-(endindex-startindex));

//    xp->init(startindex, endindex, n, data->ntraj, data->fieldID); //@@
    this->init_spatial(n);
    nxp= xp->getNParam();
    totalnxp= xp->getNAllParam();

//    hippo_Print(nxp);
//    hippo_Print(totalnxp);

    isiID= new traj_type[n];
    allocated.isiID= true;
    for (int t=0; t < n; t++) {
        if (data->isi[t] > 0 && data->isi[t] < max_isi) {
            isiID[t]=0;
        } else {
            isiID[t] =-1;
        }
    }
    isi->init(startindex, endindex, n, data->isi, 1, isiID);
    nisi= isi->getNParam();

//    hippo_Print(nisi);

    isi->setMinX(0.001); //@@ temp. hack!  should be min_isi in model
    isi->setMaxX(max_isi);

    cerr << "Model initialized:";
    cerr << "  valid timesteps= " << startindex << " to " << endindex << endl;
}

VDynamics1d* IsiModel::get1dSpatialModel() const 
{
    // first ensure that the model is compatible with 1-d
    hippo_Assert(xp->getDim()==1 | xp->getDim()==0, "PosPhase_Isi requires 2-dim spatial model.");

    return dynamic_cast<VDynamics1d*>(xp);
}

VDynamics2d* IsiModel::get2dSpatialModel() const 
{
    // first ensure that the model is compatible with 2-d
    hippo_Assert(xp->getDim()==2 | xp->getDim()==0, "PosPhase_Isi requires 2-dim spatial model.");

    return dynamic_cast<VDynamics2d*>(xp);
}

bool IsiModel::getParam(int t, double z[], int a[]) const
{
    bool txp= xp->getParam(t, z, a);
    bool tisi= isi->getParam(t, z+nxp, a+nxp);
    if(tisi) {
        for(int j=0;  j< nisi; ++j)  a[nxp+j]+= totalnxp;
    }
    return txp>=0 || tisi>=0;
}

void IsiModel::getAllParam(int t, double _x[]) const
{
    xp->getAllParam(t, _x);
//    hippo_Print(totalnxp)
    isi->getAllParam(t, _x+ totalnxp);
}

void IsiModel::setParam(int t, double z[])
{
    hippo_Empty;
}

void IsiModel::setInitialEst(int t)
{ 
    double r= data->getMeanRate();
    xp->setValues(t,r);  
    isi->setValues(t,1);  
}

double IsiModel::eval(int t) const
{
    double xpval= xp->evalAtTime(t);
    double isival= isi->evalAtTime(t);
//    hippo_Print(xpval);
//    hippo_Print(isi);
//    hippo_Print(isival);
//    hippo_Print(isi->eval(0, 0));
//    hippo_Print(data->isi[t]);
//    hippo_Print(isiID[t]);
//    hippo_Exit;
    return xpval*isival;
}

double IsiModel::evalAtIsi(int t, double _isi){
    double xpval= xp->evalAtTime(t);
    data->isi[t]= _isi;
    isiID[t]= (_isi < max_isi) ? 0 : -1;
    double isival= isi->evalAtTime(t);
    return xpval*isival;
}

void IsiModel::evalNoMove(int t, Grad &g) const
{
    g.reset();
    xp->evalNoMove(t, g);
    double f= g.f, logf= g.logf;
    g+= nxp; 
    isi->evalNoMove(t,g); 
    g-= nxp;
    g.f*= f; g.logf+= logf;
}

void IsiModel::eval(int t, Grad &g) const 
{
    hippo_ErrorMsg("BUG needs FIX");
    g.reset();
    xp->evalAtTime(t, g);
    double f= g.f, logf= g.logf;
    g+= nxp; 
    isi->evalAtTime(t,g); 
    g-= nxp;
    g.f*= f; g.logf+= logf;
}

double IsiModel::evalGrad(int t, double diff[]) const
{
    double xpval= xp->evalGradAtTime(t, diff);
    double isival= isi->evalGradAtTime(t, diff+nxp);
    return xpval*isival;
}

double IsiModel::evalNoMove(int t, double diff[]) const
{
    double xpval= xp->evalNoMove(t, diff);
//    hippo_Print(xpval);
//    hippo_Print(diff[0]);
    double isival= isi->evalNoMove(t, diff+nxp);
//    hippo_Print(isival);
//    hippo_Print(*(diff+nxp));
//    exit(0);
    return xpval*isival;
}

//double IsiModel::evalLogGrad(int t, double diff[]) const
//{
//    double xpval= xp->evalLogGrad(t, diff);
//    double isival= isi->evalLogGrad(data->isi[t],isiID[t],diff+nxp);
//    return xpval+isival;
//}

void IsiModel::updateModel(int tStart, int tEnd, double diff[], int fixing)
{
    if(fixing != 0) xp->updateDynamics(tStart, tEnd, diff);
    if(fixing != 1) isi->updateDynamics(tStart, tEnd, diff+nxp);
}

void IsiModel::allocForConvergence() 
{
    convCopy= true;
    isi->allocForConverge();
    xp->allocForConverge();
}

void IsiModel::copyForConverge() 
{
    hippo_Assert(convCopy, "must allocate copies before copying for convergence");
    xp->copyForConverge();
    isi->copyForConverge();
}

bool IsiModel::hasConverged() 
{
    hippo_Assert(convCopy, "must allocate copies before testing convergence");
    return  xp->hasConverged() && isi->hasConverged();
}

IsiModel* IsiModel::link_copy()
{
//    hippo_Mark;
    AdaptModel *am= alloc();
    IsiModel *im= reinterpret_cast<IsiModel*>(am);
    im->data= data;
    im->xp= xp->link_copy();
    im->isi= isi->link_copy();
    im->totalnxp= totalnxp;
    im->nxp= nxp;
    im->nisi= nisi;
    im->max_isi= max_isi;
    im->isiID= isiID;
    im->allocated.isiID= false;
    im->startindex= startindex;
    im->endindex= endindex;
    return im;
}

void IsiModel::setMexInput(const mxArray *in)
{
    mxArray *modelPtr;

    MX_FieldScalar(in, "max_isi", max_isi, double);
//    hippo_Print(max_isi);

    modelPtr= mxGetField(in,0,"spatial");
    hippo_Assert(modelPtr, "spatial model not defined in matlab input");
    xp= allocDynamics(modelPtr, data);
//    hippo_Print(xp);
//    hippo_Print(xp->getName());

//    hippo_Mark;
    modelPtr= mxGetField(in,0,"isi");
    hippo_Assert(modelPtr, "isi model not defined in matlab input");
    isi= allocDynamics1d(modelPtr, data);
//    hippo_Print(isi);
//    hippo_Mark;
}

mxArray* IsiModel::getMexOutput() 
{
    mxArray *out= AdaptModel::getMexOutput();

    mxAddField(out, "name");
    mxArray *mxtmp= mxCreateString(getName());
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);

    mxAddField(out, "max_isi");
    mxtmp= mxCreateScalarDouble(max_isi);
    hippo_Assert(mxtmp, "couldn't allocate space for variable max_isi");
    mxSetField(out, 0, "max_isi", mxtmp);

    mxAddField(out, "spatial");
    mxSetField(out, 0, "spatial", xp->getMexOutput());
    mxAddField(out, "isi");
    mxSetField(out, 0, "isi", isi->getMexOutput());

    return out;
}
