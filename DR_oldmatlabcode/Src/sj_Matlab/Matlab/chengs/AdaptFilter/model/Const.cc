/* $Id: Const.cc,v 1.11 2008/10/23 21:24:10 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
#include "Const.h"

/* ************************************************************
                          class Const
   ************************************************************ */

using namespace AFilter;

Const::Const(Const &c)
//    : VDynamics1d(c), VDynamics2d(c)
    : VDynamics1d(c)
{
    constant= c.constant;
    conv= c.conv;
    convper= c.convper;
}

Const::Const(TData *_data) : AdaptModel(_data) 
{ 
    nFct= 1; constant=1; 
};

Const::~Const()
{
}


// Initialize as model for conditional intensity
void Const::init() 
{ 
    hippo_Assert(AdaptModel::data, "no data object defined");
    int n= 0;
    while(n < AdaptModel::data->getNTimesteps() && (!mxIsFinite(AdaptModel::data->posarray[n]) || !mxIsFinite(AdaptModel::data->phasearray[n]) || !mxIsFinite(AdaptModel::data->isi[n]) ))
        n++;
    AdaptModel::startindex= n;
    while(n < AdaptModel::data->getNTimesteps() && mxIsFinite(AdaptModel::data->posarray[n]) && 
            mxIsFinite(AdaptModel::data->phasearray[n]) && mxIsFinite(AdaptModel::data->isi[n]) ) 
        n++;
    AdaptModel::endindex= n-1;
    if((double)(AdaptModel::endindex-AdaptModel::startindex)/(double)AdaptModel::data->getNTimesteps() < 0.95)
        hippo_Error("too many invalid timesteps", 
                    AdaptModel::data->getNTimesteps()-(AdaptModel::endindex-AdaptModel::startindex));
    T= AdaptModel::endindex-AdaptModel::startindex+1;
};

// Initialize as 1-d function 
void Const::init(int _startindex, int _endindex, int _T, double *_x, 
                 int _nFct, traj_type *_id) 
{
    hippo_Mark;
    nFct= _nFct;
    VDynamics1d::init(_startindex, _endindex, _T, _x, _nFct, _id);
//    hippo_Print(this);
};

// Initialize as 2-d function 
void Const::init(int _startindex, int _endindex, int _T, double *_x, double *y,
                 int _nFct, traj_type *_id) {
    nFct= _nFct;
    VDynamics2d::init(_startindex, _endindex, _T, _x, y,_nFct, _id);
};

void Const::reestimate(AdaptModel *v)
{
    hippo_Empty;
}

void Const::setMexInput(const mxArray *in) 
{
//    VDynamics2d::setMexInput(in);
    VDynamics1d::setMexInput(in);

    MX_FieldScalarDefault(in, "constant", constant, -1, double);

    hippo_Assert(constant>0, "constant was not assigned in input");

    MX_FieldScalarDefault(in, "conv", conv, 2, double);
    MX_FieldScalarDefault(in, "convper", convper, .1, double);
//    hippo_Print(this);
//    hippo_Print(constant);
};

mxArray*  Const::getMexOutput() 
{
//    mxArray *out= VDynamics2d::getMexOutput();
    mxArray *out= VDynamics1d::getMexOutput();
    mxAddField(out, "name");
    
    mxArray *mxtmp;
    
    mxtmp= mxCreateString("Const");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);
    

    mxtmp= mxCreateScalarDouble(constant);
    hippo_Assert(mxtmp, "couldn't allocate space for variable conv");
    mxAddField(out, "constant");
    mxSetField(out, 0, "constant", mxtmp);

    mxtmp= mxCreateScalarDouble(conv);
    hippo_Assert(mxtmp, "couldn't allocate space for variable conv");
    mxAddField(out, "conv");
    mxSetField(out, 0, "conv", mxtmp);

    mxtmp= mxCreateScalarDouble(convper);
    hippo_Assert(mxtmp, "couldn't allocate space for variable conv");
    mxAddField(out, "convper");
    mxSetField(out, 0, "convper", mxtmp);

    hippo_Print(constant);

    return out;
};

