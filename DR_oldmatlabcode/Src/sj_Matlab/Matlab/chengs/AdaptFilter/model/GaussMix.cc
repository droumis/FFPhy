/* $Id: GaussMix.cc,v 1.5 2008/08/24 19:47:04 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
//#include "AscentFilter.h"
#include "GaussMix.h"

/* ************************************************************
                          class GaussMix
   ************************************************************ */

using namespace AFilter;

static const double _2_pi= 0.63661977236758;

GaussMix::GaussMix(TData *_data)
{
    nMix= 0;
    comp= 0;
    alloc= false;
    T= 0;
    mult= 1;
    data= _data;
}

GaussMix::~GaussMix()
{
    dealloc();
}

void GaussMix::dealloc()
{
    if( alloc ) {
	if(comp) delete[] comp;
    }
    startindex= endindex= T= 0;
    alloc= false;
}

void GaussMix::init(int _startindex, int _endindex, double *x, double *y, traj_type *id)
{
    VDynamics::init(_startindex, _endindex, _T, _nid, _id);
    hippo_Message("empty");
}

void GaussMix::setInitialEst() const
{
    // initialize parameter values
    hippo_Message("no function needed at this point");
}

double GaussMix::eval(int t, double x, double y, double isi, traj_type traj)  const
{
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
    return mult*eval(t, x, y, traj);
}

double GaussMix::eval(int t, double x, double y, traj_type traj) const
{
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
    double lambda=0;
    for (int n= 0; n < nMix; n++) {
	lambda+= comp[n].eval(t, x, y, traj);
    }
    return lambda;
}

// double GaussMix::evalPartial(int t, double x, double y, 
// 				  traj_type traj, double partial[]) const
// {
//    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
//     double lambda=0;
//     for (int n= 0; n < nMix; n++) {
// 	lambda+= comp[n].evalPartial(t, x, y, traj, partial + (nCompParam*n));
//     }
//     return lambda;
// }

double GaussMix::evalGrad(int t, double x, double y, 
			    traj_type traj, double partial[]) const
{
    const double tiny= 1e-15;
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
    double lambda=0;
    for (int n= 0; n < nMix; n++) {
	lambda+= comp[n].evalGrad(t, x, y, traj, partial+(nCompParam*n));
    }

    if (lambda > tiny) {
	for (int i= 0; i < nCompParam*nMix; i++) partial[i]/= lambda;
    } else {
	hippo_Print(lambda);
	for (int i= 0; i < nCompParam*nMix; i++) {
	    partial[i]/= tiny;
	}
    }
    
    return lambda; 
}

void GaussMix::updateDynamics(int tStart, int tEnd, double diff[])
{
    for (int n= 0; n < nMix; n++) {
	comp[n].updateDynamics(tStart, tEnd, diff+(nCompParam*n));
    }
}

// TFilter*  GaussMix::allocFilter(char algo[]) 
// {
//     if (strcmp(algo,"SteepestAscent")==0) {
// 	return new AFilter::AscentFilter(data, this);
//     } else {
// 	MX_Error("Algorithm \'%s\' has not been implemented for AFilter::GaussMix", algo);
// 	return 0;
//     }
// };  

void  GaussMix::setMexInput(const mxArray *in) 
{
    VDynamics2d::setMexInput(in);
    mxArray *mxtmp= mxGetField(in,0,"GaussComp");
    if (mxtmp) {
	nMix= mxGetN(mxtmp)* mxGetM(mxtmp);
	comp= new Gauss[nMix](data);
	for (int n= 0; n < nMix; n++) {
	    comp[n].setMexInput(mxGetCell(mxtmp,n));
	}
	T= comp[0].getNTimesteps();
    } else {
	hippo_Message("no Gaussian components given, hope that's ok");
    }
    
    MX_FieldScalarDefault(in, "mult", mult, 1, double);

//    hippo_Print(T);
};

mxArray*  GaussMix::getMexOutput() 
{
    mxArray *out= VDynamics2d::getMexOutput();
    mxAddField(out, "name");
    mxAddField(out, "GaussComp");
    mxAddField(out, "mult");
    
    mxArray *mxtmp= mxCreateString("GaussMix");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);
    
    mxtmp= mxCreateCellMatrix(1, nMix);
    hippo_Assert(mxtmp, "couldn't allocate space for Gaussian components");
    for (int n= 0; n < nMix; n++) {
	mxSetCell(mxtmp, n, comp[n].getMexOutput());
    }
    mxSetField(out, 0, "GaussComp", mxtmp);
    
    mxtmp=mxCreateDoubleScalar(mult);
    mxSetField(out, 0, "mult", mxtmp);
    
    return out;
};

