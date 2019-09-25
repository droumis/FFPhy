/* $Id: Uniform.cc,v 1.3 2006/07/14 21:37:06 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
#include "Uniform.h"

/* ************************************************************
                          class Uniform
   ************************************************************ */

using namespace AFilter;

Uniform::Uniform()
    : T(0), a(0)
{
}

Uniform::Uniform(Uniform &u)
{
    hippo_Empty;
}

Uniform::~Uniform()
{
    if(a) delete[] a;
}

void  Uniform::setMexInput(const mxArray *in) 
{
    VDynamics2d::setMexInput(in);
    MX_FieldAssign(in, "a", a, T);
};

mxArray*  Uniform::getMexOutput() 
{
    mxArray *out= VDynamics2d::getMexOutput();
    mxAddField(out, "name");
    mxAddField(out, "a");

    mxArray *mxtmp;
    double *ptr;
    
    mxtmp= mxCreateString("Uniform");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, a, T* sizeof(double));
    mxSetField(out, 0, "a", mxtmp);

    return out;
};

