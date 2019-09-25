/* $Id: VStat.cc,v 1.4 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, Fri Sep  3 12:58:42 PDT 2004
   
   program description
*/

#include "VStat.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class VStat
   ************************************************************ */

using namespace AFilter;

VStat::VStat(const TData *_data, const AdaptModel *m, const StatLimits *_lim)
    : data(_data), model(m), lim(_lim)
{
    stat= 0;
    nana= 0;
    ndim= 0;
    nsteps= 0;
    label= 0;
}

VStat::~VStat()
{
    dealloc();
}

void VStat::dealloc()
{
    if (stat) {
	for (int a= 0; a < nana; a++) delete[] stat[a];
	delete[] stat;    
    }
    if (nsteps) delete[] nsteps;
    if (label) mxFree(label);
}

void VStat::init(int a, int d, int *s)
{
    dealloc();

    nana= a;
    ndim= d;
    nsteps= new int[nana];
    stat= new double*[nana];

    for (int a= 0; a < nana; a++) {
        nsteps[a]= s? s[a] : 1;
        stat[a]= new double[ndim*nsteps[a]];
//    hippo_Print(nsteps[a]);
    }

}


double VStat::get(int a, int dim, int step) const
{
    hippo_Assert(a >= 0 && a < nana, "analysis number out of bounds");
    hippo_Assert(dim >= 0 && dim < ndim, "dimension number out of bounds");
    hippo_Assert(step >= 0 && step < nsteps[a], "step number out of bounds");
    hippo_Assert(stat, "object in empty");

    return stat[a][dim*nsteps[a]+step];
}

void VStat::set(int a, int dim, int step, double x)
{
//    printf("ana= %d, dim= %d, steps= %d, x= %f\n", a,dim,steps,x);
//    cout << nsteps;
//    printf(", nsteps= %d", nsteps[a]); 
    hippo_Assert(a >= 0 && a < nana, "analysis number out of bounds");
    hippo_Assert(dim >= 0 && dim < ndim, "dimension number out of bounds");
    hippo_Assert(step >= 0 && step < nsteps[a], "step number out of bounds");
    hippo_Assert(stat, "object in empty");

    stat[a][dim*nsteps[a]+step]= x;
}

void VStat::setMexInput(const mxArray *in)
{
    mxArray *mxtmp= mxGetField(in, 0, "label");
    if(mxtmp) {
        MX_StringAllocAssign(mxtmp, label);
    } else if (mxtmp) {
        mxFree(label);
        label= 0;
    }
}

mxArray* VStat::getMexOutput()
{
    mxArray *out= mxCreateCellMatrix(nana, 1);
    hippo_Assert(out, "could allocate cell matrix for output");
    mxArray *mxtmp;
    double *dblPtr;
    for (int a=0; a < nana; a++) {
        MX_CreateDoubleAssign(mxtmp, nsteps[a], ndim, dblPtr);
        memcpy(dblPtr, stat[a], ndim*nsteps[a]* sizeof(double));
        mxSetCell(out, a, mxtmp);
    }
    return out;
}

mxArray* AFilter::collectMexOutput(int n, VStat *stats[])
{
    const char *fields[n];
    for (int i= 0; i < n; i++) {
        if (stats[i]->label)
            fields[i]= stats[i]->label;
        else
            fields[i]= stats[i]->getName();
    }
    mxArray *out= mxCreateStructMatrix(1,1, n,  fields);
    for (int i= 0; i < n; i++) {
        mxSetFieldByNumber(out, 0, i, stats[i]->getMexOutput());
    }
    return out;
}


