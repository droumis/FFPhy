/* $Id: StatLimits.cc,v 1.5 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <string.h>

#include "StatLimits.h"

/* ************************************************************
                          class StatLimits
   ************************************************************ */

using namespace AFilter;

StatLimits::StatLimits(const mxArray *in)
{
    zeroOut();
    if (in) setMexInput(in);
}
StatLimits::StatLimits(StatLimits &)
{
    //   error;
}

StatLimits::~StatLimits()
{
    dealloc();
}

void StatLimits::zeroOut()
{
    valid= false;
    nana=0;
    nvar= nsteps= 0;
    xname= 0;
    varnames= 0;
    xmin= xmax= 0;
    varmin= varmax= 0;
    stepBegin= stepEnd= 0;
}

void StatLimits::dealloc() 
{
    if (nvar) delete[] nvar;
    if (nsteps) delete[] nsteps;
    if (xmin) {
        delete[] xmin;
        // the elements were allocated by matlab
    }
    if (xmax) {
        delete[] xmax;
        // the elements were allocated by matlab
    }
    if (varmin) {
        for (mwSize i= 0; i < nana; i++) delete[] varmin[i];
        delete[] varmin;
    }
    if (varmax) {
        for (mwSize i= 0; i < nana; i++) delete[] varmax[i];
        delete[] varmax;
    }
    if (xname) { delete[] xname; }
    if (varnames) {
        for (mwSize i= 0; i < nana; i++) delete[] varnames[i];
        delete[] varnames;
    }
    if (stepBegin) delete[] stepBegin;
    if (stepEnd) delete[] stepEnd;

    zeroOut();
}

int StatLimits::name2index(const char varname[], int ana) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    for (int i= 0; i< nvar[ana]; i++) {
        if (strcmp(varname, varnames[ana][i])== 0) return i;
    }
    return -1;
}

bool StatLimits::isConstraint(const char varname[], int ana) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    return (name2index(varname, ana)!= -1);
}

double StatLimits::getXMin(mwSize ana, int step) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    hippo_Assert(ana < nana, "analysis requested out of bounds");
    hippo_Assert(step < nsteps[ana], "timestep requested out of bounds");

    return xmin[ana][step];
}

double StatLimits::getXMax(mwSize ana, int step) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    hippo_Assert(ana < nana, "analysis requested out of bounds");
    hippo_Assert(step < nsteps[ana], "timestep requested out of bounds");
    return xmax[ana][step];
}

double StatLimits::getVarMin(const char varname[], int ana) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    int ind= name2index(varname, ana);
    if (ind == -1) {
        hippo_Error("Bounds not defined for variable ", varname);
    }
    return varmin[ana][ind];
}

double StatLimits::getVarMax(const char varname[], int ana) const
{
    hippo_Assert(valid, "initialize StatLimit object properly");
    int ind= name2index(varname, ana);
    hippo_Assert(ind != -1, "requested variable not found");
    return varmax[ana][ind];
}

void StatLimits::setMexInput(const mxArray *in)
{
    dealloc(); // in case object was allocated before
    hippo_Assert(!mxIsEmpty(in), "limits object is empty");
    hippo_Assert(mxIsStruct(in), "limits object is not a structure");
    if(!mxIsStruct(in)) {
        valid= false;
        return;
    }
    
    const mxArray *mxa;
    const mxArray *mxx;
    MX_GetField(in, mxa, "a");
    hippo_Assert(mxa, "Field 'a' not defined in limits");
    MX_GetField(in, mxx, "x");
    hippo_Assert(mxx, "Field 'x' not defined in limits");

    nana= mxGetNumberOfElements(mxa);
    hippo_Assert(mxGetNumberOfElements(mxx)==nana, "a and x have to have same length")
    nvar= new int[nana];
    nsteps= new int[nana];
    xmin= new double *[nana];
    xmax= new double *[nana];
    varmin= new double *[nana];
    varmax= new double *[nana];
    varnames= new char const **[nana];
    xname= new char const *[nana];

    const mxArray *tmpana, *field;
    // read in the contraints
    for(mwSize ana= 0; ana < nana; ana++) {
        tmpana= mxGetCell(mxa, ana);
        nvar[ana]= mxGetNumberOfFields(tmpana);
        varmin[ana]= new double[nvar[ana]];
        varmax[ana]= new double[nvar[ana]];
        varnames[ana]= new char const *[nvar[ana]];
        for (int i= 0; i < nvar[ana]; i++) {
            varnames[ana][i]= mxGetFieldNameByNumber(tmpana, i);
            field= mxGetFieldByNumber(tmpana, 0, i);
            if (mxGetNumberOfElements(field)== 1) {
                varmin[ana][i]= mxGetNaN();
                varmax[ana][i]= mxGetScalar(field);
            } else {
                varmin[ana][i]= mxGetPr(field)[0];
                varmax[ana][i]= mxGetPr(field)[1];
            }
        }
    }

    // read in the independent variable
    for(mwSize ana= 0; ana < nana; ana++) {
        tmpana= mxGetCell(mxx, ana);
        xname[ana]= mxGetFieldNameByNumber(tmpana, 0);
//        cerr << ana << ": " << xname[ana] << endl;
        field= mxGetFieldByNumber(tmpana, 0, 0);
        int m= mxGetM(field);
        int n= mxGetN(field);
        if (m==1 || n==1) {
            xmin[ana]= 0;
            xmax[ana]= mxGetPr(field);
            nsteps[ana]= n > m ? n : m;
        } else if (n==2) {
            nsteps[ana]= n;
            xmin[ana]= mxGetPr(field);
            xmax[ana]= mxGetPr(field)+n;
        } else {
            hippo_ErrorMsg("Time array has incorrect format. Must be vector or 2xn matrix.");
        }
    }

    valid= true;
}

mxArray* StatLimits::getMexOutput()
{
    hippo_Empty;
    return 0;
}

