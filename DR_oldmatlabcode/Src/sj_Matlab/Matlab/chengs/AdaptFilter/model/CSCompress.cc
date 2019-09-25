/* $Id: CSCompress.cc,v 1.5 2008/08/24 19:47:03 chengs Exp $
   
   Sen Cheng, Wed Oct  6 08:39:54 PDT 2004
   implementation of class CSCompress
*/

#include <stdlib.h>
#include <string.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSCompress.h"

/* ************************************************************
                          class CSCompress
   ************************************************************ */

using namespace AFilter;

//#define OTYPE double
//#define OTYPE_ID mxDOUBLE_CLASS
//#define OTYPE float
//#define OTYPE_ID mxSINGLE_CLASS
#define OTYPE short unsigned int
#define OTYPE_ID mxUINT16_CLASS
static const OTYPE maxOTYPE= 65535;


CSCompress::CSCompress()
{
    // really lame check, need to do better
    mxAssert(sizeof(OTYPE)== 2, "ohoh, storage size for variables don't match");

    outputInterval= 0;
    nInterval= 0;
    tParam= 0; tParamAllocated= 0;
}


CSCompress::CSCompress(const CSStorage &cs)
    : CSStorage()
{
    // really lame check, need to do better
    mxAssert(sizeof(OTYPE)== 2, "ohoh, storage size for variables don't match");
//    mxAssert(cs.startindex== cs.tIni, "initial time should be set to startindex");

    tParam= 0; tParamAllocated= 0;
    param= 0;

    startindex= cs.startindex;
    endindex= cs.endindex;
    T= cs.T;

    super= cs.super;
    nFct= cs.nFct;
    nUpdate= cs.nUpdate;
    nPerFct= cs.nPerFct;
    nTotal= cs.nTotal;

    outputInterval= super->getOutputInterval();
    nInterval= (int) floor((endindex-startindex)/outputInterval)+1;
//    hippo_Print(outputInterval);
//    hippo_Print(nInterval);
//    hippo_Print(nTotal);
//    hippo_Print(super->getName());

	hippo_CAlloc(tParam, nInterval*nTotal, double);
    tParamAllocated= true;

    // evaluate cardinal splines function at requested intervals 
    double *tmp= tParam;
    for (int t= startindex; t<= endindex; t+= outputInterval) {
        cs.getAllParam(t, tmp);
        tmp+= nTotal;
    }
    calcConversion();

}

CSCompress::~CSCompress()
{
    if(tParamAllocated && tParam) free(tParam);
}

void CSCompress::init(int _startindex, int _endindex, CSFunction *f)
{
    CSStorage::init(_startindex, _endindex, f);

    hippo_Assert(nInterval== floor((endindex-startindex)/outputInterval)+1, 
            "size of input parameter vector inconsistent");

}   

void CSCompress::moveTo(int t) const
{
    hippo_Assert(T != -1, "Parameters was not properly initialized");

    int index= getTimestep(t);
    if(index== tCurr) return;

    updateParam(index);
    tCurr= index;
}

void CSCompress::calcMinMax(const VOp *op) const
{
    min= DBL_MAX; max= -DBL_MAX;
    double x;
    for(int i=0; i<nInterval*nTotal; i++) {
        x= op->op(tParam[i]); 
        if(x > max) max= x;
        if(x < min) min= x;
    }
    min-= tiny;
    max+= tiny;
}

void CSCompress::calcConversion() const
{
    convertMin= DBL_MAX;
    convertMax= -DBL_MAX;
    for(int i=0; i<nInterval*nTotal; i++) {
        if(tParam[i] < convertMin) convertMin= tParam[i];
        if(tParam[i] > convertMax) convertMax= tParam[i];
    }
    convertMax+= tiny;
    convertMin-= tiny;
    convertFac= (fabs(convertMax-convertMin) < 1) ? (maxOTYPE-1) : (maxOTYPE-1)/(convertMax-convertMin);
}

void CSCompress::getAllParam(int t, double z[]) const
{
    double *tmp= tParam + getTimestep(t)*nTotal;
    for (traj_type i= 0; i < nFct; i++) {
        memcpy(z, tmp, nPerFct[i]*sizeof(double));
        tmp+=nPerFct[i];
        z+=nPerFct[i];
    }
}

void CSCompress::setAllParam(int t, double z[])
{
    int index= getTimestep(t);
    double *tmp= tParam + index*nTotal;
    for (traj_type i= 0; i < nFct; i++) {
        memcpy(tmp, z, nPerFct[i]*sizeof(double));
        tmp+=nPerFct[i];
        z+=nPerFct[i];
    }
    if(index== tCurr) updateParam(index);
}

void CSCompress::setAllParam(int t, double z)
{
    int index= getTimestep(t);
    double *tmp= tParam + index*nTotal;
    for(int i=0; i < nTotal; i++) tmp[i]= z;
    if(index== tCurr) updateParam(index);
}


//void CSCompress::updateDynamics(int tStart, int tEnd, double diff[])
//{
//    hippo_Empty;
//}

//double CSCompress::getIntegral(int timestep, traj_type traj) const
//{
//    moveTo(timestep);
//    return obj[traj]->integral();
//}

void CSCompress::setMexInput(const mxArray *in)
{
    mxArray *mxtmp;
    MX_FieldScalar(in, "outputInterval", outputInterval, int);

    MX_FieldScalar(in, "convertMin", convertMin, double);
    MX_FieldScalar(in, "convertFac", convertFac, double);

    if(tParamAllocated && tParam) {
        free(tParam);
        tParamAllocated= false;
    }
    mxtmp= mxGetField(in, 0, "param");
    if(mxtmp) {
        OTYPE *inptr= (OTYPE*) mxGetData(mxtmp);
        nTotal= mxGetM(mxtmp);
        nInterval= mxGetN(mxtmp);
        hippo_CAlloc(tParam,nInterval*nTotal,double);
        tParamAllocated= true;
        double *tmp= tParam;
        for (int i= 0; i < nInterval*nTotal; i++) {
            *tmp= (*inptr) / convertFac + convertMin;
            inptr++;
            tmp++;
        }
    }
    convertMax= -DBL_MAX;
    for(int i=0; i<nInterval*nTotal; i++) 
        if(tParam[i] > convertMax) convertMax= tParam[i];
    convertMax+= tiny;
}    

mxArray* CSCompress::getMexOutput() const
{
    mxArray *out= mxCreateStructMatrix(1,1, 0, 0);

    mxAddField(out, "convertMin");
    mxAddField(out, "convertFac");

    mxArray *mxtmp;
    mxtmp= mxCreateScalarDouble(convertMin);
    hippo_Assert(mxtmp, "couldn't allocate space for variable min");
    mxSetField(out, 0, "convertMin", mxtmp);
    mxtmp= mxCreateScalarDouble(convertFac);
    hippo_Assert(mxtmp, "couldn't allocate space for variable convertFac");
    mxSetField(out, 0, "convertFac", mxtmp);

    mxtmp= mxCreateNumericMatrix(nTotal, nInterval, OTYPE_ID, mxREAL);
    MX_AllocAssert(mxtmp);
    OTYPE *outptr= (OTYPE*) mxGetData(mxtmp);
    for (int i= 0; i < nInterval*nTotal; i++) {
        hippo_Assert((tParam[i]-convertMin)*convertFac <= maxOTYPE, 
                "param value outside range (shouldn't get here");
        outptr[i]= (OTYPE) ((tParam[i]-convertMin) * convertFac);
    }
    mxAddField(out, "param");
    mxSetField(out, 0, "param", mxtmp);

    mxtmp= mxCreateScalarDouble(outputInterval);
    hippo_Assert(mxtmp, "couldn't allocate space for variable outputInterval");
    mxAddField(out, "outputInterval");
    mxSetField(out, 0, "outputInterval", mxtmp);

    return out;
}

#undef OTYPE
#undef OTYPE_ID
