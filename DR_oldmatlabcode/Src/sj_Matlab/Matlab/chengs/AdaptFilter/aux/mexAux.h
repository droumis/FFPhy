#ifndef MEXAUX_H
#define MEXAUX_H

/* $Id: mexAux.h,v 1.2 2006/11/07 18:02:14 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "mex.h"
#include <stdio.h>

#define MX_Warn(formatString,arg)					\
do{									\
    char tmpstring[500];						\
        snprintf(tmpstring, 500, formatString " in %s, line %d\n", arg,	\
                __FILE__, __LINE__);					\
        mexWarnMsgTxt(tmpstring);					\
} while(0);

#define MX_Error(formatString,arg)					\
do{									\
    char tmpstring[500];						\
        snprintf(tmpstring, 500, formatString " in %s, line %d\n", arg,	\
                __FILE__, __LINE__);					\
        mexErrMsgTxt(tmpstring);					\
} while(0);

#define MX_CAlloc(ptr,N,TYPE)						\
(ptr)= (TYPE *) mxCalloc((N), sizeof(TYPE));			\
    if ((ptr)== NULL) {							\
        char tmpstring[100];						\
            snprintf(tmpstring, 100, "memory allocation failed in %s, line %d\n"	\
                    , __FILE__, __LINE__);					\
            mexErrMsgTxt(tmpstring);					\
    } 

#define MX_CreateDouble(mxArrayPtr,M,N)					\
(mxArrayPtr)= mxCreateDoubleMatrix(M, N, mxREAL);			\
    if ((mxArrayPtr)== NULL) {						\
        char tmpstring[100];						\
            snprintf(tmpstring, 100, "memory allocation failed in %s, line %d\n"	\
                    , __FILE__, __LINE__);					\
            mexErrMsgTxt(tmpstring);					\
    } 

#define MX_CreateDoubleAssign(mxArrayPtr,M,N,dblptr)			\
     MX_CreateDouble(mxArrayPtr,M,N) \
    (dblptr)= mxGetPr(mxArrayPtr);


#define MX_AssignDouble(mxArrayPtr,fieldName,M,N, values)	\
    do {\
        double *_dblptr_;\
        mxArray *_mxtmp_;\
        MX_CreateDoubleAssign(_mxtmp_,(M),(N),_dblptr_); \
        memcpy(_dblptr_, values, (M)*(N)* sizeof(double));\
        mxSetField(mxArrayPtr, 0, fieldName, _mxtmp_);\
    } while(0)

#define MX_AllocAssert(ptr)						\
if ((ptr)== NULL) {							\
    char tmpstring[100];						\
        snprintf(tmpstring, 100, "memory allocation failed in %s, line %d\n"	\
                , __FILE__, __LINE__);					\
        mexErrMsgTxt(tmpstring);					\
}

#define MX_MemAssert(ptr)						\
if ((ptr)== NULL) {							\
    char tmpstring[100];						\
        snprintf(tmpstring, 100, "memory not assigned in %s, line %d\n"	\
                , __FILE__, __LINE__);					\
        mexErrMsgTxt(tmpstring);					\
}


#define MX_Assign(mxArrayPtr,ptr)					\
(ptr)= mxGetPr(mxArrayPtr);						\
    if ((ptr)== NULL) {							\
        char tmpstring[100];						\
            snprintf(tmpstring, 100, "memory not assigned in %s, line %d\n"	\
                    , __FILE__, __LINE__);					\
            mexErrMsgTxt(tmpstring);					\
    }

#define MX_StringAllocAssign(mxPtr,ptr)					\
do {								\
    int buflen= (mxGetM(mxPtr)*mxGetN(mxPtr)*sizeof(mxChar)) +1;	\
        MX_CAlloc(ptr,buflen,char);					\
        int res= mxGetString(mxPtr, ptr, buflen); \
        mxAssert(res== 0, "string assignment failed");\
        res= 0;     \
} while(0);

#define MX_GetField(mxIn, mxOut, field)					\
{									\
    mxOut= mxGetField(mxIn,0,field);			\
        MX_Assert(mxOut, "mxArray has no field named \""		\
                field "\"");						\
}

#define MX_FieldAssign(mxPtr, field, ptr, n)				\
{									\
    mxArray *_FieldAssignTmp= mxGetField(mxPtr,0,field);		\
        MX_Assert(_FieldAssignTmp, "mxArray has no field named \""	\
                field "\"");						\
        MX_Assign(_FieldAssignTmp, ptr);				\
        n= mxGetM(_FieldAssignTmp)*mxGetN(_FieldAssignTmp);		\
}

#define MX_FieldScalarDefault(mxPtr, field, x, defaultval, type)	\
{									\
    mxArray *_FieldAssignTmp= mxGetField(mxPtr,0,field);		\
        if (_FieldAssignTmp) {						\
            MX_Assert(mxGetN(_FieldAssignTmp)==1 &&			\
                    mxGetM(_FieldAssignTmp)==1,			\
#field "does not contain only one number");	\
                    x= (type) mxGetScalar(_FieldAssignTmp);			\
        } else x= defaultval;						\
}

#define MX_FieldScalar(mxPtr, field, x, type)				\
{									\
    mxArray *_FieldAssignTmp= mxGetField(mxPtr,0,field);		\
        MX_Assert(_FieldAssignTmp, "mxArray has no field named \""	\
                field "\"");						\
        x= (type) mxGetScalar(_FieldAssignTmp);			\
}

#define MX_GetNumber(mxPtr,x,type)				\
MX_Assert(mxGetN(mxPtr)==1 && mxGetM(mxPtr)==1, #mxPtr	\
        "does not only contain one number");		\
        x= (type) mxGetScalar(mxPtr);

#ifdef DEBUG
#define MX_Assert(cond, string)						\
if (!(cond)) {							\
    char tmpstring[100];						\
        snprintf(tmpstring, 100, \
                "assertion failed in %s, line %d\nerror: %s\n"		\
                , __FILE__, __LINE__, string);				\
        mexErrMsgTxt(tmpstring);					\
}
#define MX_Assert2(cond, format, arg)						\
if (!(cond)) {							\
    char tmpstring[100];						\
        snprintf(tmpstring, 100, format "in %s, line %d\n"		\
                , arg, __FILE__, __LINE__);				\
        mexErrMsgTxt(tmpstring);					\
}

#else // DEBUG not defined
#define MX_Assert(cond, string) 
#define MX_Assert2(cond, format, arg) 
#endif // DEBUG

#endif   // MEXAUX_H
