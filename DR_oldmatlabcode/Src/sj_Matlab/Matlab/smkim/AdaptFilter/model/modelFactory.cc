/* $Id: modelFactory.cc,v 1.9 2008/08/24 19:47:05 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "modelFactory.h"
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"

#include "../model/Pos_Isi.h"
#include "../model/PosPhase_Isi.h"
#include "../model/PosXY_Isi.h"

#include "../model/Gauss.h"
//#include "../model/GaussMix.h"
#include "../model/Const.h"
//#include "../model/Uniform.h"
#include "../model/CSplines.h"
#include "../model/CSplinesNorm.h"
#include "../model/CSplines2d.h"

/* ************************************************************
   ************************************************************ */

using namespace AFilter;

AdaptModel* AFilter::allocModel(const mxArray *opts, TData *data, bool warn)
{
    AdaptModel *model;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "model.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     

    if (strcmp(name,"PosPhase_Isi")==0) {
        model= new AFilter::PosPhase_Isi(data);
    } else if (strcmp(name,"Pos_Isi")==0) {
        model= new AFilter::Pos_Isi(data);
    } else if (strcmp(name,"PosXY_Isi")==0) {
        model= new AFilter::PosXY_Isi(data);
    } else if (strcmp(name,"Const")==0) {
        model= new AFilter::Const(data);
    } else {
        if (warn) {
            MX_Error("unknown model: '%s' ", name);
            model= 0;
        } else
            return 0;
    }
    model->setMexInput(opts);
    model->init();

    mxFree(name);
    return model;
}

VDynamics1d* AFilter::allocDynamics1d(const mxArray *opts, TData *data, bool warn)
{
    VDynamics1d *dynamics;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "dynamics.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     


    if (strcmp(name,"Const")==0) {
        dynamics= new AFilter::Const();
//        hippo_Print(dynamics);
//    } else if (strcmp(name,"Uniform")==0) {
//        dynamics= new AFilter::Uniform();
    } else if (strcmp(name,"CSplines")==0) {
        dynamics= new AFilter::CSplines();
    } else if (strcmp(name,"CSplinesNorm")==0) {
        dynamics= new AFilter::CSplinesNorm();
    } else {
        if (warn) {
            MX_Error("unknown 1-d dynamics model: %s ", name);
            dynamics= 0;
        } else 
            return 0;
    }
    dynamics->setData(data);
    dynamics->setMexInput(opts);

//    hippo_Print(dynamics);
    mxFree(name);
    return dynamics;
}

VDynamics2d* AFilter::allocDynamics2d(const mxArray *opts, TData *data, bool warn)
{
    VDynamics2d *dynamics;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "dynamics.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     


    if  (strcmp(name,"Gauss")==0) {
        dynamics= new AFilter::Gauss();
//     } else if  (strcmp(name,"GaussMix")==0) {
// 	dynamics= new AFilter::GaussMix();
    } else if (strcmp(name,"Const")==0) {
        dynamics= new AFilter::Const();
//    } else if (strcmp(name,"Uniform")==0) {
//        dynamics= new AFilter::Uniform();
    } else if (strcmp(name,"CSplines2d")==0) {
        dynamics= new AFilter::CSplines2d();
    } else {
        if (warn) {
            MX_Error("unknown 2-d dynamics model: %s ", name);
            dynamics= 0;
        } else
            return 0;
    }
    dynamics->setData(data);
    dynamics->setMexInput(opts);

//    hippo_Print(dynamics);
//    hippo_Print(dynamics->getMinX());

    mxFree(name);
    return dynamics;
}

VDynamics* AFilter::allocDynamics(const mxArray *opts, TData *data, bool warn)
{
    VDynamics *dyn= allocDynamics1d(opts, data, false);

    if (!dyn) dyn= allocDynamics2d(opts, data, false);

    if (!dyn) {
        if (warn) {
            char *name;
            mxArray *mxTmp= mxGetField(opts,0,"name");
            hippo_Assert(mxTmp, "dynamics.name was not defined");
            MX_StringAllocAssign(mxTmp,name);     
            MX_Error("unknown dynamics model: %s ", name);
            dyn= 0;
        } else return 0;
    }

    return dyn;
}
