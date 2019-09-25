// $Id: visModel.cc,v 1.2 2008/10/23 21:02:33 chengs Exp $
#include <FL/Fl.H>
#include <iostream>

#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"
#include "../aux/TData.h"
#include "../model/modelFactory.h"
//#include "VarNIsiVis.h"
//#include "TwoVarNIsiVis.h"
#include "GUIMain.h"
#include "ModelVis.h"

using namespace AFilter;

GUIMain *gui;

// usage: visModel(data, model, colormap)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //  Check numbers of arguments 
    if (  (!( nlhs <= 1)) || ((nrhs != 3))) {
        mexErrMsgTxt("incorrect number of parameters, type \"help visModel\" for usage (not implemented)\n");
    }

    cerr << "\nVisualizing model ...\n";

    //  Creating data objects
    TData *data = new TData(prhs[0]);
    char *name;
    mxArray *mxTmp= mxGetField(prhs[1],0,"name");
    hippo_Assert(mxTmp, "model.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     

    ModelVis *mv= new ModelVis(data);
//    if (strcmp(name,"PosPhase_Isi")==0) {
//        mv= new TwoVarNIsiVis(data);
//    } else if (strcmp(name,"Pos_Isi")==0) {
//        mv= new VarNIsiVis(data);
//    } else  {
//        MX_Error("unknown model: '%s' ", name);
//    }
    mv->setMexInput(prhs[1],prhs[2]);

    gui= new GUIMain;
    gui->setModel(mv);
    gui->setMexInput(prhs[2]);
    Fl::visual(FL_DOUBLE|FL_INDEX);
    gui->show();
    Fl::run();
    Fl::flush();

    delete gui->mainWindow; // actually, fltk should take care of this, but it dosn't 
    delete gui;
    delete mv;
    delete data;
//    hippo_Mark;
}
