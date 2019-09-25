// $Id: mexMain.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include <FL/Fl.H>
#include "afViewUI.h"
#include <iostream>

#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"
#include "../aux/TData.h"
#include "../model/modelFactory.h"

using namespace AFilter;

AFViewUI *avui;

// usage: vis(data, isimodel, xpmodel, visopts)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //  Check numbers of arguments 
    if (  (!( nlhs <= 1)) || ((nrhs != 3))) {
        mexErrMsgTxt("incorrect number of parameters, type \"help vis\" for usage (not implemented)\n");
    }

    //  Creating data objects
    TData *data = new TData(prhs[0]);
    AdaptModel *am= allocModel(prhs[1], data);
    PosPhase_Isi *model= static_cast<PosPhase_Isi*>(am);

    AFDataModel *dataModel = new AFDataModel(data,model);
    dataModel->setColormap(prhs[2]);

    // Creating GUI
    avui = new AFViewUI();

    avui->setDataModel(dataModel);

    Fl::visual(FL_DOUBLE|FL_INDEX);

    // Showing window
    avui->show();

    Fl::run();
    Fl::flush();

    delete avui;
    delete am;
    delete data;
}
