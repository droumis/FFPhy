#ifndef POSPHASE_ISI_H
#define POSPHASE_ISI_H

/* $Id: PosPhase_Isi.h,v 1.8 2008/08/24 19:47:04 chengs Exp $
   
   Sen Cheng, Thu Aug 21 12:17:48 PDT 2008
*/

#include "../model/IsiModel.h"

/* ************************************************************
                          class PosPhase_Isi
   ************************************************************ */
namespace AFilter {

class PosPhase_Isi : public IsiModel {
protected:
    inline virtual AdaptModel* alloc() const { return new PosPhase_Isi; }

    virtual void init_spatial(int n) {
        VDynamics2d * dyn2= get2dSpatialModel();

        dyn2->init(startindex, endindex, n, 
                data->posarray, data->phasearray, data->ntraj, data->fieldID);
    }

public:

    PosPhase_Isi(TData *_data=0) : IsiModel(_data) { };

    virtual ~PosPhase_Isi() { };

    virtual const char* getName() const { return "PosPhase_Isi"; };
};


} // namespace AFilter

#endif   // POSPHASE_ISI_H
