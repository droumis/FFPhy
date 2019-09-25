#ifndef PosXY_Isi_H
#define PosXY_Isi_H

/* $Id: PosXY_Isi.h,v 1.3 2008/08/24 19:47:04 chengs Exp $
   Sen Cheng, Thu Aug 21 12:15:58 PDT 2008
*/

#include "../model/IsiModel.h"

/* ************************************************************
                          class PosXY_Isi
   ************************************************************ */
namespace AFilter {

class PosXY_Isi : public IsiModel {
protected:
    virtual AdaptModel* alloc() const { return new PosXY_Isi; }
    virtual void init_spatial(int n) {
        VDynamics2d * dyn2= get2dSpatialModel();
        dyn2->init(startindex, endindex, n, 
                data->getXPtr(), data->getYPtr(), data->ntraj, data->fieldID);
    }
public:
    PosXY_Isi(TData *_data=0)  : IsiModel(_data){ };
    virtual ~PosXY_Isi() { };
    virtual const char* getName() const { return "PosXY_Isi"; };
};


} // namespace AFilter

#endif   // PosXY_Isi_H
