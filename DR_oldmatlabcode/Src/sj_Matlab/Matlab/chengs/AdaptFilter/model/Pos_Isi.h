#ifndef POS_ISI_H
#define POS_ISI_H

/* $Id: Pos_Isi.h,v 1.10 2008/08/24 19:47:04 chengs Exp $
   
*/

#include "../model/IsiModel.h"

/* ************************************************************
                          class Pos_Isi
   ************************************************************ */
namespace AFilter {

class Pos_Isi : public IsiModel {
protected:
    virtual AdaptModel* alloc() const { return new Pos_Isi; }
    virtual void init_spatial(int n) {
        VDynamics1d * dyn1= get1dSpatialModel();
        dyn1->init(startindex, endindex, n, 
                data->posarray, data->ntraj, data->fieldID);
    }
public:
    Pos_Isi(TData *_data=0) : IsiModel(_data) { };
    virtual ~Pos_Isi() { };
    virtual const char* getName() const { return "Pos_Isi"; };
};


} // namespace AFilter

#endif   // POS_ISI_H
