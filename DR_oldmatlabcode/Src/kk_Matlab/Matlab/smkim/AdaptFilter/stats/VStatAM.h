#ifndef VSTATAM_H
#define VSTATAM_H

/* $Id: VStatAM.h,v 1.2 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, Fri Sep  3 12:58:58 PDT 2004

   program description
*/

#include "VStat.h"
#include "StatLimits.h"
#include "../aux/defs.h"
#include "../aux/TData.h"
#include "../model/AdaptModel.h"

/* ************************************************************
                          class VStatAM
   ************************************************************ */

namespace AFilter {


class VStatAM : public VStat
{
private:
protected:
    const AdaptModel* model;
public:
    VStatAM(const TData*, const AdaptModel*, const StatLimits*);
    VStatAM(VStatAM &);
    virtual ~VStatAM();

    virtual void init_obj();

    virtual void run(const TData*, const AdaptModel*, const StatLimits*);

};


} // namespace AFilter

#endif   // VSTATAM_H
